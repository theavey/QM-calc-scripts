#! /usr/bin/env python3
"""
Run QM calculations at multiple levels consecutively, using queuing system

Automate Multi-Level calculations

This should help with running a set of QM calculations at several levels (e.g.,
increasing basis set size), while intelligently using the queuing system such
as Sun Grid Engine.
It can receive the signal from the queuing system that the job will be killed
soon and consequently submit a continuation of the job, using the
intermediate files to speed up subsequent calculations.
"""

########################################################################
#                                                                      #
# This script/module was written by Thomas Heavey in 2018.             #
#        theavey@bu.edu     thomasjheavey@gmail.com                    #
#                                                                      #
# Copyright 2018 Thomas J. Heavey IV                                   #
#                                                                      #
# Licensed under the Apache License, Version 2.0 (the "License");      #
# you may not use this file except in compliance with the License.     #
# You may obtain a copy of the License at                              #
#                                                                      #
#    http://www.apache.org/licenses/LICENSE-2.0                        #
#                                                                      #
# Unless required by applicable law or agreed to in writing, software  #
# distributed under the License is distributed on an "AS IS" BASIS,    #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or      #
# implied.                                                             #
# See the License for the specific language governing permissions and  #
# limitations under the License.                                       #
#                                                                      #
########################################################################

import json
import logging
import MDAnalysis as mda
import numpy as np
import os
import paratemp
from paratemp.geometries import XYZ
import pathlib
import re
import shutil
import signal
import subprocess
import sys
import threading
import thtools
import time
from . import tools

if not sys.version_info >= (3, 6):
    raise ValueError('Python >= 3.6 is required')


class Calc(object):

    def __init__(self, status=None, base_name=None, ind=None, top=None,
                 traj=None, criteria=None, react_dist=None, ugt_dicts=None):
        """

        :param str status: The path to the status file to be read for a
            calculation restart. If this is not a restarted job, this should
            be None (the default).
        :param str base_name:
        :param int ind:
        :type top: pathlib.Path or str
        :param top:
        :type traj: pathlib.Path or str
        :param traj:
        :param dict criteria: The criteria for selecting frames from the
            trajectory.
            This is a dict with distance names (or other columns that will
            be in `Universe.data`) as the keys and the values being a
            List-like of min and max values.
            For example, `{'c1_c2': (1.5, 4.0), 'c1_c3': (2.2, 5.1)}` will
            select frames where 'c1_c2' is between 1.5 and 4.0 and 'c1_c3'
            is between 2.2 and 5.1.
        :type react_dist: str or float
        :param react_dist: Distance to set between the two reacting atoms (
            with indices 20 and 39).
            If this argument as given evaluates to False, no movement/changes
            to the geometry will be made.
        :param List[dict] ugt_dicts:
        :return:
        """
        if status is not None:
            self._status = StatusDict(status)
            self._json_name = status
            self.args = self.status['args']
            a = self.args
            base_name = a['base_name']
            ind = a['ind']
            top = a['top']
            traj = a['traj']
            criteria = a['criteria']
            react_dist = a['react_dist']
            ugt_dicts = a['ugt_dicts']
            self._base_name = self.status['base_name']
        else:
            self.args = {
                'base_name': base_name,
                'ind': ind, 'top': top, 'traj': traj,
                'criteria': criteria,
                'react_dist': react_dist,
                'ugt_dicts': ugt_dicts}
            self.check_args()
            self._base_name = '{}-ind{}'.format(base_name, ind)
            self._json_name = '{}.json'.format(self._base_name)
            self._status = StatusDict(self._json_name)
        self.top = top
        self.traj = traj
        self.criteria = criteria
        self.react_dist = react_dist
        self.ugt_dicts = ugt_dicts
        self.log = logging.getLogger('{}.log'.format(self._base_name))
        self.mem, self.node = None, None
        self.scratch_path: pathlib.Path = None
        self.last_scratch_path: pathlib.Path = None
        self.n_slots, self.last_node = None, None
        self.cwd_path: pathlib.Path = None
        self.output_scratch_path: pathlib.Path = None

    @property
    def status(self):
        return self._status

    def check_args(self):
        for key in self.args:
            if self.args[key] is None:
                raise ValueError(f'Argument "{key}" cannot be None')

    def _startup_tasks(self):
        node = os.environ['HOSTNAME'].split('.')[0]
        self.node = node
        scratch_path = pathlib.Path('/net/{}/scratch/theavey'.format(node))
        scratch_path.mkdir(exist_ok=True)
        self.scratch_path = scratch_path
        self.mem = thtools.job_tools.get_node_mem()
        n_slots = int(os.environ['NSLOTS'])
        self.n_slots = n_slots
        self.log.info('Running on {} using {} cores and up to {} GB '
                      'mem'.format(node, n_slots, self.mem))
        if self.status:
            self.last_node = self.status['current_node']
            self.status['last_node'] = self.last_node
            node_list = self.status['node_list']
            self.last_scratch_path = pathlib.Path(self.status[
                                                      'current_scratch_dir'])
            self.status['last_scratch_dir'] = str(self.last_scratch_path)
        else:
            self.status['args'] = self.args
            self.status['base_name'] = self._base_name
            node_list = []
        self.status['node_list'] = node_list + [node]
        self.status['current_node'] = node
        self.status['current_scratch_dir'] = str(scratch_path)
        self.cwd_path = pathlib.Path('.').resolve()
        self.status['cwd'] = str(self.cwd_path)
        self.log.info('Submitted from {} and will be running in {}'.format(
            self.cwd_path, self.scratch_path))

    def run_calc(self):
        """


        :return:
        """
        rerun = True if self.status else False
        self._startup_tasks()
        if rerun:
            self.log.info('loaded previous status file: {}'.format(
                self._json_name))
            self.resume_calc()
        else:
            self.log.warning('No previous status file found. '
                             'Starting new calculation?')
            self.new_calc()

    def _make_rand_xyz(self):
        u = paratemp.Universe(self.top, self.traj)
        u.read_data()
        frames = u.select_frames(self.criteria, 'QM_frames')
        select = np.random.choice(frames)
        self.status['source_frame_num'] = select
        system = u.select_atoms('all')
        xyz_name = self._base_name + '.xyz'
        with mda.Writer(xyz_name, system.n_atoms) as w:
            u.trajectory[select]
            w.write(system)
        self.log.info('Wrote xyz file from frame {} to {}'.format(select,
                                                                  xyz_name))
        if self.react_dist:
            paratemp.copy_no_overwrite(xyz_name, xyz_name+'.bak')
            self.log.info('Copied original geometry to {}'.format(
                xyz_name+'.bak'))
            xyz = XYZ(xyz_name)
            diff = xyz.coords[20] - xyz.coords[39]
            direction = diff / np.linalg.norm(diff)
            xyz.coords[20] = xyz.coords[39] + self.react_dist * direction
            xyz.write(xyz_name)
            self.log.info('Moved reactant atoms 20 and 39 to a distance of {} '
                          'and wrote this to {}'.format(self.react_dist,
                                                        xyz_name))
            self.status['original_xyz'] = xyz_name + '.bak'
        self.status['starting_xyz'] = xyz_name
        return pathlib.Path(xyz_name).resolve()

    def new_calc(self):
        xyz_path = self._make_rand_xyz()
        bn = self._base_name
        com_name = f'{bn}-lvl0.com'
        tools.use_gen_template(
            out_file=com_name,
            xyz=str(xyz_path),
            job_name=bn,
            checkpoint=f'{bn}.chk',
            rwf=f'{bn}.rwf',
            nproc=self.n_slots, mem=self.mem,
            **self.ugt_dicts[0]
        )
        self.log.info('Wrote Gaussian input for '
                      f'first job to {com_name}')
        self.status['g_in_0'] = com_name
        chk_ln_path = pathlib.Path(f'{bn}-running.chk')
        chk_ln_path.symlink_to(self.scratch_path.joinpath(f'{bn}.chk'))
        self.log.info(f'Linked checkpoint file as {chk_ln_path.resolve()}')
        killed = self._run_gaussian(com_name)
        if killed:
            self.status['calc_cutoff'] = True
            self.resub_calc()
            self.log.info('Resubmitted. Exiting this job')
            chk_ln_path.unlink()
            sys.exit('Exiting because job timeout')  # probably don't do this
        else:
            self.status['calc_cutoff'] = False
            self._check_normal_completion(self.output_scratch_path)
        self._copy_back_files(com_name, killed)

    def _run_gaussian(self, com_name):
        out_name = com_name.replace('com', 'out')
        com_path: pathlib.Path = self.cwd_path.joinpath(com_name)
        if not com_path.exists():
            raise FileNotFoundError('Gaussian input {} not found in '
                                    '{}'.format(com_name, self.cwd_path))
        out_path: pathlib.Path = self.scratch_path.joinpath(out_name)
        self.output_scratch_path = out_path
        signal.signal(signal.SIGUSR2, self._signal_catch_time)
        signal.signal(signal.SIGUSR1, self._signal_catch_done)
        cl = ['g16', ]
        killed = False
        with com_path.open('r') as f_in, out_path.open('w') as f_out:
            self.log.info('Starting Gaussian with input {} and writing '
                          'output to {}'.format(com_path, out_path))
            proc = subprocess.Popen(cl, stdin=f_in, stdout=f_out,
                                    cwd=self.scratch_path)
            self.log.info('Started Gaussian; waiting for it to finish or '
                          'timeout')
            try:
                threading.Thread(target=self._check_proc, args=(proc,))
                signal.pause()
            except self.TimesUp:
                killed = True
                proc.terminate()  # Should be within `with` clause?
                self.log.info('Gaussian process terminated because of SIGUSR2')
            except self.GaussianDone:
                self.log.info('Gaussian process completed')
        return killed

    def _signal_catch_time(self, signum, frame):
        self.log.warning('Caught SIGUSR2 signal! Trying to quit Gaussian '
                         'and resubmit continuation calculation')
        raise self.TimesUp

    def _signal_catch_done(self, signum, frame):
        self.log.warning('/ Caught SIGUSR1 signal! Likely, this was because '
                         'Gaussian process exited')
        raise self.GaussianDone

    def _check_proc(self, proc):
        while proc.poll() is None:
            time.sleep(15)
        self.log.warning('Gaussian process completed. Sending SIGUSR1')
        os.kill(os.getpid(), signal.SIGUSR1)

    def _copy_back_files(self, com_name: str, killed: bool):
        if not killed:
            out_path = pathlib.Path(com_name.replace('com', 'out'))
        else:
            outs = list(self.cwd_path.glob(com_name[:-4]+'-*.out'))
            outs.sort()
            outs.sort(key=len)
            new_out = re.sub(r'(\d+)\.out',
                             lambda m: f'{int(m.group(1))+1)}.out',
                             outs[-1])
            out_path = pathlib.Path(new_out)
        paratemp.copy_no_overwrite(self.output_scratch_path, out_path)
        chk_name = f'{self._base_name}.chk'
        shutil.copy(self.output_scratch_path.joinpath(chk_name), chk_name)

    def _check_normal_completion(self, filepath):
        pass

    def resub_calc(self):
        pass

    def resume_calc(self):
        # TODO rwf/chk needs to be copied over
        pass

    class TimesUp(Exception):
        pass

    class GaussianDone(Exception):
        pass


class StatusDict(dict):
    """
    A dict subclass that writes the dict to disk every time a value gets set

    Note, any other action on the dict will not currently trigger a write to
    disk.

    The dict will be written in JSON to the path given at instantiation. If
    there is already a file at that path, it will be read and used as
    the initial definition of the dict. Otherwise, it will be instantiated as
    an empty dict.

    Keys in dictionary:

    * args: a Dict of the arguments given when starting the calculation.

    * current_node: Name of the node on which the job is currently running.
    This is set during :func:`Calc._startup_tasks`. It should be formatted as
    'scc-xxx' (e.g., 'scc-na1').

    * last_node: This is set during :func:`Calc._startup_tasks` if it's a
    restarted calculation. It will be set from the last 'current_node'.

    * node_list: This is a list of past and current nodes on which this
    calculation has run.

    * current_scratch_dir: str of the absolute path to the scratch directory
    on the current node.

    * base_name: the base name for this calculation including the index of
    this calculation.

    * cwd: str of the absolute path from which the current calculation was
    submitted.

    * last_scratch_dir: str of the absolute path to the scratch directory
    from which the last job was run, if this is not a new calculation.

    * source_frame_num: The index of the frame in the trajectory that was
    used to create the initial configuration.

    * original_xyz: str of the name of file with the coordinates as they were
    taken from the trajectory, before moving the reacting atoms to the
    correct distance. This will not be set if no distance correction is made.

    * starting_xyz: str of the name of the file with the coordinates for
    starting the calculation, before any optimization.

    * g_in_0: str of the name of the file with the initial input to Gaussian.
    * calc_cutoff: bool of whether the job finished or if it was cutoff
    because of running out of time.

    """
    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()
        if pathlib.Path(path).is_file():
            d = json.load(open(path, 'r'))
            super(StatusDict, self).__init__(d)
        else:
            super(StatusDict, self).__init__()

    def __setitem__(self, key, value):
        super(StatusDict, self).__setitem__(key, value)
        json.dump(self, open(self.path, 'w'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--base_name', type=str, required=True,
                        help='base name for this calculation, likely not '
                             'including any index')
    parser.add_argument('-i', '--index', type=int,
                        help='index of this calculation')
    parser.add_argument('-c', '--top', type=str,
                        help='topology/structure file (e.g., .gro, .xyz)')
    parser.add_argument('-f', '--trajectory', type=str,
                        help='trajectory file (e.g., .xtc, .trr, .dcd)')

    def parse_crit(kvv):
        k, vv = kvv.split('=')
        vs = (int(v) for v in vv.split(','))
        return k, vs

    parser.add_argument('-s', '--criteria', action='append',
                        type=parse_crit, metavar='key=min,max',
                        help='criteria for selection of possible frames from '
                             'the trajectory. To provide more than one '
                             'criterion, use this argument multiple times')
    parser.add_argument('-d', '--react_dist', type=float, default=False,
                        help='Distance to set between atoms 20 and 39, '
                             'in angstroms. If this evaluates to False, '
                             'no changes to the geometry will be made')
    parser.add_argument('-g', '--ugt_dicts', type=str, required=True,
                        help='path to json file that parses to a list of '
                             'dicts of arguments for use_gen_template in '
                             'order to create inputs to Gaussian')
    parser.add_argument('--restart', default=None,
                        help='Path to status file for resuming an already '
                             'started calculation')
    args = parser.parse_args()
    if args.restart is not None:
        calc = Calc(status=args.restart)
    else:
        calc = Calc(base_name=args.base_name,
                    ind=args.index,
                    top=args.top,
                    traj=args.trajectory,
                    criteria=dict(args.criteria),
                    react_dist=args.react_dist,
                    ugt_dicts=json.load(open(args.ugt_dicts, 'r'))
                    )
    calc.run_calc()
