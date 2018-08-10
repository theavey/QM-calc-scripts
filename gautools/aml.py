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
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
import signal
import subprocess
import sys
import threading
import thtools
import time
from . import tools


class Calc(object):

    def __init__(self, base_name, ind, top, traj, criteria, ugt_dicts):
        """

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
        :param List[dict] ugt_dicts:
        :return:
        """
        # TODO: get ind (index) from SGE_TASK_ID?
        # do this if not given, and save it to args for passing on?
        self.args = base_name, ind, top, traj, criteria, ugt_dicts
        self.top = top
        self.traj = traj
        self.criteria = criteria
        self.ugt_dicts = ugt_dicts
        # TODO: save args to disk?
        # Could then give these defaults, and if it's a continuation, just load
        # all the arguments from disk. Need to think about how this will
        # be used, exactly
        self._base_name = '{}-ind{}'.format(base_name, ind)
        self.log = logging.getLogger('{}.log'.format(self._base_name))
        self._json_name = '{}.json'.format(self._base_name)
        self._status = StatusDict(self._json_name)
        self.mem, self.node, self.scratch_path = None, None, None
        self.last_scratch_path, self.n_slots, self.last_node = None, None, None
        self.cwd_path = None

    @property
    def status(self):
        return self._status

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
        self._startup_tasks()
        if self.status:
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
        self.status['original_xyz'] = xyz_name

    def new_calc(self):
        self._make_rand_xyz()
        bn = self._base_name
        com_name = bn + '-lvl0.com'
        tools.use_gen_template(
            out_file=com_name,
            xyz=bn + '.xyz',
            job_name=bn,
            checkpoint=bn + '.chk',
            **self.ugt_dicts[0]
        )
        self.log.info('Wrote Gaussian input for '
                      'first job to {}'.format(com_name))
        self.status['g_in_0'] = com_name
        self._run_gaussian(com_name)

    def _run_gaussian(self, com_name):
        # TODO rwf/chk needs to be copied over
        out_name = com_name.replace('com', 'out')
        com_path: pathlib.Path = self.cwd_path.joinpath(com_name)
        if not com_path.exists():
            raise FileNotFoundError('Gaussian input {} not found in '
                                    '{}'.format(com_name, self.cwd_path))
        out_path: pathlib.Path = self.scratch_path.joinpath(out_name)
        signal.signal(signal.SIGUSR2, self._signal_catch_time)
        signal.signal(signal.SIGUSR1, self._signal_catch_done)
        cl = ['g16',
              '-m="{}GB"'.format(self.mem),
              '-c="0-{}"'.format(self.n_slots-1), ]
        killed = False
        with com_path.open('r') as f_in, out_path.open('w') as f_out:
            self.log.info('Starting Gaussian with input {} and writing '
                          'output to {}'.format(com_name, out_name))  # TODO fix
            # TODO link checkpoint and such (maybe not in this function)
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
        if killed:
            self.status['calc_cutoff'] = True
            self.resub_calc()
            self.log.info('Resubmitted. Exiting this job')
            # TODO copy back file(s) (maybe not here?)
            # TODO remove links to files (maybe not here?)
            sys.exit('Exiting because job timeout')  # probably don't do this
        else:
            self.status['calc_cutoff'] = False
            # TODO copy back files (maybe not here?)
            self._check_normal_completion(out_path)
        pass

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

    def _check_normal_completion(self, filepath):
        pass

    def resub_calc(self):
        pass

    def resume_calc(self):
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
