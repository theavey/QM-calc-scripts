#! /usr3/graduate/theavey/anaconda_envs/py3.6/bin/python
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
from gautools import tools
import functools
import filecmp

if not sys.version_info >= (3, 6):
    raise ValueError('Python >= 3.6 is required')

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.WARNING)
formatter = logging.Formatter('%(asctime)s - %(name)s - '
                              '%(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)


obabel_module_lines = ('\n'
                       'module load wxwidgets/3.0.2\n'
                       'module load openbabel/2.4.1\n')

signal_catch_lines = ('PY_PID=$!\n'
                      'trap "kill -n 12 $PY_PID" usr2\n'
                      'wait\n\n')


def log_exception(f):
    @functools.wraps(f)
    def log_exc(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception:
            log.exception(f'An exception was raised in {f.__name__}!')
            raise

    return log_exc


class Calc(object):
    """
    An object for consecutive Gaussian calculation using a queuing system

    Environment variables that need to be defined:
        * HOSTNAME
        * JOB_ID
        * NSLOTS
        * SGE_STDOUT_PATH
    An exception will be raised if any of these are not defined when doing
    :func:`Calc.run_calc()`.
    """

    def __init__(self, status=None, base_name=None, ind=None,
                 geometry=None, top=None, traj=None,
                 criteria=None, react_dist=None, ugt_dicts=None):
        """

        :param str status: The path to the status file to be read for a
            calculation restart. If this is not a restarted job, this should
            be None (the default).
        :param str base_name:
        :param int ind:
        :type geometry: pathlib.Path or str
        :param geometry: File from which to get the starting coordinates. This
            argument will take priority if top and traj are also given.
            Currently, this must be an XYZ file, but it shouldn't be too hard
            to implement using openbabel to convert to an XYZ.
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
            # TODO this could be done (better?) with a classmethod
            self.rerun = True
            self._status = StatusDict(status)
            self._json_name = status
            try:
                self.args = self.status['args']
            except KeyError:
                self.log.error('Attempting to do a calculation restart, '
                               f'but it seems that the status file {status} '
                               'was empty or at least did not have "args" in '
                               'it.')
                raise
            a = self.args
            base_name = a['base_name']
            ind = a['ind']
            try:
                geometry = a['geometry']
            except KeyError:
                geometry = None
            top = a['top']
            traj = a['traj']
            criteria = a['criteria']
            react_dist = a['react_dist']
            ugt_dicts = a['ugt_dicts']
            self._base_name = self.status['base_name']
        else:
            self.rerun = False
            self.args = {
                'base_name': base_name,
                'ind': ind, 'geometry': geometry,
                'top': top, 'traj': traj,
                'criteria': criteria,
                'react_dist': react_dist,
                'ugt_dicts': ugt_dicts}
            self.check_args()
            self._base_name = '{}-ind{}'.format(base_name, ind)
            self._json_name = '{}.json'.format(self._base_name)
            self._status = StatusDict(self._json_name)
        self.geometry = geometry
        self.top = top
        self.traj = traj
        self.criteria = criteria
        self.react_dist = react_dist
        self.ugt_dicts = ugt_dicts
        self.log = logging.getLogger(self.__class__.__name__)
        self.log.setLevel(logging.DEBUG)
        f_handler = logging.FileHandler('{}.log'.format(self._base_name))
        f_handler.setLevel(logging.DEBUG)
        f_formatter = logging.Formatter('%(asctime)s - %(name)s - '
                                        '%(levelname)s - %(message)s')
        f_handler.setFormatter(f_formatter)
        self.log.addHandler(f_handler)
        self.log.addHandler(handler)
        self.log.debug('Initializing the log')
        self.mem, self.node = None, None
        self.scratch_path: pathlib.Path = None
        self.last_scratch_path: pathlib.Path = None
        self.n_slots, self.last_node = None, None
        self.cwd_path: pathlib.Path = None
        self.output_scratch_path: pathlib.Path = None
        self.chk_ln_path: pathlib.Path = None
        self.h_rt: str = None
        self.stdout_file: str = None
        self.resub_cl: list[str] = None
        self.job_id: str = None
        self.next_job_id: str = None
        self.resubmitted: bool = False

    def check_args(self):
        args = self.args.copy()
        geom = args.pop('geometry')
        top, traj = args.pop('top'), args.pop('traj')
        crit = args.pop('criteria')
        if (geom is None and
                (top is None or traj is None or crit is None)):
            raise ValueError('either geometry or top, traj, and criteria must '
                             'be given')
        for key in self.args:
            if self.args[key] is None:
                raise ValueError(f'Argument "{key}" cannot be None')

    @property
    def status(self):
        return self._status

    @property
    def current_lvl(self):
        try:
            return self.status['current_lvl']
        except KeyError:
            raise AttributeError('Could not find current level')

    @current_lvl.setter
    def current_lvl(self, value):
        if not (isinstance(value, int) or value is None):
            raise TypeError('current_level must be an int or None.'
                            f'Given type {type(value)}')
        self.status['current_lvl'] = value

    @property
    def between_levels(self):
        # Might need to look for rwf file otherwise should start the
        # calculation again
        self.log.debug('Checking to see if calculation left off between '
                       'calculation levels')
        try:
            between = self.status['between_levels']
        except KeyError:
            self.log.warning('No key in status for determining if between '
                             'levels currently')
            lvl = self.current_lvl
            try:
                if self.status[f'g_in_{lvl}'] == self.status['g_in_curr']:
                    out_path = pathlib.Path(
                        self.status['g_in_curr']).with_suffix('.out')
                    if out_path.exists():
                        self._advance_level()
                        between = True
                    else:
                        between = False
                else:
                    between = True
            except KeyError:
                between = True
            self.status['between_levels'] = between
        return between

    @between_levels.setter
    def between_levels(self, value):
        if not isinstance(value, bool):
            raise TypeError(f'between_levels must be a bool, given type '
                            f'{type(value)}')
        self.status['between_levels'] = value

    @property
    def cleaned_up(self):
        try:
            cleaned = self.status['cleaned_up']
        except KeyError:
            self.log.warning('Could not find "cleaned_up" in status. Assuming '
                             'dirty')
            cleaned = False
            self.status['cleaned_up'] = cleaned
        return cleaned

    @cleaned_up.setter
    def cleaned_up(self, value):
        if not isinstance(value, bool):
            raise TypeError(f'cleaned_up must be a bool, given type '
                            f'{type(value)}')
        self.status['cleaned_up'] = value

    def _startup_tasks(self):
        """
        Some startup tasks to set variables for later use

        This requires the environment variables HOSTNAME and NSLOTS be set.
        :return: None
        """
        self.log.debug('Running some introductory tasks and setting variables')
        try:
            node = os.environ['HOSTNAME'].split('.')[0]
        except KeyError:
            self.log.exception('Could not find HOSTNAME!')
            raise
        self.node = node
        scratch_path = pathlib.Path('/net/{}/scratch/theavey'.format(node))
        scratch_path.mkdir(exist_ok=True)
        self.scratch_path = scratch_path
        try:
            self.job_id = os.environ['JOB_ID']
        except KeyError:
            self.log.exception('Could not find JOB_ID!')
            raise
        self.mem = thtools.job_tools.get_node_mem()
        try:
            n_slots = int(os.environ['NSLOTS'])
        except KeyError:
            self.log.exception('Could not find NSLOTS!')
            raise
        self.n_slots = n_slots
        self.h_rt = self._get_h_rt()
        self.log.info(f'Running on {node} using {n_slots} cores and up to '
                      f'{self.mem} GB mem for {self.h_rt} seconds')
        try:
            self.stdout_file = os.environ['SGE_STDOUT_PATH']
            self.log.debug(f'Using stdout path: {self.stdout_file}')
        except KeyError:
            self.log.exception('Could not find SGE_STDOUT_PATH!')
            raise
        if self.rerun:
            self.last_node = self.status['current_node']
            self.status['last_node'] = self.last_node
            node_list = self.status['node_list']
            self.last_scratch_path = pathlib.Path(self.status[
                                                      'current_scratch_dir'])
            self.status['last_scratch_dir'] = str(self.last_scratch_path)
            self._get_chk_ln_path()
            self._get_output_scratch_path()
        else:
            self.status['args'] = self.args
            self.status['base_name'] = self._base_name
            node_list = []
        self.status['node_list'] = node_list + [node]
        self.status['current_node'] = node
        self.status['current_scratch_dir'] = str(scratch_path)
        self.status['job_id'] = self.job_id
        self.cwd_path = pathlib.Path('.').resolve()
        self.status['cwd'] = str(self.cwd_path)
        self.log.info('Submitted from {} and will be running in {}'.format(
            self.cwd_path, self.scratch_path))
        self._make_resub_sh_and_cl()

    def _get_output_scratch_path(self):
        self.log.debug('Getting path to scratch output')
        try:
            self.output_scratch_path = pathlib.Path(
                self.status['output_scratch_path'])
        except KeyError:
            self.output_scratch_path = self.last_scratch_path.joinpath(
                self.status['g_in_curr']).with_suffix('.out')

    def _get_chk_ln_path(self):
        self.log.debug('Getting path to linked chk file')
        try:
            self.chk_ln_path = pathlib.Path(self.status['chk_ln_path'])
        except KeyError:
            self.chk_ln_path = pathlib.Path(
                f'{self._base_name}-running.chk').resolve()

    @log_exception
    def run_calc(self):
        """
        The primary function to start (or restart) running a calculation

        :return: None
        """
        self.log.debug('Welcome. Just starting to run this calculation')
        self._startup_tasks()
        if self.rerun:
            self.log.info('loaded previous status file: {}'.format(
                self._json_name))
            self.resume_calc()
        else:
            self.log.warning('No previous status file found. '
                             'Starting new calculation?')
            self.new_calc()

    def _make_rand_xyz(self):
        self.log.debug('Making XYZ file to start calculation')
        import tables
        u = paratemp.Universe(self.top, self.traj, verbosity=0)
        while True:
            try:
                u.read_data()
                break
            except tables.HDF5ExtError:
                self.log.warning(
                    'HDF5ExtError raised. Likely because trying to read the '
                    'store at the same time as another process. Waiting 5 '
                    'seconds and trying again.')
                time.sleep(5)
                continue
        frames = u.select_frames(self.criteria, 'QM_frames')
        select = np.random.choice(frames)
        self.status['source_frame_num'] = int(select)
        system: mda.AtomGroup = u.select_atoms('all')
        xyz_name = self._base_name + '.xyz'
        with mda.Writer(xyz_name, system.n_atoms) as w:
            u.trajectory[select]
            for frag in u.atoms.fragments:
                mda.lib.mdamath.make_whole(frag)
                # This should at least make the molecules whole if not
                # necessarily in the correct unit cell together.
            w.write(system)
        self.log.info(f'Wrote xyz file from frame {select} to {xyz_name}')
        return pathlib.Path(xyz_name).resolve()

    def _move_reactant_atoms(self, xyz_path):
        self.log.debug('Moving reactant atoms (20 and 39) to '
                       f'{self.react_dist}')
        xyz_name = str(xyz_path)
        bak_name = xyz_name + '.bak'
        paratemp.copy_no_overwrite(xyz_name, bak_name)
        self.status['original_xyz'] = bak_name
        self.log.info(f'Copied original geometry to {bak_name}')
        xyz = XYZ(xyz_name)
        diff = xyz.coords[19] - xyz.coords[38]
        direction = diff / np.linalg.norm(diff)
        xyz.coords[19] = xyz.coords[38] + self.react_dist * direction
        xyz.write(xyz_name)
        self.log.info(f'Wrote updated xyz file to {xyz_name}')

    def new_calc(self):
        self.log.debug('Setting up a new calculation')
        self.current_lvl = 0
        if self.geometry is None:
            xyz_path = self._make_rand_xyz()
        else:
            self.log.debug(f'Using provided geometry from {self.geometry}')
            xyz_path = pathlib.Path(self.geometry).resolve()
        if self.react_dist:
            self._move_reactant_atoms(xyz_path)
        self.status['starting_xyz'] = str(xyz_path)
        if not xyz_path.exists():
            raise FileNotFoundError('Could not find start geometry that was '
                                    f'supposed to be at {xyz_path}')
        com_name = self._make_g_in(xyz_path)
        self._setup_and_run(com_name)

    def _setup_and_run(self, com_name):
        self.log.debug('Starting setup to run Gaussian')
        bn = self._base_name
        chk_ln_path = pathlib.Path(f'{bn}-running.chk').resolve()
        self.chk_ln_path = chk_ln_path
        self.status['chk_ln_path'] = str(chk_ln_path)
        chk_ln_path.symlink_to(self.scratch_path.joinpath(f'{bn}.chk'))
        self.log.info(f'Linked checkpoint file as {chk_ln_path}')
        if not self.resubmitted:
            self.resub_calc()
        self.status['manual_input'] = None
        self.status['manual_restart'] = None
        self.status['g_in_curr'] = com_name
        self.cleaned_up = False
        self.between_levels = False
        self.status['calc_cutoff'] = None
        self.status['gaussian_failed'] = None
        killed = self._run_gaussian(com_name)
        self.status['calc_cutoff'] = killed
        if killed:
            self.log.info('Exited from function running Gaussian because '
                          'SIGUSR2')
        else:
            try:
                self._check_normal_completion(self.output_scratch_path)
                self.log.info(f'Seemed to correctly finish level '
                              f'{self.current_lvl} calculation. Moving on to '
                              f'next level')
                self._advance_level()
            finally:
                self._copy_and_cleanup()
            self._next_calc()

    def _advance_level(self):
        self.log.debug(f'Advancing from {self.current_lvl}')
        self.between_levels = True
        self.current_lvl += 1

    def _make_g_in(self, xyz_path):
        self.log.debug(f'Making new Gaussian input from {xyz_path}')
        bn = self._base_name
        lvl = self.current_lvl
        com_name = f'{bn}-lvl{lvl}.com'
        try:
            ugt_dict = self.ugt_dicts[lvl]
        except IndexError:
            self.log.warning('Seems that there are no more calculation '
                             'levels to complete')
            raise self.NoMoreLevels
        tools.use_gen_template(
            out_file=com_name,
            xyz=str(xyz_path),
            job_name=bn,
            checkpoint=f'{bn}.chk',
            rwf=f'{bn}.rwf',
            nproc=self.n_slots, mem=self.mem,
            **ugt_dict
        )
        self.log.info('Wrote Gaussian input for '
                      f'level {lvl} job to {com_name}')
        self.status[f'g_in_{lvl}'] = com_name
        return com_name

    def _run_gaussian(self, com_name):
        self.log.debug('Doing final setup to run Gaussian')
        out_name = com_name.replace('com', 'out')
        com_path: pathlib.Path = self.cwd_path.joinpath(com_name)
        if not com_path.exists():
            raise FileNotFoundError('Gaussian input {} not found in '
                                    '{}'.format(com_name, self.cwd_path))
        out_path: pathlib.Path = self.scratch_path.joinpath(out_name)
        self.output_scratch_path = out_path
        self.status['output_scratch_path'] = str(out_path)
        old_sigusr1 = signal.signal(signal.SIGUSR1, self._signal_catch_done)
        old_sigusr2 = signal.signal(signal.SIGUSR2, self._signal_catch_time)
        cl = ['g16', ]
        killed = False
        with com_path.open('r') as f_in, out_path.open('w') as f_out:
            self.log.info('Starting Gaussian with input {} and writing '
                          'output to {}'.format(com_path, out_path))
            self.status['running'] = True
            proc = subprocess.Popen(cl, stdin=f_in, stdout=f_out,
                                    cwd=str(self.scratch_path))
            self.log.info('Started Gaussian; waiting for it to finish or '
                          'timeout')
            try:
                thread = threading.Thread(target=self._check_proc, args=(proc,))
                thread.start()
                signal.pause()
            except self.TimesUp:
                killed = True
                proc.terminate()  # Should be within `with` clause?
                self.log.info('Gaussian process terminated because of SIGUSR2')
            except self.GaussianDone:
                self.log.info('Gaussian process completed')
            finally:
                signal.signal(signal.SIGUSR1, old_sigusr1)
                signal.signal(signal.SIGUSR2, old_sigusr2)
        self.status['running'] = False
        return killed

    def _signal_catch_time(self, signum, frame):
        self.log.warning(f'Caught {signal.Signals(signum).name} signal! '
                         'Trying to quit Gaussian')
        raise self.TimesUp

    def _signal_catch_done(self, signum, frame):
        self.log.warning(f'Caught {signal.Signals(signum).name} signal! '
                         f'Likely, this was because Gaussian process exited')
        raise self.GaussianDone

    def _check_proc(self, proc):
        self.log.debug('Started process to check on Gaussian completion')
        while proc.poll() is None:
            time.sleep(15)
        self.log.warning('Gaussian process no longer running. Sending SIGUSR1')
        os.kill(os.getpid(), signal.SIGUSR1)

    def _copy_and_cleanup(self):
        self.log.debug('Attempting to copy back files and unlink chk file')
        com_name: str = self.status['g_in_curr']
        cc = self.status['calc_cutoff']
        killed: bool = True if cc is None else cc
        if killed:
            scratch_path = pathlib.Path(self.status['last_scratch_dir'])
        else:
            scratch_path = self.scratch_path
        if not killed:
            out_path = pathlib.Path(com_name.replace('com', 'out'))
            try:
                if self.status['gaussian_failed'] is True:
                    out_path = self._make_unique_output_path(
                        f'{out_path.stem}-failed')
            except KeyError:
                self.log.info('No key "gaussian_failed". Assuming it did '
                              'not...')
        else:
            out_path = self._make_unique_output_path(com_name[:-4])
        try:
            paratemp.copy_no_overwrite(str(self.output_scratch_path),
                                       str(out_path))
            self.log.debug(f'Copied back output file to {out_path}')
        except FileExistsError:
            if filecmp.cmp(str(self.output_scratch_path), str(out_path),
                           shallow=False):
                self.log.debug("Don't need to copy back output as it's already "
                               f"at {out_path}")
            else:
                self.log.error('Output files differ; unable to copy back.\n'
                               f'New: {self.output_scratch_path}\n'
                               f'Existing: {out_path}')
                raise
        if self.chk_ln_path.exists():
            self.chk_ln_path.unlink()
            self.log.debug(f'Unlinked checkpoint run file: {self.chk_ln_path}')
        chk_name = f'{self._base_name}.chk'
        scr_chk_path = scratch_path.joinpath(chk_name)
        if scr_chk_path.exists():
            shutil.copy(str(scr_chk_path), chk_name)
            self.log.debug(f'Copied back checkpoint file to {chk_name}')
        else:
            self.log.debug(f'chk file not found at {scr_chk_path} so not '
                           f'copied back')
        self.cleaned_up = True

    def _make_unique_output_path(self, com_base_name):
        outs = [str(p) for p in self.cwd_path.glob(com_base_name + '-*.out')]
        if not outs:
            new_out = f'{com_base_name}-1.out'
        else:
            outs.sort()
            outs.sort(key=len)
            new_out = re.sub(r'(\d+)\.out',
                             lambda m: '{}.out'.format(int(m.group(1)) + 1),
                             outs[-1])
        return pathlib.Path(new_out).resolve()

    def _check_normal_completion(self, filepath):
        self.log.debug('Attempting to check for completion status of Gaussian')
        output = subprocess.check_output(['tail', '-n', '1', str(filepath)],
                                         universal_newlines=True)
        if 'normal termination' not in output.lower():
            self.log.error(f'Abnormal termination of Gaussian job in output: '
                           f'{filepath}')
            self._qdel_next_job()
            self.status['gaussian_failed'] = True
            raise self.GaussianError('Gaussian did not finish normally. '
                                     f'See output: {filepath}')
        self.status['gaussian_failed'] = False
        self.log.info(f'Normal termination of Gaussian job! Output at '
                      f'{filepath}')

    def resub_calc(self):
        self.log.info(f'resubmitting job with the following commandline:\n'
                      f'{self.resub_cl}')
        proc = subprocess.run(self.resub_cl,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              universal_newlines=True)
        self.log.info(f'The following was returned from qsub:\n{proc.stdout}')
        if proc.returncode:
            self.log.exception('Resubmission of calculation failed with '
                               f'returncode {proc.returncode}')
        match = re.search(r'(\d+)\s\("(\w.*)"\)', proc.stdout)
        if match:
            self.next_job_id = match.group(1)
            self.resubmitted = True
        else:
            self.log.warning('Could not find submitted job id from qsub '
                             'command. Will not be able to cancel it if this '
                             'Calc is completed')

    def _make_resub_sh_and_cl(self):
        """
        Make command line for a calculation for resuming in another job

        Requires SGE_STDOUT_PATH and JOB_ID for running `qstat`
        :return: None
        """
        self.log.debug('Setting up for calculation resubmission')
        arg_d = dict(pe=f'omp {self.n_slots}', M='theavey@bu.edu', m='eas',
                     l=f'h_rt={self.h_rt}', N=self._base_name, j='y',
                     o=self.stdout_file, notify='', hold_jid=self.job_id)
        resub_dir_path = self.cwd_path.joinpath('resub_scripts')
        resub_dir_path.mkdir(exist_ok=True)
        sub_sh_path = resub_dir_path.joinpath(f'resub-{self._base_name}.sh')
        curr_file = pathlib.Path(__file__).resolve()
        with sub_sh_path.open('w') as sub_sh:
            sub_sh.write('#!/bin/bash -l\n\n')
            for key in arg_d:
                sub_sh.write(f'#$ -{key} {arg_d[key]}\n')
            sub_sh.write(obabel_module_lines)
            sub_sh.write(f'\n{curr_file} --restart {self._json_name} &\n\n')
            sub_sh.write(signal_catch_lines)
        self.log.info(f'Wrote resubmission script to {sub_sh_path}')
        self.resub_cl = ['qsub', str(sub_sh_path)]

    def _get_h_rt(self):
        """
        Find the amount of time requested for the currently running job

        Requires JOB_ID and uses `qstat`.
        :return:
        """
        self.log.debug('Attempting to find requested job run time')
        cl = ['qstat', '-j', self.job_id]
        output: str = subprocess.check_output(cl, universal_newlines=True)
        for line in output.splitlines():
            m = re.search(r'h_rt=(\d+)', line)
            if m:
                self.log.debug(f'Found required info: {m.group(0)}')
                return m.group(1)
        self.log.error('Could not find requested run time! Assuming 24 hours')
        return '24:00:00'

    def resume_calc(self):
        self.log.debug('Attempting to resume calculation')
        if not self.cleaned_up:
            self._copy_and_cleanup()
        try:
            manual_input = self.status['manual_input']
        except KeyError:
            manual_input = None
        try:
            manual_restart = self.status['manual_restart']
        except KeyError:
            manual_restart = None
        if manual_input is not None:
            self._setup_and_run(manual_input)
        elif self.between_levels:
            self._next_calc()
        else:
            if manual_restart is not None:
                com_name = manual_restart
            else:
                com_name = self._update_g_in_for_restart()
            self._copy_in_restart()
            self._setup_and_run(com_name)

    def _copy_in_restart(self):
        self.log.debug('Copying rwf and chk files to scratch for restart')
        bn = self._base_name
        old_rwf_path = self.last_scratch_path.joinpath(f'{bn}.rwf')
        if not old_rwf_path.exists():
            mes = f'Could not find old rwf file at {old_rwf_path}'
            self.log.error(mes)
            raise FileNotFoundError(mes)
        old_chk_path = self.last_scratch_path.joinpath(f'{bn}.chk')
        if not old_chk_path.exists():
            mes = f'Could not find old chk file at {old_chk_path}'
            self.log.error(mes)
            raise FileNotFoundError(mes)
        try:
            shutil.copy(str(old_rwf_path), str(self.scratch_path))
            shutil.copy(str(old_chk_path), str(self.scratch_path))
            self.log.info(f'Copied rwf and chk files from last scratch '
                          f'directory: {self.last_scratch_path}\nto node '
                          f'scratch dir: {self.scratch_path}')
        except shutil.SameFileError:
            self.log.info('Working on the same node; no need to copy files')

    def _update_g_in_for_restart(self):
        self.log.debug('Updating Gaussian input for restart')
        com_name = self.status['g_in_curr']
        lines = open(com_name, 'r').readlines()
        paratemp.copy_no_overwrite(com_name, com_name+'.bak')
        with open(com_name, 'w') as f_out:
            for line in lines:
                if '%mem=' in line:
                    line = f'%mem={self.mem}GB\n'
                elif line.startswith('#'):
                    line = '# Restart\n'
                f_out.write(line)
        os.remove(pathlib.Path(com_name+'.bak'))
        self.log.info(f'Updated Gaussian input to do a calculation restart '
                      f'and to use all the memory on this node')
        return com_name

    def _next_calc(self):
        self.log.debug('Moving on to next level calculation')
        out_path = pathlib.Path(self.status['g_in_curr']).with_suffix('.out')
        xyz_path_str = self._create_opt_xyz(out_path)
        try:
            com_name = self._make_g_in(xyz_path_str)
        except self.NoMoreLevels:
            self.log.info('No more calculation levels to complete! Completed '
                          f'all {self.current_lvl} levels')
            self._qdel_next_job()
            return None
        except FileExistsError:
            self.log.warning(f'Gaussian input file for level '
                             f'{self.current_lvl} already exists! This file '
                             f'will be used to start a new calculation')
            com_name = f'{self._base_name}-lvl{self.current_lvl}.com'
        self._setup_and_run(com_name)
        # This will get nested, but likely no more than twice (unless the
        # optimizations are very quick). This shouldn't be an issue,
        # and should never get near the recursion limit unless something goes
        # very wrong.

    def _qdel_next_job(self):
        self.log.debug('Deleting the re-submitted job from the queue')
        if self.next_job_id is not None:
            cl = ['qdel', self.next_job_id]
            output = subprocess.check_output(cl, stderr=subprocess.STDOUT)
            self.log.info('Cancelled job resubmission. qdel said: '
                          f'{output}')
        else:
            self.log.warning('Do not know job id of resubmission so '
                             'unable to delete it.')

    def _create_opt_xyz(self, out_path: pathlib.Path):
        self.log.debug('Converting output to xyz file for next level')
        xyz_path_str = str(out_path.with_suffix('.xyz'))
        cl = ['obabel', str(out_path), '-O',
              xyz_path_str]
        proc = subprocess.run(cl, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              universal_newlines=True)
        if proc.returncode or \
                '1 molecule converted' not in proc.stdout.lower():
            mes = (f'obabel failed to convert {out_path} to an xyz file. '
                   f'It said: {proc.stdout}')
            self.log.error(mes)
            raise subprocess.CalledProcessError(proc.returncode, cmd=cl,
                                                output=mes)
        self.log.info(f'Converted optimized structure to xyz file: '
                      f'{xyz_path_str}')
        return xyz_path_str

    class SignalMessage(Exception):
        pass

    class TimesUp(SignalMessage):
        pass

    class GaussianDone(SignalMessage):
        pass

    class NoMoreLevels(Exception):
        pass

    class GaussianError(Exception):
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
    * g_in_curr: str of the name of the currently running or most recent
        Gaussian input

    * current_lvl: int of current level of calculation running (max is len(
        ugt_dicts))

    * calc_cutoff: bool of whether the job finished or if it was cutoff
        because of running out of time.

    * cleaned_up: bool of whether linked files and old outputs have been
        cleaned up and copied back to the starting directory

    * chk_ln_path: str of path to where checkpoint file is linked in
        submission directory

    * output_scratch_path: str of path to where output is in the scratch
        directory

    * job_id: str of the job number from the sun grid system

    * between_levels: bool of if between levels (not in the middle of a
        calculation). Useful for figuring out where to restart.

    * gaussian_failed: bool of if Gaussian terminated abnormally. Will be None
        while Gaussian is running and before it gets checked.

    """
    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()
        if pathlib.Path(path).is_file():
            d = json.load(open(path, 'r'))
            super(StatusDict, self).__init__(d)
        else:
            super(StatusDict, self).__init__()
        self.log = logging.getLogger(self.__class__.__name__)
        self.temp_path = self.path.with_suffix('.json.new')

    def __setitem__(self, key, value):
        try:
            super(StatusDict, self).__setitem__(key, value)
            json.dump(self, open(self.temp_path, 'w'), indent=4)
            os.rename(str(self.temp_path), str(self.path))
        except Exception:
            self.log.exception('Exception raised when trying to write status '
                               'file!')
            raise


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--base_name', type=str,
                        help='base name for this calculation, likely not '
                             'including any index')
    parser.add_argument('-i', '--index', type=int,
                        help='index of this calculation')
    # TODO use mutually exclusive arguments here
    parser.add_argument('-x', '--xyz', type=str, default=None,
                        help='manual input geometry (will take precedence '
                             'over topology/trajectory frame selection)')
    parser.add_argument('-c', '--top', type=str, default=None,
                        help='topology/structure file (e.g., .gro, .xyz)')
    parser.add_argument('-f', '--trajectory', type=str, default=None,
                        help='trajectory file (e.g., .xtc, .trr, .dcd)')

    def parse_crit(kvv):
        k, vv = kvv.split('=')
        vs = tuple((float(v) for v in vv.split(',')))
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
    parser.add_argument('-g', '--ugt_dicts', type=str,
                        help='path to json file that parses to a list of '
                             'dicts of arguments for use_gen_template in '
                             'order to create inputs to Gaussian')
    parser.add_argument('--restart', default=None,
                        help='Path to status file for resuming an already '
                             'started calculation')
    p_args = parser.parse_args()
    if p_args.restart is not None:
        calc = Calc(status=p_args.restart)
    else:
        ugt_dicts = json.load(open(p_args.ugt_dicts, 'r'))
        criteria = None if p_args.criteria is None else dict(p_args.criteria)
        calc = Calc(base_name=p_args.base_name,
                    ind=p_args.index,
                    geometry=p_args.xyz,
                    top=p_args.top,
                    traj=p_args.trajectory,
                    criteria=criteria,
                    react_dist=p_args.react_dist,
                    ugt_dicts=ugt_dicts
                    )
    calc.run_calc()
