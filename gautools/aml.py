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
import thtools
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
        self._base_name = '{}-{}'.format(base_name, ind)
        self.log = logging.getLogger('{}.log'.format(self._base_name))
        self._json_name = '{}.json'.format(self._base_name)
        self._status = StatusDict(self._json_name)
        self.mem, self.node, self.last_node = None, None, None

    @property
    def status(self):
        return self._status

    def _startup_tasks(self):
        self.mem = thtools.job_tools.get_node_mem()
        node = os.environ['HOSTNAME'].split('.')[0]
        self.node = node
        n_slots = int(os.environ['NSLOTS'])
        self.log.info('Running on {} using {} cores and up to {} GB '
                      'mem'.format(node, n_slots, self.mem))
        if self.status:
            self.last_node = self.status['current_node']
            self.status['last_node'] = self.last_node
            node_list = self.status['node_list']
        else:
            node_list = []
        self.status['node_list'] = node_list + [node]
        self.status['current_node'] = node

    def run_calc(self):
        """


        :return:
        """
        self._startup_tasks()
        signal.signal(signal.SIGUSR2, self.resub_calc)
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
        com_name = bn + '-0.com'
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
        pass

    def resub_calc(self, signum, frame):
        pass

    def resume_calc(self):
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
