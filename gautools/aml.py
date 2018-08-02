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
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
import signal


class Calc(object):

    def __init__(self, base_name, ind, top, traj, criteria, ugt_dicts):
        """

        :param str base_name:
        :param int ind:
        :type top: pathlib.Path or str
        :param top:
        :type traj: pathlib.Path or str
        :param traj:
        :param dict criteria:
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
        self.status = json.load(open(self._json_name, 'r')) if pathlib.Path(
            self._json_name).is_file() else dict()

    def run_calc(self, base_name, ind, top, traj, criteria, ugt_dicts):
        """

        :param str base_name:
        :param int ind:
        :type top: pathlib.Path or str
        :param top:
        :type traj: pathlib.Path or str
        :param traj:
        :param dict criteria:
        :param List[dict] ugt_dicts:
        :return:
        """
        # TODO: get ind (index) from SGE_TASK_ID?
        # do this if not given, and save it to args for passing on?
        args = base_name, ind, top, traj, criteria, ugt_dicts
        # TODO: save args to disk?
        # Could then give these defaults, and if it's a continuation, just load
        # all the arguments from disk. Need to think about how this will
        # be used, exactly
        _base_name = '{}-{}'.format(base_name, ind)
        log = logging.getLogger('{}.log'.format(_base_name))
        _json_name = '{}.json'.format(_base_name)
        status = json.load(open(_json_name, 'r')) if pathlib.Path(
            _json_name).is_file() else dict()
        # TODO: catch job-ending warning here?
        signal.signal(signal.SIGUSR2, self.resub_calc)
        if status:
            log.info('loaded previous status file: {}'.format(_json_name))
            self.resume_calc(status, *args)
        else:
            log.warning('No previous status file found. '
                        'Starting new calculation?')
            self.new_calc(status, *args)

    def resume_calc(self, status, base_name, ind, top, traj, criteria,
                    ugt_dicts):
        pass

    def new_calc(self, status, base_name, ind, top, traj, criteria, ugt_dicts):
        pass

    def resub_calc(self, signum, frame):
        pass
