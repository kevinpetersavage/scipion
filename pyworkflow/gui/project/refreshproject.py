#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import threading
import time


class RefreshProject(threading.Thread):
    INITIAL_AUTOREFRESH = 3

    def __init__(self, threadName, refreshCallback, project):
        super(RefreshProject, self).__init__(name=threadName)
        self.threadID = threadName
        self._autoRefreshCounter = RefreshProject.INITIAL_AUTOREFRESH
        self.cancelRefresh = False
        self.refreshCallback = refreshCallback
        self.refreshing = False
        self.project = project
        self._lock = threading.Lock()

    def run(self):

        from pyworkflow.utils import envVarOn

        # Cancel automatic refresh if evn variable is on.
        if envVarOn('DO_NOT_AUTO_REFRESH'):
            return

        while True:

            # If paused
            # self._lock.acquire()

            # Wait
            self.waitUntilNextRefresh()

            # If need to cancel refresh (due to external refresh during sleeptime)
            if self.cancelRefresh:
                print "Manual refresh detected"
                self.cancelRefresh = False

                # self._lock.release
                continue
            else:

                # Double the refresh time up to 30 minutes
                self._autoRefreshCounter = min(self._autoRefreshCounter * 2, 1800)

                # Refresh project
                print 'Automatic refresh'
                self._refresh(True, False)

            # self._lock.release()

    def waitUntilNextRefresh(self):

        waitingTime = 0
        keepWaiting = True

        print "Waiting %s seconds to star next refresh." % self._autoRefreshCounter

        while not self.cancelRefresh and keepWaiting:

            # Sleep in second intervals
            time.sleep(self.INITIAL_AUTOREFRESH)
            waitingTime += self.INITIAL_AUTOREFRESH
            keepWaiting = waitingTime < self._autoRefreshCounter

    def _refresh(self, refresh, checkPids):

        self.refreshing = True

        # Refresh project
        runs = self.project.getRunsGraph(refresh=refresh, checkPids=checkPids)

        self.refreshing = False

        # Call the callback
        self.refreshCallback(runs)

    """Public method to be called at any refresh need, it will reset the autorefresh"""
    def refresh(self, refresh=True, checkPids=False):

        # reset autorefresh counter
        self._autoRefreshCounter = self.INITIAL_AUTOREFRESH
        print "Autorefresh counter reset by external refresh request."

        # if refresh is happening ...
        if self.refreshing:
            # ... do nothing
            # this may not be running with same parameters (checkPids or refresh).
            print "refreshing ongoing"
            return

        else:
            # Cancel refresh
            self.cancelRefresh = True
            self._refresh(refresh, checkPids)

    def pause(self):
        self._autoRefreshCounter = self.INITIAL_AUTOREFRESH
        # self._lock.acquire()


    def resume(self):
        # self._lock.release()
        pass