#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module is responsible for running workflows from a json description.
"""

from pyworkflow.manager import Manager
import sys

def run(project_name, workflow):
    manager = Manager()
    project = manager.createProject(project_name)
    protocols = project.loadProtocols(workflow)
    for id, protocol in protocols.iteritems():
        project.launchProtocol(protocol, wait=True)

        if protocol.isFailed():
            print "\n>>> ERROR running protocol %s" % protocol.getRunName()
            print "    FAILED with error: %s\n" % protocol.getErrorMessage()
            raise Exception("ERROR launching protocol.")

        if not protocol.isFinished():
            print "\n>>> ERROR running protocol %s" % protocol.getRunName()
            raise Exception("ERROR: Protocol not finished")


if __name__ == '__main__':
    project_name = sys.argv[1]
    workflow = sys.argv[2]
    run(project_name, workflow)