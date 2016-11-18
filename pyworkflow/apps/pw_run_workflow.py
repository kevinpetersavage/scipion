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
import time


def run(name, workflow, launch_timeout):
    manager = Manager()
    project = manager.createProject(name)
    protocols = project.loadProtocols(workflow)

    graph = project.getGraphFromRuns(protocols.values())
    nodes = graph.getRoot().iterChildsBreadth()

    protocols_by_id = {protocol.strId():protocol for protocol in protocols.values()}

    for node in nodes:
        protocol = protocols_by_id[node.getName()]
        parents = collect_parents(node, protocols_by_id)

        launch_when_ready(parents, project, protocol, launch_timeout)


def collect_parents(node, protocols_by_id):
    return {protocols_by_id[parent_node.getName()]
            for parent_node in node.getParents() if not parent_node.isRoot()}


def launch_when_ready(parents, project, protocol, launch_timeout):
    start = time.time()
    while time.time() < start + launch_timeout:
        update_protocols(parents, project)

        errors = protocol.validate()
        if errors:
            print errors
            time.sleep(1)
        else:
            project.launchProtocol(protocol)
            break

        check_protocol(protocol)


def update_protocols(parents, project):
    for parent in parents:
        project._updateProtocol(parent)


def check_protocol(protocol):
    if protocol.isFailed():
        print "\n>>> ERROR running protocol %s" % protocol.getRunName()
        print "    FAILED with error: %s\n" % protocol.getErrorMessage()
        raise Exception("ERROR launching protocol.")


if __name__ == '__main__':
    project_name = sys.argv[1]
    workflow = sys.argv[2]
    if len(sys.argv) > 3:
        timeout = int(sys.argv[3])
    else:
        timeout = 60*60
    run(project_name, workflow, timeout)