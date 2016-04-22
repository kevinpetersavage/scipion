# **************************************************************************
# *
# * Authors:  Jose Luis Vilas Prieto (jlvilas@cnb.csic.es), april 2016
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
# from pyworkflow.em.viewer import DataView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

from pyworkflow.em.metadata import MetaData
from pyworkflow.em.metadata.constants import MDL_TILT_AXIS_X, MDL_TILT_AXIS_Y, MDL_TILT_AXIS_Z, MDL_ROTATION_ANGLE
import matplotlib.pyplot as plt
from protocol_validate_tilt_pairs import XmippProtValidateTiltPairs
from matplotlib.figure import Figure

class XmippValidationTiltPairsViewer(XmippViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer volume strain'
    _targets = [XmippProtValidateTiltPairs]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        import os
        fnCmd = self.protocol._getExtraPath('scatter_md.xmd')
        print fnCmd
        if os.path.exists(fnCmd):
            self.plotScatterTiltAxis(fnCmd)
            plt.figure()
            self.plothistogram(fnCmd)
#         fnCmd = self.protocol._getPath('scatter_md.xmd')
#         if os.path.exists(fnCmd):
#             self.plotscatter(fnCmd)
            
    def plotScatterTiltAxis(self, fnCmd):
        import numpy as np
        md = MetaData()
        md.read(fnCmd)
        k = 0
        lenmd = md.size()
        tiltaxis = np.zeros(shape=(lenmd,3))
        print tiltaxis
        for idx in md:
            tiltaxis0_ = md.getValue(MDL_TILT_AXIS_X, idx)
            tiltaxis1_ = md.getValue(MDL_TILT_AXIS_Y, idx)
            tiltaxis2_ = md.getValue(MDL_TILT_AXIS_Z, idx)

            tiltaxis[k][0] = tiltaxis0_
            tiltaxis[k][1] = tiltaxis1_
            tiltaxis[k][2] = tiltaxis2_
            k = k + 1
        self.sphere(tiltaxis[:][0], tiltaxis[:][1], tiltaxis[:][2])


    def plothistogram(self, fnCmd):
        md = MetaData()
        md.read(fnCmd)
        k = 0
        theta = []
        for idx in md:
            theta_ = md.getValue(MDL_ROTATION_ANGLE, idx)
            theta.append(theta_)
            k = k + 1;
        print 'plotting...'
        plt.hist(theta)
        plt.title("Tilt Angle Histogram")
        plt.xlabel("Tilt Angle")
        plt.ylabel("Counts")
        plt.show()
        
        print 'plotted!'
        
        
    def sphere(self, data_xx, data_yy, data_zz):
        import numpy as np
        # Create a sphere
        r = 1
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
        x = r*sin(phi)*cos(theta)
        y = r*sin(phi)*sin(theta)
        z = r*cos(phi)
        
        #Import data
        xx = data_xx
        yy = data_yy
        zz = data_zz
        
        #Set colours and render
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.plot_surface(
            x, y, z,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)
        
        ax.scatter(xx,yy,zz,color="k",s=20)
        
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])
        ax.set_zlim([-1,1])
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.show()