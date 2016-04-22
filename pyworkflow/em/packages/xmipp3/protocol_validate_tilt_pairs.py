# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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

from os.path import join
from pyworkflow.protocol.params import (PointerParam, StringParam)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import makePath
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
import pyworkflow.em.metadata as md


PROJECTION_MATCHING = 0
SIGNIFICANT = 1


class XmippProtValidateTiltPairs(ProtAnalysis3D):
    """    
    Validate if a set of Tilt Pairs have been well assigned.
    """
    _label = 'validate_Tilt_Pairs'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
       
       
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputUntilt', PointerParam, pointerClass='SetOfParticles', 
                      label="Input Untilt particles",  
                      help='Select the input untilt particles')
        form.addParam('inputTilt', PointerParam, pointerClass='SetOfParticles', 
                      label="Input Tilt particles",  
                      help='Select the input tilt particles')  
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='Symmetry group, only c1,..., cn, or d1, ... dn.'
                       'See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        deps = [] # store volumes steps id to use as dependencies for last step
        
        #convert step split and stores particletiltpair in two sets, untilted, and tilted
        convertId = self._insertFunctionStep('convertInputStep')
        
        fnuntilted = self._getExtraPath("untilted/untiltedpaticles.xmd")
        fntilted = self._getExtraPath("tilted/tiltedpaticles.xmd")
         
        volStepId = self._insertFunctionStep('validationStep', fnuntilted, fntilted, prerequisites=[convertId])
        deps.append(volStepId)
#         
#         self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self):
        """ Write the input images as a Xmipp metadata file.
        """
        makePath(self._getExtraPath("untilted"))
        makePath(self._getExtraPath("tilted"))
        
        untiltedSet = self.inputUntilt.get()
        tiltedSet = self.inputTilt.get()

        writeSetOfParticles(untiltedSet, self._getExtraPath("untilted/untiltedpaticles.xmd"))
        writeSetOfParticles(tiltedSet, self._getExtraPath("tilted/tiltedpaticles.xmd"))

    
    def anglesStep(self):
        print 'acabado'
    
    def validationStep(self, fnuntilt, fntilt):
        params =  ' --untilt %s' % fnuntilt
        params += ' --tilt %s' % fntilt
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --odir %s' % self._getExtraPath()

        self.runJob('xmipp_validation_tilt_pairs', params)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputTilt.get() and not self.inputTilt.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
    
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            size = 0
            for i, vol in enumerate(self._iterInputVols()):
                size +=1
            summary.append("Volumes to validate: *%d* " % size)
            summary.append("Angular sampling: %s" % self.angularSampling.get())
            summary.append("Significance value: %s" % self.significanceNoise.get())

        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and significant value of %f' % (self.angularSampling.get(), self.alpha.get()))
        return messages
    
    def _citations(self):
        return ['Not yet']
