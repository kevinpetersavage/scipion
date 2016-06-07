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

from itertools import izip
from os.path import split, splitext, join

from pyworkflow.protocol.params import (PointerParam, EnumParam, IntParam, FloatParam, BooleanParam, StringParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.em.metadata import MetaData, getBlocksInMetaDataFile
from convert import readSetOfVolumes, getImageLocation
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from pyworkflow.em.data import Volume
from shutil import copyfile
import pyworkflow.em as em
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles



class XmippProtAddNoise(ProtRefine3D):
    """    
    Given two sets of particles (untilted and tilted) the protocol assings angles to the 
    untilted particles and then determines by means of the tilt axis and tilt angle the
     orientation of the tilted particles.
    """
    _label = 'Add noise'
    
    def __init__(self, *args, **kwargs):
        ProtRefine3D.__init__(self, *args, **kwargs)
        #self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
       
        form.addParam('input', PointerParam, pointerClass='Volume', 
                      label="Input", 
                      help='Select a volume or image.')
        
        form.addParam('noiseType', EnumParam, choices=['Gaussian', 'Student', 'Uniform'],
                      default = 0,
                      label="Noise Type")
        
        form.addParam('gaussianStd', FloatParam, default=0.08, 
                      condition='noiseType == %d' % 0,
                      label="Standard Deviation", 
                      help='Please, introduce the standard deviation value.'
                      'Mean value can be change in advanced mode.')
        
        form.addParam('gaussianMean', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      condition='noiseType == %d' % 0,
                      label="Mean", 
                      help='Please, introduce the mean value (default = 0).')
        
        form.addParam('studentDf', FloatParam, default=1, 
                      condition='noiseType == %d' % 1,
                      label="Degree of Freedom", 
                      help='Please, introduce the Degree of Freedom.'
                      'Mean value can be change in advanced mode.')
        
        form.addParam('studentStd', FloatParam, default=0.08, 
                      condition='noiseType == %d' % 1,
                      label="Standard Deviation", 
                      help='Please, introduce the standard deviation value.'
                      'Mean value can be change in advanced mode.')
        
        form.addParam('studentMean', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      condition='noiseType == %d' % 1,
                      label="Mean", 
                      help='Please, introduce the mean value (default = 0).')
        
        form.addParam('uniformMin', FloatParam, default=0, 
                      condition='noiseType == %d' % 2,
                      label="Minimum Value", 
                      help='Please, introduce the minimum value. (default = 0)')
        
        form.addParam('uniformMax', FloatParam, default=1, 
                      condition='noiseType == %d' % 2,
                      label="Maximum Value", 
                      help='Please, introduce the maximum value (default = 1).')
        
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------


    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        deps = []
         
        self._insertFunctionStep('addNoiseStep')
# 
# 
        self._insertFunctionStep('createOutputStep')


    def addNoiseStep(self):
        if self.noiseType.get() == 0:
            Noise_Type = 'gaussian'
            Noise_Params = '%f %f' % (self.gaussianStd.get(), self.gaussianMean.get())
        if self.noiseType.get() == 1:
            Noise_Type = 'student'
            Noise_Params = '%f %f %f' % (self.studentDf.get(), self.studentStd.get(), self.studentMean.get())
        if self.noiseType.get() == 2:
            Noise_Type = 'uniform'
            Noise_Params = '%f %f' % (self.uniformMin.get(), self.uniformMax.get())
        
        args = "-i %s --type %s %s -o %s" % (self._getExtraPath('init_vol.vol'), Noise_Type, Noise_Params,
                                              self._getExtraPath('Noisy_volume.vol'))
        self.runJob('xmipp_transform_add_noise', args, numberOfMpi=1)


    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        volume_init = self.input.get().getFileName()
        copyfile(volume_init,self._getExtraPath('init_vol.vol'))


    def createOutputStep(self):
        #Output Volume
        
        fnOutVol = self._getExtraPath('init_vol.vol')
        
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate( self.input.get().getSamplingRate() )
        
        readSetOfVolumes(fnOutVol, volumesSet)
        
        self._defineOutputs(outputVolume=volumesSet)
        self._defineSourceRelation(self.input.get(), volumesSet)

    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        
        validateMsgs = []
        if self.input.get() and not self.input.hasValue():
            validateMsgs.append('Please provide input volume.')  
        return validateMsgs

    def _summary(self):
        summary = []

        if self.noiseType.get() == 0:
            Noise_Type = 'gaussian'
        if self.noiseType.get() == 1:
            Noise_Type = 'student'
        if self.noiseType.get() == 2:
            Noise_Type = 'uniform'

        if  (not hasattr(self,'outputVolume')):
            summary.append("Output volume not ready yet.")
        else:
            summary.append("Volume with %s noise has been obtained" % (Noise_Type) )
            #
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolume')):
            messages.append('Noisy volume has been obtained')
        return messages
    
    def _citations(self):
        return ['Do not apply']
    
    def getSummary(self):
        summary = []
        summary.append("Particles analyzed:")
        #summary.append("Particles picked: %d" %coordsSet.getSize())
        return "\n"#.join(summary)