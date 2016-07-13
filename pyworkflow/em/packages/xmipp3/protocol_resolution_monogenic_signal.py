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

from pyworkflow.protocol.params import (PointerParam, EnumParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from pyworkflow.em.data import Volume
from convert import readSetOfVolumes
from shutil import copyfile




class XmippProtMonoRes(ProtRefine3D):
    """    
    Given a map the protocol assings local resolutions to each pixel of the map.
    """
    _label = 'MonoRes - Monogenic Resolution'
    
    def __init__(self, *args, **kwargs):
        ProtRefine3D.__init__(self, *args, **kwargs)
        #self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
       
        form.addParam('inputVolume', PointerParam, pointerClass='Volume', 
                      label="Volume", 
                      help='Select a volume for determining its local resolucion.')
        
        form.addParam('Mask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask determines the those points of the sphere where the macromolecle is')

        line = form.addLine('Resolution Range (A)', 
                      help="If the user knows the range of resolutions or only a"
                      " range of frequency needs to be analysed", expertLevel=LEVEL_ADVANCED)
        line.addParam('minRes', FloatParam, default=1, label='Min')
        line.addParam('maxRes', FloatParam, default=100, label='Max')
        
        form.addParam('freqStep', FloatParam, label="Frequency Step",  default=0.01, expertLevel=LEVEL_ADVANCED,
                      help='The resolution is computed at frequencies between 0 and 0.5 px/A. this parameter determines'
                      'the sweep step')
        
        form.addParam('premask', BooleanParam, default=False,
                      label="Is the volume in a circular mask?",
                      help='Sometimes the volume is in an sphere, then this option ought to be selected')
        form.addParam('circularRadius', FloatParam, label="Radius Circular Mask",  condition = 'premask', 
                      default=50, 
                      help='This is the radius of the circular mask.')

#         form.addParam('gaussianfilter', FloatParam, label="Variance Gaussian filter",  default=0.5, expertLevel=LEVEL_ADVANCED,
#                       help='The map of resolution is filtered by means of a gaussian filter with' 
#                       'standard deviation provided by this parameter.')
        
        

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------


    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        deps = []
                 
        fnVol = self._getExtraPath('input_volume.vol')

        MS = self._insertFunctionStep('resolutionMonogenicSignalStep', fnVol, prerequisites=[convertId])
        
        #FilterSt = self._insertFunctionStep('medianFilterStep',  prerequisites=[MS])
        
        self._insertFunctionStep('createOutputStep', prerequisites=[MS])
       

    def convertInputStep(self):
        """ Read the input volume.
        """
        # Get the converted input micrographs in Xmipp format

        path_vol = self._getExtraPath() + '/input_volume.vol'
        vol_ = self.inputVolume.get().getFileName()
        copyfile(vol_,path_vol)
   
    def resolutionMonogenicSignalStep(self, fnVol):

        params =  ' --vol %s' % fnVol
        params +=  ' --mask %s' % self.Mask.get().getFileName()
        params +=  ' -o %s' % self._getExtraPath('MGresolution.vol')
        params +=  ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
        
        params +=  ' --stepW %f' % self.freqStep.get()
        params +=  ' --minRes %f' % self.minRes.get()
        params +=  ' --maxRes %f' % self.maxRes.get()
        if self.premask.get() is True:
            params +=  ' --circular_mask %f' % self.circularRadius.get()


        self.runJob('xmipp_resolution_monogenic_signal', params)
        
#     def medianFilterStep(self):
# 
#         params =  ' -i %s' % self._getExtraPath('MGresolution.vol')
#         params +=  ' --median'
#         params +=  ' -o %s' % self._getExtraPath('outputresolution.vol')
# 
#         self.runJob('xmipp_transform_filter', params)
    
    def createOutputStep(self):
        volume_path = self._getExtraPath('MGresolution.vol')
        
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate(self.inputVolume.get().getSamplingRate())
                
        readSetOfVolumes(volume_path, volumesSet)
        
        self._defineOutputs(outputVolume=volumesSet)
        self._defineSourceRelation(self.inputVolume, volumesSet)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        
        validateMsgs = []
        if not self.inputVolume.get().hasValue():
            validateMsgs.append('Please provide input volume.')  
        if not self.Mask.get().hasValue():
            validateMsgs.append('Please provide input mask.')         
        return validateMsgs

    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputParticles')):
            summary.append("Output tilpairs not ready yet.")
        else:
            summary.append("Three-uples of Tilt pairs angles assigned: %d" %self.outputParticles.__len__())
            #
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputParticles')):
            messages.append('An angular assignment of untilted and tilted particles is carried out [Publication: Not yet]')
        return messages
    
    def _citations(self):
        return ['Not yet']
    
    def getSummary(self):
        summary = []
        summary.append("Particles analyzed:")
        #summary.append("Particles picked: %d" %coordsSet.getSize())
        return "\n"#.join(summary)