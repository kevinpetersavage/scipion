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

from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (IntParam, PointerParam, FloatParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.utils.path import makePath, replaceBaseExt
import xmipp
from convert import writeSetOfCoordinates, writeSetOfParticles, readSetOfParticles
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D


class XmippProtAngularTiltPairAssignment(ProtRefine3D):
    """    
    Given two sets of particles (untilted and tilted) the protocol assings angles to the 
    untilted particles and then determines by means of the tilt axis and tilt angle the
     orientation of the tilted particles.
    """
    _label = 'Angular tiltpair assignment '
    
    def __init__(self, *args, **kwargs):
        ProtRefine3D.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
#         form.addParam('tiltpair', PointerParam, pointerClass='MicrographsTiltPair',
#                       label="Micrograph tilt pair",  
#                       help='Select micrographs tilt pair.')
        
        form.addParam('tilpairparticles', PointerParam, pointerClass='ParticlesTiltPair', 
                      label="Set of particles tilt pairs", 
                      help='Select a set of particles tilt pairs.')
        
        form.addParam('refvol', PointerParam, pointerClass='Volume', 
                      label="Initial Volumen",   
                      help='Select an initial volumen for assigning angles.')
        
        form.addParam('samplingangle', FloatParam, default=5,
                      label="Angular sampling",  
                      help='The volume is projected, in many directions.'
                      'The angular sampling defines the distance between two '
                      'adjacents directions')   
        
        form.addParam('maxshift', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Maximum shift",  
                      help='Maximum shift for aligning images. If the obtained shift in the alignment'
                      'is higher than the maximum shift, the align between those images'
                      'will not be considered. If maxshift = -1, '
                      'there will not be considered a maximum shift')
        
#         form.addParam('numberIter', FloatParam, default=4, expertLevel=LEVEL_ADVANCED,
#                       label="Number of iterations",  
#                       help='Number of performed iterations in the reconstruction process')        
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        #self.micsFn = self._getPath('input_micrographs.xmd')
        # Convert input into xmipp Metadata format
        convertId=self._insertFunctionStep('convertInputStep')
        
        deps = []
            
        parUntilted = self.tilpairparticles.get().getUntilted()
        coordUntilted = parUntilted.getCoordinates()
        parTilted = self.tilpairparticles.get().getTilted()
        coordTilted = parTilted.getCoordinates()
        
        Unpath, Unname = split(parUntilted.getFileName())
        Unname, ext = splitext(Unname)
        Tpath, Tname = split(parTilted.getFileName())
        Tname, ext = splitext(Tname)
              
        
        fnUntilt = 'particles@'+self._getExtraPath("untilted/")+Unname+'.xmd'
        fnTilt = 'particles@'+self._getExtraPath("tilted/")+Tname+'.xmd'
        
        fnvol = self.refvol.get().getFileName()
        
        for mic_untilted in coordUntilted.iterMicrographs():
            for mic_tilted in coordTilted.iterMicrographs():
                if mic_untilted.getObjId() == mic_tilted.getObjId():
                    micName_untilted = mic_untilted.getFileName()
                    micName_tilted = mic_tilted.getFileName()
                    posFn_untilted = join(self._getExtraPath('untilted'), 
                                        replaceBaseExt(micName_untilted, "pos"))
                    posFn_tilted = join(self._getExtraPath('tilted'), 
                                        replaceBaseExt(micName_tilted, "pos"))
                    fnposUntilt = 'particles@'+posFn_untilted
                    fnposTilt = 'particles@'+posFn_tilted
                    stepId = self._insertFunctionStep('angularTiltAssignmentStep',
                                                      fnposUntilt, fnposTilt, 
                                                      fnvol, self._getExtraPath(), 
                                                      prerequisites=[convertId])
                 
                    deps.append(stepId)
            
        self._insertFunctionStep('createOutputStep', Unname, Tname, prerequisites=deps)
        

    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        makePath(self._getExtraPath("untilted"))
        makePath(self._getExtraPath("tilted"))
        
        U_set = self.tilpairparticles.get().getUntilted().getCoordinates()
        T_set = self.tilpairparticles.get().getTilted().getCoordinates()
        
        U_set_2 = self.tilpairparticles.get().getUntilted()
        T_set_2 = self.tilpairparticles.get().getTilted()
        
        writeSetOfCoordinates(self._getExtraPath("untilted"), U_set)
        writeSetOfCoordinates(self._getExtraPath("tilted"), T_set)
        
        writeSetOfParticles(U_set_2, self._getExtraPath("untilted/untilted_particles.xmd"))
        writeSetOfParticles(T_set_2, self._getExtraPath("tilted/tilted_particles.xmd"))
        
        

    
    def angularTiltAssignmentStep(self,fnposUntilt, fnposTilt, fnvol, Unpath):

        self.estimateTiltAxis(fnposUntilt, fnposTilt)
        
        fn = self._getPath()+'/input_micrographs.xmd'
        
        md = xmipp.MetaData(fn)

        alpha_t = md.getValue(xmipp.MDL_ANGLE_Y,1)
        alpha_u = md.getValue(xmipp.MDL_ANGLE_Y2,1)
        tiltmic = md.getValue(xmipp.MDL_ANGLE_TILT,1)
        
        unpath = self._getExtraPath("untilted/untilted_particles.xmd")
        tpath = self._getExtraPath("tilted/tilted_particles.xmd")
        params =  ' --untiltparticles %s' % unpath        
        params += ' --tiltparticles %s' % tpath
        params += ' --vol %s' % fnvol
        params += ' --alphaT %f' % alpha_t
        params += ' --tiltmic %f' % tiltmic
        params += ' --alphaU %f' % alpha_u
        params += ' --angular_sampling %f' % self.samplingangle.get()
        params += ' --maxshift %f' % self.maxshift.get()
        params += ' --odir %s' % Unpath


        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        self.runJob('xmipp_angular_tilt_pair_assignment', params, numberOfMpi=nproc,numberOfThreads=nT)
        

    def estimateTiltAxis(self, fnposUntilt, fnposTilt):

        params =  ' --untilted %s' % fnposUntilt        
        params += ' --tilted %s' % fnposTilt
        params += ' -o %s' % self._getPath()+'/input_micrographs.xmd'

        self.runJob('xmipp_angular_estimate_tilt_axis', params)
        
          
    def createOutputStep(self, Unname, Tname):
       
#         #Defining paths with metadatafiles       
        unpath = self._getExtraPath('untilted_particles_angular_assignment.xmd')
        tpath = self._getExtraPath('tilted_particles_angular_assignment.xmd')
        AllTpath = self._getExtraPath("All_tilted_assignment.xmd")

        md_U = xmipp.MetaData(unpath)
        md_T = xmipp.MetaData(tpath)

        md_U.unionAll(md_T)
        md_U.write(AllTpath)

        #Output of untilt particles - angular assignment
        output_untilt = self._createSetOfParticles(suffix='untilted')
        
        inputUntilted = self.tilpairparticles.get().getUntilted()
        output_untilt.copyInfo(inputUntilted)
        readSetOfParticles(unpath, output_untilt)

        self._defineOutputs(outputParticlesUntilted=output_untilt)
        self._defineSourceRelation(inputUntilted, output_untilt)
  
        #Output of All particles(untilted and tilted) - angular assignment
        output_all = self._createSetOfParticles(suffix='all')
        output_all.copyInfo(inputUntilted)
        readSetOfParticles(AllTpath, output_all)

        self._defineOutputs(outputParticlesAll=output_all)
        self._defineSourceRelation(inputUntilted, output_all)
        

        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        
        validateMsgs = []
        if self.tilpairparticles.get().getUntilted() and not self.tilpairparticles.hasValue():
            validateMsgs.append('Please provide input particles.')  
        if self.tilpairparticles.get().getTilted() and not self.tilpairparticles.hasValue():
            validateMsgs.append('Please provide input particles.')         
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