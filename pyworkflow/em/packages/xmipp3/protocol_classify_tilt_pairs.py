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

from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.utils.path import makePath, removeBaseExt
from pyworkflow.em.data import SetOfParticles
from pyworkflow.em.data_tiltpairs import TiltPair, CoordinatesTiltPair,ParticlesTiltPair

from pyworkflow.em.metadata.constants import (MDL_ANGLE_Y, MDL_IMAGE, MDL_ANGLE_Y2, MDL_ANGLE_TILT, MDL_PARTICLE_ID, 
                                              MD_APPEND, MDL_MICROGRAPH_ID)
from pyworkflow.em.metadata import MetaData, getBlocksInMetaDataFile
from convert import readSetOfVolumes, getImageLocation

from convert import readAnglesFromMicrographs, readSetOfParticles
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs

from convert import writeSetOfParticles, writeSetOfCoordinates

TYPE_COORDINATES = 0
TYPE_PARTICLES = 1


class XmippProtClassifyTiltPairs(XmippProtParticlePickingPairs):
    """    
    From two sets of points (tilted and untilted) the protocol determines
    the affine transformation between these sets.
    """
    _label = 'Classify tiltpair'

    def __init__(self, *args, **kwargs):
        XmippProtParticlePickingPairs.__init__(self, *args, **kwargs)
        #self.stepsExecutionMode = params.STEPS_PARALLEL

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('Volume1', PointerParam, pointerClass='Volume',
                      label='Volume 1',
                      help='Select a volume.')
        
        form.addParam('Volume2', PointerParam, pointerClass='Volume',
                      label='Volume 2',
                      help='Select a volume.')

        form.addParam('tilpairparticles', PointerParam, pointerClass='ParticlesTiltPair', 
                      label="Set of particles tilt pairs", 
                      help='Select a set of particles tilt pairs.')
        
#         form.addParam('untiltedSet', PointerParam,
#                       pointerClass='SetOfCoordinates,SetOfParticles',
#                       label="Untilted input",
#                       help='Select the untilted input set, it can be either '
#                            'coordinates or particles (that contains coordinates.')
# 
#         form.addParam('tiltedSet', PointerParam,
#                       pointerClass='SetOfCoordinates,SetOfParticles',
#                       label="Tilted input",
#                       help='Select the tilted input set, it can be either '
#                            'coordinates or particles (that contains coordinates. '
#                            'It should be of the same type of the input untilted.')
#         
#         form.addParam('assign', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
#                       label="Perform an angular assignment?",
#                       help='The angular assignment can be performed by the protocol,'
#                       'Then, the classification is performed in a second step. Nevertheless,'
#                       'the angular assignment takes time')
        
        form.addParam('samplingAngle', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Sampling angle",
                      help='Sampling angle for angular assignment')

        form.addParam('maxShift', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Maximum Shift",
                      help='Maximum shift')
        
        form.addParam('padFactor', FloatParam, default=2, expertLevel=LEVEL_ADVANCED,
                      label="Padding Factor",
                      help='Padding Factor')
        
        

        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self.micsFn = self._getPath()
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        deps = []
        
#         if self.assign.get() is True:
        stepId  = self._insertFunctionStep('angularTiltAssignmentStep', 1,prerequisites=[convertId])
        stepId2 = self._insertFunctionStep('angularTiltAssignmentStep', 2,prerequisites=[stepId])
        stepId3 = self._insertFunctionStep('classifyStep', prerequisites=[stepId2])
        
        deps.append(stepId3)

        self._insertFunctionStep('createOutputStep', prerequisites=deps)


    def convertInputStep(self):
        """ Read the input metadatata. """
#         # Get the converted input micrographs in Xmipp format
#         makePath(self._getExtraPath("untilted"))
#         makePath(self._getExtraPath("tilted"))
# 
#         uSet = self.tilpairparticles.get().getUntilted()
#         tSet = self.tilpairparticles.get().getTilted()
#         self.imgsFnUntilted=self._getExtraPath('Untilted_particles.xmd')
#         self.imgsFnTilted=self._getExtraPath('Tilted_particles.xmd')
#         
#         writeSetOfParticles(uSet, self.imgsFnUntilted)
#         writeSetOfParticles(tSet, self.imgsFnTilted)
#         
        
        """ Read the input metadatata."""
        # Get the converted input micrographs in Xmipp format
        lastMicId = -1
        path_u = self._getExtraPath() + '/untilted_particles_in_blocks_input.xmd'
        path_t = self._getExtraPath() + '/tilted_particles_in_blocks_input.xmd'
        path_u_all = self._getExtraPath() + '/all_untilted_particles_input.xmd'
        path_t_all = self._getExtraPath() + '/all_tilted_particles_input.xmd'
        self.path_angles = self._getExtraPath() + '/Mic_angles.xmd'
        sangles = self.tilpairparticles.get().getCoordsPair().getAngles()
        md_angle = MetaData()
        for tpair in self.tilpairparticles.get().iterItems(orderBy='_untilted._micId'):
            micId = tpair.getUntilted().getMicId()
            if micId != lastMicId:
                if lastMicId != -1:
                    md_u.write("mic"+str(micId-1) +"@" + path_u, MD_APPEND)
                    md_t.write("mic"+str(micId-1) +"@" + path_t, MD_APPEND)
                
                #print micId
                angles = sangles[micId]
                (angleY, angleY2, angleTilt) = angles.getAngles()
                objId_angles = md_angle.addObject()
                md_angle.setValue(MDL_ANGLE_Y, float(angleY), objId_angles)
                md_angle.setValue(MDL_ANGLE_Y2, float(angleY2), objId_angles)
                md_angle.setValue(MDL_ANGLE_TILT, float(angleTilt), objId_angles)
                md_angle.setValue(MDL_MICROGRAPH_ID, long(micId), objId_angles)
                    
                md_u = MetaData()
                md_t = MetaData()
 
                lastMicId = micId
            objIdu = md_u.addObject()
            objIdt = md_t.addObject()
            #print getImageLocation(tpair.getUntilted())
            md_u.setValue(MDL_IMAGE, getImageLocation(tpair.getUntilted()), objIdu)
            md_t.setValue(MDL_IMAGE, getImageLocation(tpair.getTilted()), objIdt)
            md_u.setValue(MDL_PARTICLE_ID, long(tpair.getObjId()), objIdt)
            md_t.setValue(MDL_PARTICLE_ID, long(tpair.getObjId()), objIdt)
            md_u.setValue(MDL_MICROGRAPH_ID, long(micId), objIdt)
            md_t.setValue(MDL_MICROGRAPH_ID, long(micId), objIdt)
        

        md_angle.write(self.path_angles)
        md_u.write("mic"+str(micId) +"@" + path_u, MD_APPEND)
        md_t.write("mic"+str(micId) +"@" + path_t, MD_APPEND)
        
        mdblocks = getBlocksInMetaDataFile(path_u)

        md_u_aux = MetaData()
        md_t_aux = MetaData()
        md_u_all = MetaData()
        md_t_all = MetaData()

        #Creating two metadata(untilted and tilted) with all particles. They merges all blocks
        for blck in mdblocks:
            md_u_aux.read(blck+"@" + path_u)
            md_t_aux.read(blck+"@" + path_t)

            md_u_all.unionAll(md_u_aux)
            md_t_all.unionAll(md_t_aux)

        md_u_all.write(path_u_all)
        md_t_all.write(path_t_all)
        
        self.imgsFnUntilted=path_u_all
        self.imgsFnTilted=path_t_all
        
        
#         uSet = self.untiltedSet.get()
#         tSet = self.tiltedSet.get()
#         self.imgsFnUntilted=self._getExtraPath('Untilted_particles.xmd')
#         self.imgsFnTilted=self._getExtraPath('Tilted_particles.xmd')
#         
#         writeSetOfParticles(self.untiltedSet.get(), self.imgsFnUntilted)
#         writeSetOfParticles(self.tiltedSet.get(), self.imgsFnTilted)
        

    def classifyStep(self):

        params =  ' --md_Untilted1 %s' % self._getExtraPath('vol1_untilted.xmd')
        params += ' --md_Tilted1 %s' % self._getExtraPath('vol1_tilted.xmd')
        params += ' --md_Untilted2 %s' % self._getExtraPath('vol2_untilted.xmd')
        params += ' --md_Tilted2 %s' % self._getExtraPath('vol2_tilted.xmd')
        params += ' --odirMdUntilted_Vol1 %s' % self._getExtraPath('Untilted_assignment_vol1.xmd')
        params += ' --odirMdTilted_Vol1 %s' % self._getExtraPath('Tilted_assignment_vol1.xmd')
        params += ' --odirMdUntilted_Vol2 %s' % self._getExtraPath('Untilted_assignment_vol2.xmd')
        params += ' --odirMdTilted_Vol2 %s' % self._getExtraPath('Tilted_assignment_vol2.xmd')
        
        self.runJob('xmipp_classify_tilt_pairs', params, numberOfMpi=1)
        
    
    def angularTiltAssignmentStep(self, volume2assign):
        
        if (volume2assign == 1):
            RefVolume = self.Volume1.get().getFileName()
            fnbase = 'vol1' 
        else:
            RefVolume = self.Volume2.get().getFileName()
            fnbase = 'vol2'
            
        params =  ' --untiltparticles %s' % self.imgsFnUntilted
        params += ' --tiltparticles %s' % self.imgsFnTilted
        params += ' --odir %s' % self._getExtraPath()
        params += ' --untiltassignment %s' % self._getExtraPath(fnbase+'_untilted.xmd')
        params += ' --tiltassignment %s' % self._getExtraPath(fnbase+'_tilted.xmd')
        params += ' --angular_sampling %f' % self.samplingAngle.get()
        params += ' --maxshift %f' % self.maxShift.get()
        params += ' --pad %f' % self.padFactor.get()
        params += ' --initvol %s' % RefVolume
        params += ' --micangles %s' % self.path_angles

        self.runJob('xmipp_angular_tilt_pair_assignment', params)        
    


    def createOutputStep(self):
        USet = self._createSetOfParticles("UntiltedVol1")
        Upath = self._getExtraPath('Untilted_assignment_vol1.xmd')
        readSetOfParticles(Upath, USet)
        USet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        self._defineOutputs(outputUntiltedParticles1=USet)
        self._defineSourceRelation(self.tilpairparticles.get(), USet)
         
        TSet = self._createSetOfParticles("TiltedVol1")
        Tpath = self._getExtraPath('Tilted_assignment_vol1.xmd')
        readSetOfParticles(Tpath, TSet)
        TSet.setSamplingRate((self.tilpairparticles.get().getTilted().getSamplingRate()))
        self._defineOutputs(outputTiltedParticles2=TSet)
        self._defineSourceRelation(self.tilpairparticles.get(), TSet)

        #Define output ParticlesTiltPair 
        outputset_TiltPair_vol1 = ParticlesTiltPair(filename=self._getPath('particles_pairs_vol1.sqlite'))
        outputset_TiltPair_vol1.setTilted(TSet)
        outputset_TiltPair_vol1.setUntilted(USet)
        for imgU, imgT in izip(USet, TSet):
            outputset_TiltPair_vol1.append(TiltPair(imgU, imgT))

        #outputset_TiltPair_vol1.setCoordsPair(self.inputCoordinatesTiltedPairs.get())
        self._defineOutputs(outputParticlesTiltPair_vol1=outputset_TiltPair_vol1)
        self._defineSourceRelation(self.tilpairparticles, outputset_TiltPair_vol1)
        

        USet = self._createSetOfParticles("UntiltedVol2")
        Upath = self._getExtraPath('Untilted_assignment_vol2.xmd')
        readSetOfParticles(Upath, USet)
        USet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        self._defineOutputs(outputUntiltedParticles3=USet)
        self._defineSourceRelation(self.tilpairparticles.get(), USet)
        
        TSet = self._createSetOfParticles("TiltedVol2")
        Tpath = self._getExtraPath('Tilted_assignment_vol2.xmd')
        readSetOfParticles(Tpath, TSet)
        TSet.setSamplingRate((self.tilpairparticles.get().getTilted().getSamplingRate()))
        self._defineOutputs(outputTiltedParticles4=TSet)
        self._defineSourceRelation(self.tilpairparticles.get(), TSet)
        
        #Define output ParticlesTiltPair 
        outputset_TiltPair_vol2 = ParticlesTiltPair(filename=self._getPath('particles_pairs_vol2.sqlite'))
        outputset_TiltPair_vol2.setTilted(TSet)
        outputset_TiltPair_vol2.setUntilted(USet)
        for imgU, imgT in izip(USet, TSet):
            outputset_TiltPair_vol2.append(TiltPair(imgU, imgT))

        #outputset_TiltPair_vol2.setCoordsPair(self.inputCoordinatesTiltedPairs.get())
        self._defineOutputs(outputParticlesTiltPair_vol2=outputset_TiltPair_vol2)
        self._defineSourceRelation(self.tilpairparticles, outputset_TiltPair_vol2)
        
        

    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        errors = []
        uSet = self.tilpairparticles.get().getUntilted()
        tSet = self.tilpairparticles.get().getTilted()

        if (uSet is not None and tSet is not None and
            uSet.getClassName() != tSet.getClassName()):
            errors.append('Both untilted and tilted inputs should be of the '
                          'same type. ')

        return errors

    def _methods(self):
        messages = []
        if hasattr(self,'outputCoordinatesTiltPair'):
            messages.append('The assignment has been performed using and '
                            'affinity transformation [Publication: Not yet]')
        return messages

    def _citations(self):
        return ['Not yet']

