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
                                              MD_APPEND, MDL_MICROGRAPH_ID, MDL_ANGLE_ROT, MDL_ANGLE_PSI, MDL_MAXCC, MDL_SHIFT_X, MDL_SHIFT_Y)
from pyworkflow.em.metadata import MetaData, getBlocksInMetaDataFile
from convert import readSetOfVolumes, getImageLocation

from convert import readAnglesFromMicrographs, readSetOfParticles
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs

from convert import writeSetOfParticles, writeSetOfCoordinates



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

        form.addParam('tilpairparticles', PointerParam, pointerClass='ParticlesTiltPair', 
                      label="Set of particles tilt pairs", 
                      help='Select a set of particles tilt pairs.')

        form.addParam('inputVolume', BooleanParam, default=True, 
                      label="Is there a reference volume?",
                      help='The particles will be correlated with this volume, then, two references volumes will be generated')
                
        form.addParam('TwoVolumes', BooleanParam, default=False, condition = 'inputVolume', 
                      label="Use previous volumes?",
                      help='The angular assignment can be performed by the protocol,'
                      'Then, the classification is performed in a second step. Nevertheless,'
                      'the angular assignment takes time')

        form.addParam('Volume1', PointerParam, pointerClass='Volume', condition = 'inputVolume',
                      label='Volume 1',
                      help='Select a volume.')

        form.addParam('Volume2', PointerParam, pointerClass='Volume', condition = 'TwoVolumes', 
                      label='Volume 2',
                      help='Select a volume.')

        form.addParam('iters', FloatParam, pointerClass='ParticlesTiltPair', 
                      label="Iterations", 
                      help='Number of iterations for clasifying.')
        
        form.addParam('randomVolumes', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Perform an angular assignment?",
                      help='The angular assignment can be performed by the protocol,'
                      'Then, the classification is performed in a second step. Nevertheless,'
                      'the angular assignment takes time')
        
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
        
        if self.inputVolume.get() is True:
            if self.TwoVolumes.get() is False:
                stepIdTV  = self._insertFunctionStep('angularTiltAssignmentStep', 1, 0, prerequisites=[convertId])
                stepIdTV2 = self._insertFunctionStep('InitialVolumesStep', prerequisites=[stepIdTV])
                stepIdTV3 = self._insertFunctionStep('medianFilterStep', 1, 0, prerequisites=[stepIdTV2])
                stepIdTV4 = self._insertFunctionStep('medianFilterStep', 2, 0, prerequisites=[stepIdTV3])
#             else:
#                 stepId  = self._insertFunctionStep('angularTiltAssignmentStep', 1,prerequisites=[convertId])
#                 stepId2 = self._insertFunctionStep('angularTiltAssignmentStep', 2,prerequisites=[stepId])
#                 stepId3 = self._insertFunctionStep('classifyStep', prerequisites=[stepId2])

        for iterNum in range(int(self.iters.get())):
            stepId  = self._insertFunctionStep('angularTiltAssignmentStep', 1, iterNum+1, prerequisites=[convertId])
            stepId2 = self._insertFunctionStep('angularTiltAssignmentStep', 2, iterNum+1, prerequisites=[stepId])
            stepId3 = self._insertFunctionStep('classifyStep', iterNum+1, prerequisites=[stepId2])
            stepId4 = self._insertFunctionStep('reconstructFourierStep', 1, iterNum+1, prerequisites=[stepId3])
            stepId5 = self._insertFunctionStep('reconstructFourierStep', 2, iterNum+1, prerequisites=[stepId4])
            stepId6 = self._insertFunctionStep('medianFilterStep', 1, iterNum+1, prerequisites=[stepId5])
            stepId7 = self._insertFunctionStep('medianFilterStep', 2, iterNum+1, prerequisites=[stepId6])
            

        self._insertFunctionStep('createOutputStep', prerequisites=[stepId7])


    def convertInputStep(self):
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
                print angleY
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

     
    def InitialVolumesStep(self):
        
        """ Read the input metadatata."""
        # Get the converted input micrographs in Xmipp format
        TSet = self.tilpairparticles.get().getTilted()
        NumElems =  TSet.getSize()
        
        md_TSet = MetaData()
        md_TSet.read(self._getExtraPath('tilted_assignment_vol1_iter_0.xmd'))

        md_TSet1 = MetaData()
        md_TSet2 = MetaData()
        objId = 0
        #TSetPrime.iterItems(orderBy='_maxCC').printAll()
        md_TSet.sort(MDL_MAXCC)
        
        
        for tpair in md_TSet:
            objId = objId + 1
            #objId = md_TSet.addObject()

            fnImg = md_TSet.getValue(MDL_IMAGE, objId)
            par_id = md_TSet.getValue(MDL_PARTICLE_ID, objId)
            rot = md_TSet.getValue(MDL_ANGLE_ROT, objId)
            tilt = md_TSet.getValue(MDL_ANGLE_TILT, objId)
            psi = md_TSet.getValue(MDL_ANGLE_PSI, objId)
            corr = md_TSet.getValue(MDL_MAXCC, objId)
            Sx = md_TSet.getValue(MDL_SHIFT_X, objId)
            Sy = md_TSet.getValue(MDL_SHIFT_Y, objId)
            if (objId > 0.5*NumElems):
                idx1 = md_TSet1.addObject()
                md_TSet1.setValue(MDL_IMAGE, fnImg, idx1)
                md_TSet1.setValue(MDL_PARTICLE_ID, par_id, idx1)
                md_TSet1.setValue(MDL_ANGLE_ROT, rot, idx1)
                md_TSet1.setValue(MDL_ANGLE_TILT, tilt, idx1)
                md_TSet1.setValue(MDL_ANGLE_PSI, psi, idx1)
                md_TSet1.setValue(MDL_MAXCC, corr, idx1)
                md_TSet1.setValue(MDL_SHIFT_X, Sx, idx1)
                md_TSet1.setValue(MDL_SHIFT_Y, Sy, idx1)
            else:
                idx2 = md_TSet2.addObject()
                md_TSet2.setValue(MDL_IMAGE, fnImg, idx2)
                md_TSet2.setValue(MDL_PARTICLE_ID, par_id, idx2)
                md_TSet2.setValue(MDL_ANGLE_ROT, rot, idx2)
                md_TSet2.setValue(MDL_ANGLE_TILT, tilt, idx2)
                md_TSet2.setValue(MDL_ANGLE_PSI, psi, idx2)
                md_TSet2.setValue(MDL_MAXCC, corr, idx2)
                md_TSet2.setValue(MDL_SHIFT_X, Sx, idx2)
                md_TSet2.setValue(MDL_SHIFT_Y, Sy, idx2)
        md_TSet1.write(self._getExtraPath('tilted_assignment_vol1.xmd'))
        md_TSet2.write(self._getExtraPath('tilted_assignment_vol2.xmd'))
        
        params =  '  -i %s' % self._getExtraPath('tilted_assignment_vol1.xmd')
        params += '  -o %s' % self._getExtraPath('Reference_volume1_iter0.vol')
        params += '  --sym %s' % 'c1'
        params += '  --max_resolution %f' % 0.5
        params += '  --padding %f' %self.padFactor.get()
        self.runJob('xmipp_reconstruct_fourier', params);
        

        params =  '  -i %s' % self._getExtraPath('tilted_assignment_vol2.xmd')
        params += '  -o %s' % self._getExtraPath('Reference_volume2_iter0.vol')
        params += '  --sym %s' % 'c1'
        params += '  --max_resolution %f' % 0.5
        params += '  --padding %f' %self.padFactor.get()
        self.runJob('xmipp_reconstruct_fourier', params);


    def classifyStep(self, iterNum):

        Untilted1Path = 'untilted_assignment_vol1_iter_%d.xmd' % iterNum
        Tilted1Path = 'tilted_assignment_vol1_iter_%d.xmd' % iterNum
        Untilted2Path = 'untilted_assignment_vol2_iter_%d.xmd' % iterNum
        Tilted2Path = 'tilted_assignment_vol2_iter_%d.xmd' % iterNum
        untilted1_classPath = 'Untilted_classification_vol1_iter_%d.xmd' % iterNum
        tilted1_classPath = 'Tilted_classification_vol1_iter_%d.xmd' % iterNum
        untilted2_classPath = 'Untilted_classification_vol2_iter_%d.xmd' % iterNum
        tilted2_classPath = 'Tilted_classification_vol2_iter_%d.xmd' % iterNum
        
        params =  ' --md_Untilted1 %s' % self._getExtraPath(Untilted1Path)
        params += ' --md_Tilted1 %s' % self._getExtraPath(Tilted1Path)
        params += ' --md_Untilted2 %s' % self._getExtraPath(Untilted2Path)
        params += ' --md_Tilted2 %s' % self._getExtraPath(Tilted2Path)
        params += ' --odirMdUntilted_Vol1 %s' % self._getExtraPath(untilted1_classPath)
        params += ' --odirMdTilted_Vol1 %s' % self._getExtraPath(tilted1_classPath)
        params += ' --odirMdUntilted_Vol2 %s' % self._getExtraPath(untilted2_classPath)
        params += ' --odirMdTilted_Vol2 %s' % self._getExtraPath(tilted2_classPath)
        
        self.runJob('xmipp_classify_tilt_pairs', params, numberOfMpi=1)
        
    def reconstructFourierStep(self, volume2assign, iterNum):
        
        iPath = 'Tilted_classification_vol%d_iter_%d.xmd' % (volume2assign, iterNum)
        oPath = 'Reference_volume%d_iter%d.vol' % (volume2assign, iterNum)

        params =  '  -i %s' % self._getExtraPath(iPath)
        params += '  -o %s' % self._getExtraPath(oPath)
        params += '  --sym %s' % 'c1'
        params += '  --max_resolution %f' % 0.5
        params += '  --padding %f' %self.padFactor.get()
        self.runJob('xmipp_reconstruct_fourier', params);
        
    def medianFilterStep(self, volume2assign, iterNum):
        
        iPath = 'Tilted_classification_vol%d_iter_%d.xmd' % (volume2assign, iterNum)
        oPath = 'Reference_volume%d_iter%d.vol' % (volume2assign, iterNum)

        params =  '  -i %s' % self._getExtraPath(oPath)
        params += '  -o %s' % self._getExtraPath(oPath)
        params += '  --median'

        self.runJob('xmipp_image_operate', params);
        
    
    def angularTiltAssignmentStep(self, volume2assign, iterNum):
        print iterNum
        print volume2assign
        if iterNum == 0:
            if self.inputVolume.get() is True:
                if self.TwoVolumes.get() is False:
                    print 'entro en inputVolume=True,  TwoVolume=False'
                    outputfnUntilted = self._getExtraPath('untilted_assignment_vol1_iter_0.xmd')
                    outputfnTilted = self._getExtraPath('tilted_assignment_vol1_iter_0.xmd')
                    RefVolume = self.Volume1.get().getFileName()
                else:
                    if (volume2assign == 1):
                        RefVolume = self.Volume1.get().getFileName()
                        fnbase = 'vol1' 
                    else:
                        RefVolume = self.Volume2.get().getFileName()
                        fnbase = 'vol2'
                    outputfnUntilted = self._getExtraPath(fnbase+'_untilted.xmd')
                    outputfnTilted = self._getExtraPath(fnbase+'_tilted.xmd')
            else:
                print 'TODO'
        else:
            print 'Aqui'
            fnpath1 = 'untilted_assignment_vol%d_iter_%d.xmd' % (volume2assign, iterNum)
            outputfnUntilted = self._getExtraPath(fnpath1)
            fnpath2 = 'tilted_assignment_vol%d_iter_%d.xmd' % (volume2assign, iterNum)
            outputfnTilted   = self._getExtraPath(fnpath2)
            fnpath3 = 'Reference_volume%d_iter%d.vol' %(volume2assign, iterNum-1)
            RefVolume = self._getExtraPath(fnpath3)
        
            
        params =  ' --untiltparticles %s' % self._getExtraPath() + '/all_untilted_particles_input.xmd'
        params += ' --tiltparticles %s' % self._getExtraPath() + '/all_tilted_particles_input.xmd'
        params += ' --odir %s' % self._getExtraPath()
        params += ' --untiltassignment %s' % outputfnUntilted
        params += ' --tiltassignment %s' % outputfnTilted
        params += ' --angular_sampling %f' % self.samplingAngle.get()
        params += ' --maxshift %f' % self.maxShift.get()
        params += ' --pad %f' % self.padFactor.get()
        params += ' --initvol %s' % RefVolume
        params += ' --micangles %s' % self.path_angles

        self.runJob('xmipp_angular_tilt_pair_assignment', params)        
    


    def createOutputStep(self):
        USet = self._createSetOfParticles("UntiltedVol1")
        Upath_aux = 'Untilted_classification_vol1_iter_%d.xmd', int(self.iters.get())
        Upath = self._getExtraPath(Upath_aux)
        print Upath
        readSetOfParticles(Upath, USet)
        USet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        self._defineOutputs(outputUntiltedParticles1=USet)
        self._defineSourceRelation(self.tilpairparticles.get(), USet)
         
        TSet = self._createSetOfParticles("TiltedVol1")
        Tpath_aux = 'Tilted_classification_vol1_iter_%d.xmd', int(self.iters.get())
        Tpath = self._getExtraPath(Tpath_aux)
        print Tpath
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
        Upath_aux = 'Untilted_classification_vol2_iter_%d.xmd', int(self.iters.get())
        Upath = self._getExtraPath(Upath_aux)
        print Upath
        readSetOfParticles(Upath, USet)
        USet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        self._defineOutputs(outputUntiltedParticles3=USet)
        self._defineSourceRelation(self.tilpairparticles.get(), USet)
        
        TSet = self._createSetOfParticles("TiltedVol2")
        Tpath_aux = 'Tilted_classification_vol2_iter_%d.xmd', int(self.iters.get())
        Tpath = self._getExtraPath(Tpath_aux)
        print Tpath
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

