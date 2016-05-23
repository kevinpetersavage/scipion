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

from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, BooleanParam, StringParam, STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.em.metadata import MetaData, getBlocksInMetaDataFile
from convert import readSetOfVolumes, getImageLocation
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from pyworkflow.em.metadata.constants import (MDL_ANGLE_Y, MDL_IMAGE, MDL_ANGLE_Y2, MDL_ANGLE_TILT, MDL_PARTICLE_ID, 
                                              MD_APPEND, MDL_MICROGRAPH_ID)
from pyworkflow.em.data import Volume
from shutil import copyfile
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em as em
from pyworkflow.em.packages.xmipp3.convert import readSetOfParticles



class XmippProtAngularTiltPairAssignment(ProtRefine3D):
    """    
    Given two sets of particles (untilted and tilted) the protocol assings angles to the 
    untilted particles and then determines by means of the tilt axis and tilt angle the
     orientation of the tilted particles.
    """
    _label = 'Angular tiltpair'
    
    def __init__(self, *args, **kwargs):
        ProtRefine3D.__init__(self, *args, **kwargs)
        #self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
       
        form.addParam('tilpairparticles', PointerParam, pointerClass='ParticlesTiltPair', 
                      label="Set of particles tilt pairs", 
                      help='Select a set of particles tilt pairs.')

        form.addParam('samplingangle', FloatParam, default=5,
                      label="Angular sampling",  
                      help='The volume is projected, in many directions.'
                      'The angular sampling defines the distance between two '
                      'adjacents directions')   
        
        form.addParam('thereisRefVolume', BooleanParam, default=False,
                      label="Is there a reference volume(s)?", 
                      help='You may use a reference volume to initialize the calculations.)')
        
        form.addParam('initvolume', PointerParam, pointerClass='Volume', allowsNull = True, condition="thereisRefVolume",
                      label="Initial Volumen",   
                      help='Select an initial volumen for assigning angles.')
        
        form.addParam('Resizing', BooleanParam, default=True,
                      label="Should the images be resized?", 
                      help='Resize images for speeding up the processing.)')
        
        form.addParam('maxshift', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Maximum shift",  
                      help='Maximum shift for aligning images. If the obtained shift in the alignment'
                      'is higher than the maximum shift, the align between those images'
                      'will not be considered. If maxshift = -1, '
                      'there will not be considered a maximum shift')
        
        form.addParam('padfactor', FloatParam, default=2, expertLevel=LEVEL_ADVANCED,
                      label="Padding factor",  
                      help='Padding factor for controlling Fourier accuracy'
                      'by default = 2')
        
        form.addParam('iterations', FloatParam, default=3, expertLevel=LEVEL_ADVANCED,
                      label="Number of Iterations",  
                      help='The volume is obtained by means of an iterative process.'
                      'The number of iterations can be determined'
                      'by default = 2')
        
        form.addParam('symmetry', StringParam, default='c1', expertLevel=LEVEL_ADVANCED,
                      label="Symmetry",  
                      help='Symmetry group. By default = c1')
        
#         form.addParam('downsamplingfactor', FloatParam, default=0.4, expertLevel=LEVEL_ADVANCED,
#                       label="Downsampling factor",  
#                       help='The images are downsampled in order to speed up the process'
#                       'by default = 0.4')
        
        form.addParam('maxFreq', IntParam, default=12, expertLevel=LEVEL_ADVANCED,
                      label='Max frequency of the initial volume',
                      help=' Max frequency of the initial volume in Angstroms')
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------


    def _insertAllSteps(self):        
        self.micsFn = self._getPath()
        
        # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep')
        deps = []
                 
        fnUntilt = self._getExtraPath('all_untilted_particles_input.xmd')
        fnTilt = self._getExtraPath('all_tilted_particles_input.xmd')
        fnmic = self._getExtraPath('Mic_angles.xmd')
        
        self.Xdim = self.tilpairparticles.get().getUntilted().getDimensions()[0] 
        if self.Resizing.get() == True:
            filteresId = self._insertFunctionStep('filterResizeStep', self.Xdim, fnUntilt, fnTilt, prerequisites=[convertId])
            fnUntilt = splitext(fnUntilt)[0]+'_filtered_resized.xmd'
            fnTilt = splitext(fnTilt)[0]+'_filtered_resized.xmd'
            maskId = self._insertFunctionStep('mask2InputParticles', fnUntilt, fnTilt, prerequisites=[filteresId])
        else:
            maskId = self._insertFunctionStep('mask2InputParticles', fnUntilt, fnTilt, prerequisites=[convertId])
        
        print '--------------------------------------------'
        print '--------------------------------------------'
        print '--------------------------------------------'
        print '--------------------------------------------'
        print fnUntilt
        print fnTilt
#         #TODO: resize output volume 
        for k in range(1, int(self.iterations.get())+1):
            self._insertFunctionStep('angularTiltAssignmentStep', fnUntilt, fnTilt, fnmic, k, prerequisites=[maskId])
            self._insertFunctionStep('ReconstructFourierStep', k)
# 
# 
        self._insertFunctionStep('createOutputStep')
        
    def filterResizeStep(self, Xdim, fnUntilt, fnTilt):
        ##
        fnUntilt_Noext = splitext(fnUntilt)[0]
        fnTilt_Noext = splitext(fnTilt)[0]
        
        fnUntilt_fil_stk = fnUntilt_Noext+'_filtered.stk'
        fnTilt_fil_stk = fnTilt_Noext+'_filtered.stk'
        ##
        
        fnUntilt_ReducedNoExt = splitext(fnUntilt_fil_stk)[0]
        fnTilt_ReducedNoExt = splitext(fnTilt_fil_stk)[0]
        
        maxFreq = self.maxFreq.get()
        inputSet = self.tilpairparticles.get().getUntilted()
        ts = inputSet.getSamplingRate()
        K = 0.25 * (maxFreq / ts)
        if K < 1:
            K = 1
        self.Xdim2 = Xdim / K
        if self.Xdim2 < 32:
            self.Xdim2 = 32
            K = self.Xdim / self.Xdim2
        freq = ts / maxFreq
        ts = K * ts

        paramsfilterU  =  ' -i %s' % fnUntilt
        paramsfilterU +=  ' -o %s'% fnUntilt_fil_stk
        paramsfilterU +=  ' --fourier low_pass %f' % freq
        paramsfilterU +=  ' --keep_input_columns'
        paramsfilterU +=  ' --save_metadata_stack'
        
        paramsfilterT  =  ' -i %s' % fnTilt
        paramsfilterT +=  ' -o %s'% fnTilt_fil_stk
        paramsfilterT +=  ' --fourier low_pass %f' % freq
        paramsfilterT +=  ' --keep_input_columns'
        paramsfilterT +=  ' --save_metadata_stack'
        
        paramsResizeU  =  ' -i %s' % fnUntilt_ReducedNoExt+'.xmd'
        paramsResizeU +=  ' --fourier %d' % self.Xdim2
        paramsResizeU +=  ' -o %s' % fnUntilt_ReducedNoExt+'_resized.stk'
        paramsResizeU +=  ' --keep_input_columns'
        paramsResizeU +=  ' --save_metadata_stack'
         
        paramsResizeT  =  ' -i %s' % fnTilt_ReducedNoExt+'.xmd'
        paramsResizeT +=  ' --fourier %d' % self.Xdim2
        paramsResizeT +=  ' -o %s' % fnTilt_ReducedNoExt+'_resized.stk'
        paramsResizeT +=  ' --keep_input_columns'
        paramsResizeT +=  ' --save_metadata_stack'
         
        paramsResizeVol  =  ' -i %s' % self._getExtraPath('init_vol.vol')
        paramsResizeVol +=  ' --fourier %d' % self.Xdim2
        paramsResizeVol +=  ' -o %s' % self._getExtraPath('init_vol.vol')
        paramsResizeVol +=  ' --keep_input_columns'
        paramsResizeVol +=  ' --save_metadata_stack'
        
        self.runJob('xmipp_transform_filter',paramsfilterU)
        self.runJob('xmipp_transform_filter',paramsfilterT)
        self.runJob('xmipp_image_resize',paramsResizeU)
        self.runJob('xmipp_image_resize',paramsResizeT)
        if self.thereisRefVolume.get() == True:
            self.runJob('xmipp_image_resize',paramsResizeVol)

        

    def mask2InputParticles(self, fnUntilt, fnTilt):
        if self.Resizing.get() == True:
            fnUntiltStk = splitext(fnUntilt)[0]+'.stk'
            xdim, _, _, _ = em.ImageHandler().getDimensions(fnUntiltStk)
        else:
            xdim = self.tilpairparticles.get().getUntilted().getDimensions()[0]
        
        maskArgs = "-i %s --mask circular %d -v 0" % (fnUntilt, -xdim/2)
        self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)
        maskArgs = "-i %s --mask circular %d -v 0" % (fnTilt, -xdim/2)
        self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)
        maskArgs = "-i %s --mask circular %d -v 0" % (self._getExtraPath('init_vol.vol'), -xdim/2)
        self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)


    def convertInputStep(self):
        """ Read the input metadatata.
        """
        # Get the converted input micrographs in Xmipp format
        lastMicId = -1
        path_u = self._getExtraPath() + '/untilted_particles_in_blocks_input.xmd'
        path_t = self._getExtraPath() + '/tilted_particles_in_blocks_input.xmd'
        path_u_all = self._getExtraPath() + '/all_untilted_particles_input.xmd'
        path_t_all = self._getExtraPath() + '/all_tilted_particles_input.xmd'
        path_angles = self._getExtraPath() + '/Mic_angles.xmd'
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
        

        md_angle.write(path_angles)
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
        
        if self.thereisRefVolume.get() is True:
            volume_init = self.initvolume.get().getFileName()
            copyfile(volume_init,self._getExtraPath('init_vol.vol'))
        
   
    def angularTiltAssignmentStep(self, fnUntilt, fnTilt, fnmic, iter_num):
        

        if iter_num is 1:
            if self.thereisRefVolume.get() is True:
                volume_init = self._getExtraPath('init_vol.vol')
#                 xdim = self.tilpairparticles.get().getUntilted().getDimensions()[0]
#                 maskArgs = "-i %s --mask circular %d -v 0" % (volume_init, -xdim/2)
#                 self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)
            else:
                volume_init = ''
        else:
            volume_init = self.getIterVolume(iter_num-1)
            
        params =  ' --untiltparticles %s' % fnUntilt
        params += ' --tiltparticles %s' % fnTilt
        params += ' --odir %s' % self._getExtraPath()
        params += ' --sym %s' % self.symmetry.get()
        params += ' --angular_sampling %f' % self.samplingangle.get()
        params += ' --maxshift %f' % self.maxshift.get()
        params += ' --pad %f' % self.padfactor.get()
        params += ' --initvol %s' % volume_init
        params += ' --micangles %s' % fnmic

        self.runJob('xmipp_angular_tilt_pair_assignment', params)
    
    
    def ReconstructFourierStep(self, iterNumber):
        
        input_particles = 'particles@%s/All_angular_assignment_lastiter.xmd' % self._getExtraPath()
        params =  '  -i %s' % input_particles
        params += '  -o %s' % self.getIterVolume(iterNumber)
        params += '  --sym %s' % self.symmetry.get()
        params += '  --max_resolution %f' % 0.5
        params += '  --padding %f' %self.padfactor.get()
        
        self.runJob('xmipp_reconstruct_fourier', params);
        
        xdim = self.tilpairparticles.get().getUntilted().getDimensions()[0]
        maskArgs = "-i %s --mask circular %d -v 0" % (self.getIterVolume(iterNumber), -xdim/2)
        self.runJob('xmipp_transform_mask', maskArgs, numberOfMpi=1)

    
    def createOutputStep(self):
        #Output Volume
        lastIter = self.iterations.get()
        
        volume_path = self.getIterVolume(lastIter)
        
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        
        readSetOfVolumes(volume_path, volumesSet)
        
        self._defineOutputs(outputVolume=volumesSet)
        self._defineSourceRelation(self.tilpairparticles, volumesSet)
        
        #################################################################
        USet = self._createSetOfParticles("Untilted")
        Upath = self._getExtraPath('all_untilted_particles_input_angular_assignment_lastiter.xmd')
        readSetOfParticles(Upath, USet)
        USet.setSamplingRate((self.tilpairparticles.get().getUntilted().getSamplingRate()))
        self._defineOutputs(outputParticles=USet)
        self._defineSourceRelation(self.tilpairparticles.get().getUntilted(), USet)
       
               
        TSet = self._createSetOfParticles("Tilted")
        Tpath = self._getExtraPath('all_tilted_particles_input_angular_assignment_lastiter.xmd')
        readSetOfParticles(Tpath, TSet)
        TSet.setSamplingRate((self.tilpairparticles.get().getTilted().getSamplingRate()))
        self._defineOutputs(outputParticles=TSet)
        self._defineSourceRelation(self.tilpairparticles.get().getTilted(), TSet)


    def getIterVolume(self, iterNumber):
        return self._getExtraPath('volume_iter_%03d.vol' % iterNumber)
        
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