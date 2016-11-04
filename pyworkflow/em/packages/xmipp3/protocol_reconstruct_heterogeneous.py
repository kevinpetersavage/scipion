# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
Protocol to perform high-resolution reconstructions
"""

from glob import glob
import math

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam, IntParam, EnumParam
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile
from pyworkflow.em.protocol import ProtClassify3D
from pyworkflow.em.data import Volume
from pyworkflow.em.metadata.utils import getFirstRow, getSize
from os.path import join, exists, split
from pyworkflow.em.packages.xmipp3.convert import createItemMatrix, setXmippAttributes
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
import pyworkflow.em as em

import xmipp

from protocol_reconstruct_highres import XmippProtBaseReconstructHighRes


class XmippProtReconstructHeterogeneous(ProtClassify3D, XmippProtBaseReconstructHighRes):
    """3D Reconstruction with heterogeneous datasets"""
    _label = 'significant heterogeneity'

    def __init__(self, **args):
        ProtClassify3D.__init__(self, **args)
        self.alignmentMethod = self.GLOBAL_ALIGNMENT
        self.globalMethod = self.GLOBAL_SIGNIFICANT
        self.nextLowPass = False
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParamsInput(form, multivolume=True)
        form.addParam('targetResolution', FloatParam, default=8, 
                     label='Target resolution',
                     help='Target resolution to solve for the heterogeneity')    
        form.addParam('computeDiff', BooleanParam, default=False, label="Compute the difference volumes")   
        
        form.addSection(label='Next Reference')
        form.addParam('nextSpherical', BooleanParam, label="Spherical mask?", default=True,
                      help='Apply a spherical mask of the size of the particle')
        form.addParam('nextPositivity', BooleanParam, label="Positivity?", default=True,
                      help='Remove from the next reference all negative values')
        form.addParam('nextMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass). Smooth masks are recommended.')
        form.addParam('nextReferenceScript', StringParam, label="Next reference command", default="", expertLevel=LEVEL_ADVANCED, 
                      help='A command template that is used to generate next reference. The following variables can be used ' 
                           '%(sampling)s %(dim)s %(volume)s %(iterDir)s. The command should read Spider volumes and modify the input volume.'
                           'the command should be accessible either from the PATH or provide the absolute path.\n'
                           'Examples: \n'
                           'xmipp_transform_filter -i %(volume)s --fourier low_pass 15 --sampling %(sampling)s\n' 
                           '/home/joe/myScript %(volume)s sampling=%(sampling)s dim=%(dim)s')
        form.addParam('nextRemove', BooleanParam, label="Remove reference to save space?", default=True, expertLevel=LEVEL_ADVANCED, 
                      help='Remove reference volumes once they are not needed any more.')

        form.addSection(label='Angular assignment')
        form.addParam('angularMaxShift', FloatParam, label="Max. shift (%)", default=10,
                      help='Maximum shift as a percentage of the image size')
        line=form.addLine('Tilt angle:', help='0 degrees represent top views, 90 degrees represent side views', expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMinTilt', FloatParam, label="Min.", default=0, expertLevel=LEVEL_ADVANCED)
        line.addParam('angularMaxTilt', FloatParam, label="Max.", default=90, expertLevel=LEVEL_ADVANCED)
        form.addParam('numberOfReplicates', IntParam, label="Max. Number of Replicates", default=3, 
                  expertLevel=LEVEL_ADVANCED, help="Significant alignment is allowed to replicate each image up to this number of times")

        self._defineParamsWeight(form, addJumper=False)
        self._defineParamsPostProcessing(form)
        form.addParallelSection(threads=1, mpi=8)
    
    def getNumberOfReconstructedVolumes(self):
        return len(self.inputVolumes)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def createOutput(self):
        # get last iteration
        fnIterDir=glob(self._getExtraPath("Iter*"))
        lastIter=len(fnIterDir)-1
        self.fnLastDir=self._getExtraPath("Iter%03d"%lastIter)

        fnLastImages=join(self.fnLastDir,"images.xmd")
        if not exists(fnLastImages):
            raise Exception("The file %s does not exist"%fnLastImages)
        partSet = self.inputParticles.get()
        self.Ts=self.readInfoField(self.fnLastDir,"sampling",xmipp.MDL_SAMPLINGRATE)
        self.scaleFactor=self.Ts/partSet.getSamplingRate()

        classes3D = self._createSetOfClasses3D(partSet)
        # Let use an special iterator to skip missing particles
        # discarded by classification (in case of cl2d)
        setMdIter = md.SetMdIterator(fnLastImages,
                                     sortByLabel=md.MDL_PARTICLE_ID,
                                     updateItemCallback=self._updateParticle)
        
        classes3D.classifyItems(updateItemCallback=setMdIter.updateItem,
                             updateClassCallback=self._updateClass)
        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles, classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(self.Ts)
        
        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)
        
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles, volumes)

    def _updateParticle(self, particle, row):
        particle.setClassId(row.getValue(md.MDL_REF3D))
        row.setValue(xmipp.MDL_SHIFT_X,row.getValue(xmipp.MDL_SHIFT_X)*self.scaleFactor)
        row.setValue(xmipp.MDL_SHIFT_Y,row.getValue(xmipp.MDL_SHIFT_Y)*self.scaleFactor)
        setXmippAttributes(particle, row, xmipp.MDL_SHIFT_X, 
                           xmipp.MDL_SHIFT_Y, xmipp.MDL_ANGLE_ROT, 
                           xmipp.MDL_ANGLE_TILT, xmipp.MDL_ANGLE_PSI, 
                           xmipp.MDL_MAXCC, xmipp.MDL_WEIGHT)
        createItemMatrix(particle, row, align=em.ALIGN_PROJ)

    def _updateClass(self, item):
        classId = item.getObjId()
        item.setAlignment3D()
        item.setSamplingRate(self.Ts)
        item.getRepresentative().setFileName(join(self.fnLastDir,"volume%02d.vol"%classId))
        
    def doIteration000(self, inputVolumesId):
        fnDirCurrent=self._getExtraPath('Iter000')
        makePath(fnDirCurrent)
        
        # Get volume sampling rate
        TsCurrent=max(self.TsOrig,self.targetResolution.get()/3)
        self.writeInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)

        # Copy reference volumes and window if necessary
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=long(round(Xdim*self.TsOrig/TsCurrent))
        self.writeInfoField(fnDirCurrent,"size",xmipp.MDL_XSIZE,newXdim)
        
        img = ImageHandler()
        idx = 1
        for i, vol in enumerate(self.inputVolumes):
            fnVol=join(fnDirCurrent,"volume%02d.vol"%idx)
            img.convert(vol.get(), fnVol)
            if newXdim!=vol.get().getDim()[0]:
                self.runJob('xmipp_image_resize',"-i %s --dim %d"%(fnVol,newXdim),numberOfMpi=1)
            idx+=1

        self.evaluateReconstructions(0)
        self.prepareImages(0)

    def evaluateReconstructions(self,iteration):
        self.alignReconstructedVolumes(iteration)
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        cleanPath(join(fnDirCurrent,"volumeAvg.mrc"))
     
        if self.postAdHocMask.hasValue():
            fnMask=join(fnDirCurrent,"mask.vol")
            if not exists(fnMask):
                TsCurrent = self.readInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE)
                volXdim = self.readInfoField(fnDirCurrent, "size", xmipp.MDL_XSIZE)
                self.prepareMask(self.postAdHocMask.get(), fnMask, TsCurrent, volXdim)
            for i in range(1,self.getNumberOfReconstructedVolumes()+1):
                fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
                self.runJob("xmipp_image_operate","-i %s --mult %s"%(fnVol,fnMask),numberOfMpi=1)
            cleanPath(fnMask)
        
        if self.computeDiff:
            for i in range(1,self.getNumberOfReconstructedVolumes()):
                fnVoli=join(fnDirCurrent,"volume%02d.vol"%i)
                fnVoliAdjusted=join(fnDirCurrent,"volume%02d_adjusted.vol"%i)
                for j in range(2,self.getNumberOfReconstructedVolumes()+1):
                    fnVolj=join(fnDirCurrent,"volume%02d.vol"%j)
                    fnDiff=join(fnDirCurrent,"diff%02d_%02d.vol"%(j,i))

                    copyFile(fnVoli, fnVoliAdjusted)
                    self.runJob("xmipp_volume_align","--i1 %s --i2 %s --least_squares --apply %s"%(fnVolj, fnVoliAdjusted, fnVoliAdjusted),numberOfMpi=1)
                    self.runJob("xmipp_image_operate","-i %s --minus %s -o %s"%(fnVolj, fnVoliAdjusted, fnDiff),numberOfMpi=1)
                    cleanPath(fnVoliAdjusted)
        
        if iteration>1:
            fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
            for i in range(0,self.getNumberOfReconstructedVolumes()):
                fnCurrent=join(fnDirCurrent,"images%02d.xmd"%(i+1))
                fnPrevious=join(fnDirPrevious,"images%02d.xmd"%(i+1))
                fnIntersection=self._getExtraPath("intersection.xmd")
                fnUnion=self._getExtraPath("union.xmd")
                self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s itemId -o %s"%(fnCurrent,fnPrevious,fnIntersection),numberOfMpi=1)
                self.runJob("xmipp_metadata_utilities","-i %s --set union %s itemId -o %s"%(fnCurrent,fnPrevious,fnUnion),numberOfMpi=1)
                
                sizeIntersection = float(md.getSize(fnIntersection))
                sizeUnion = float(md.getSize(fnUnion))
                
                print("Stability of class %d: %f"%(i+1,sizeIntersection/sizeUnion))
                cleanPath(fnIntersection)
                cleanPath(fnUnion)                

    def prepareImages(self,iteration):
        fnDir=self._getPath()
        TsCurrent=self.readInfoField(self._getExtraPath("Iter000"),"sampling",xmipp.MDL_SAMPLINGRATE)

        print "Preparing images to sampling rate=",TsCurrent
        Xdim=self.inputParticles.get().getDimensions()[0]
        newXdim=self.readInfoField(self._getExtraPath("Iter000"),"size",xmipp.MDL_XSIZE)
        
        fnNewParticles=join(fnDir,"images.stk")
        fnImgs = join(fnDir,"images.xmd")
        if newXdim!=Xdim:
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(self.imgsFn,fnNewParticles,newXdim),
                        numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        else:
            self.runJob("xmipp_image_convert","-i %s -o %s --save_metadata_stack %s"%(self.imgsFn,fnNewParticles,fnImgs),numberOfMpi=1)
        R=self.particleRadius.get()
        if R<=0:
            R=self.inputParticles.get().getDimensions()[0]/2
        R=min(round(R*self.TsOrig/TsCurrent*(1+self.angularMaxShift.get()*0.01)),newXdim/2)
        self.runJob("xmipp_transform_mask","-i %s --mask circular -%d"%(fnNewParticles,R),numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())

        row=getFirstRow(fnImgs)
        if row.containsLabel(xmipp.MDL_CTF_MODEL) or row.containsLabel(xmipp.MDL_CTF_DEFOCUSU):
            args="-i %s --sampling_rate %f --correct_envelope"%(fnImgs,TsCurrent)
            if self.inputParticles.get().isPhaseFlipped():
                args+=" --phase_flipped"
            self.runJob("xmipp_ctf_correct_wiener2d",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        
    def globalAssignment(self,iteration):
        fnDirPrevious=self._getExtraPath("Iter%03d"%(iteration-1))
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        makePath(fnDirCurrent)

        fnGlobal=join(fnDirCurrent,"globalAssignment")
        makePath(fnGlobal)

        TsCurrent=self.readInfoField(self._getExtraPath("Iter000"),"sampling",xmipp.MDL_SAMPLINGRATE)
        newXdim=self.readInfoField(self._getExtraPath("Iter000"),"size",xmipp.MDL_XSIZE)
        self.writeInfoField(fnGlobal,"size",xmipp.MDL_XSIZE,newXdim)
        self.writeInfoField(fnDirCurrent,"size",xmipp.MDL_XSIZE,newXdim)
        self.writeInfoField(fnDirCurrent,"sampling",xmipp.MDL_SAMPLINGRATE,TsCurrent)
        self.prepareReferences(fnDirPrevious,fnGlobal,TsCurrent,self.targetResolution.get())

        # Calculate angular step at this resolution
        ResolutionAlignment=self.targetResolution.get()
        newXdim=self.readInfoField(self._getExtraPath("Iter000"),"size",xmipp.MDL_XSIZE)
        angleStep=self.calculateAngStep(newXdim, TsCurrent, ResolutionAlignment)
        angleStep=max(angleStep,5.0)
        self.writeInfoField(fnGlobal,"angleStep",xmipp.MDL_ANGLE_DIFF,float(angleStep))
        
        # Generate projections
        fnImgs = self._getPath("images.xmd")
        fnGalleryAll = join(fnGlobal,"gallery.xmd")
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnReferenceVol=join(fnGlobal,"volumeRef%02d.vol"%i)
            fnGallery=join(fnGlobal,"gallery%02d.stk"%i)
            fnGalleryXmd=join(fnGlobal,"gallery%02d.doc"%i)
            args="-i %s -o %s --sampling_rate %f --perturb %f --sym %s --min_tilt_angle %f --max_tilt_angle %f"%\
                 (fnReferenceVol,fnGallery,angleStep,math.sin(angleStep*math.pi/180.0)/4,self.symmetryGroup,self.angularMinTilt.get(),self.angularMaxTilt.get())
            args+=" --compute_neighbors --angular_distance -1 --experimental_images %s"%fnImgs
            self.runJob("xmipp_angular_project_library",args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
            self.runJob("xmipp_metadata_utilities","-i %s -o gallery%02d@%s --mode append"%(fnGalleryXmd,i,fnGalleryAll),numberOfMpi=1)
            cleanPath(join(fnGlobal,"gallery%02d_angles.doc"%i))
            cleanPath(join(fnGlobal,"gallery%02d_sampling.xmd"%i))
            cleanPath(fnGalleryXmd)

        maxShift=round(self.angularMaxShift.get()*newXdim/100)
        args='-i %s --initgallery %s --maxShift %d --odir %s --dontReconstruct --useForValidation %d'%\
             (fnImgs,fnGalleryAll,maxShift,fnGlobal,self.numberOfReplicates.get()-1)
        self.runJob('xmipp_reconstruct_significant',args,numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            fnSignificant=join(fnGlobal,"angles_iter001_%02d.xmd"%(i-1))
            moveFile(fnSignificant, fnAngles)
            self.runJob('xmipp_metadata_utilities',"-i %s --fill ref3d constant %d"%(fnAngles,i),numberOfMpi=1)
            
            fnImages=join(fnDirCurrent,"images%02d.xmd"%i)
            fnSignificant=join(fnGlobal,"images_iter001_%02d.xmd"%(i-1))
            moveFile(fnSignificant, fnImages)
            self.runJob('xmipp_metadata_utilities',"-i %s --fill ref3d constant %d"%(fnImages,i),numberOfMpi=1)

            cleanPath(join(fnGlobal,"images_significant_iter001_%02d.xmd"%(i-1)))
        if self.saveSpace:
            self.runJob("rm -f",fnGlobal+"/gallery*",numberOfMpi=1)
            
    def weightParticles(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        fnGeneralAngles=join(fnDirCurrent,"angles.xmd")
        fnGeneralImages=join(fnDirCurrent,"images.xmd")
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            fnImages=join(fnDirCurrent,"images%02d.xmd"%i)
                
            if self.weightSSNR:
                row=getFirstRow(fnAngles)
                if row.containsLabel(xmipp.MDL_WEIGHT_SSNR):
                    self.runJob("xmipp_metadata_utilities","-i %s --operate drop_column weightSSNR"%fnAngles,numberOfMpi=1)
                self.runJob("xmipp_metadata_utilities","-i %s --set join %s particleId"%\
                            (fnAngles,self._getExtraPath("ssnrWeights.xmd")),numberOfMpi=1)

            mdAngles=xmipp.MetaData(fnAngles)
            for objId in mdAngles:
                weight=mdAngles.getValue(xmipp.MDL_WEIGHT_SIGNIFICANT,objId)
                if self.weightSSNR:
                    aux=mdAngles.getValue(xmipp.MDL_WEIGHT_SSNR,objId)
                    weight*=aux
                mdAngles.setValue(xmipp.MDL_WEIGHT,weight,objId)
            mdAngles.write(fnAngles)
            
            if i==1:
                copyFile(fnAngles,fnGeneralAngles)
                copyFile(fnImages,fnGeneralImages)
            else:
                self.runJob('xmipp_metadata_utilities',"-i %s --set union_all %s"%(fnGeneralAngles,fnAngles),numberOfMpi=1)
                self.runJob('xmipp_metadata_utilities',"-i %s --set union_all %s"%(fnGeneralImages,fnImages),numberOfMpi=1)

    def reconstruct(self, iteration):
        fnDirCurrent=self._getExtraPath("Iter%03d"%iteration)
        for i in range(1,self.getNumberOfReconstructedVolumes()+1):
            fnAngles=join(fnDirCurrent,"angles%02d.xmd"%i)
            fnVol=join(fnDirCurrent,"volume%02d.vol"%i)
            if not exists(fnVol):
                args="-i %s -o %s --sym %s --weight --thr %d"%(fnAngles,fnVol,self.symmetryGroup,self.numberOfThreads.get())
                self.runJob("xmipp_reconstruct_fourier",args,numberOfMpi=self.numberOfMpi.get()+1)
    
    def cleanDirectory(self, iteration):
        pass