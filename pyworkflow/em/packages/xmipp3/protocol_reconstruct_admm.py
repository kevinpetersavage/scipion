# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.metadata.utils import getFirstRow
import xmipp

class XmippProtReconstructADMM(ProtReconstruct3D):
    """    
    Reconstruct a volume using Xmipp_reconstruct_admm from a given set of particles.
    The alignment parameters will be converted to a Xmipp xmd file
    and used as direction projections to reconstruct.
    """
    _label = 'reconstruct admm'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('inputMask', PointerParam, label="Mask", pointerClass='VolumeMask', allowsNull=True,
                      help='The mask values must be between 0 (remove these pixels) and 1 (let them pass).')
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).\n'
                           'Param *--maxres* in Xmipp.') 
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")

        form.addParam('extraParams', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *xmipp_reconstruct_fourier* program:\n'
                      """
                      --iter () : Subtract projections of this map from the images used for reconstruction
                      """)

        form.addParallelSection(threads=1, mpi=4)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd'),
            'output_volume': self._getPath('output_volume.vol')
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        #self._insertFunctionStep('correctCTFStep')
        self._insertFunctionStep('reconstructHtbStep')
        self._insertFunctionStep('reconstructHtHStep')
        self._insertFunctionStep('reconstructAdmmStep')
#         self._insertFunctionStep('createOutputStep')
        
        
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        fnParticles = self._getFileName('input_xmd')
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, fnParticles)
        row=getFirstRow(fnParticles)
        if row.containsLabel(xmipp.MDL_CTF_DEFOCUSU) and not row.containsLabel(xmipp.MDL_CTF_K):
            self.runJob("xmipp_metadata_utilities",'-i %s --fill ctfK constant 1'%fnParticles,numberOfMpi=1)
            
        fnAligned = self._getExtraPath("imagesAligned.xmd")
        fnAlignedStk = self._getExtraPath("imagesAligned.stk")
        self.runJob("xmipp_metadata_utilities",'-i %s -o %s --operate keep_column "itemId image shiftX shiftY"'%(fnParticles,fnAligned),numberOfMpi=1)
        self.runJob("xmipp_transform_geometry",'-i %s -o %s --apply_transform --save_metadata_stack %s'%(fnAligned,fnAlignedStk,fnAligned),numberOfMpi=1)
        
        import pyworkflow.em.metadata as metadata
        from xmipp import SymList
        mdOut = metadata.MetaData()
        self.iterMd = metadata.iterRows(fnParticles, sortByLabel=metadata.MDL_ITEM_ID)
        symlist=SymList()
        symlist.readSymmetryFile(self.symmetryGroup.get())
        for row in self.iterMd:
            itemId = row.getValue(metadata.MDL_ITEM_ID)
            rot = row.getValue(metadata.MDL_ANGLE_ROT)
            tilt = row.getValue(metadata.MDL_ANGLE_TILT)
            psi = row.getValue(metadata.MDL_ANGLE_PSI)
            
            equivalentAngles = symlist.symmetricAngles(rot,tilt,psi)
            for angles in equivalentAngles:
                rowOut = metadata.Row()
                rowOut.setValue(metadata.MDL_ITEM_ID,itemId)
                rowOut.setValue(metadata.MDL_ANGLE_ROT,angles[0])
                rowOut.setValue(metadata.MDL_ANGLE_TILT,angles[1])
                rowOut.setValue(metadata.MDL_ANGLE_PSI,angles[2])
                rowOut.addToMd(mdOut)
        fnAllAngles = self._getExtraPath("all_angles.xmd")
        mdOut.write(fnAllAngles)
        self.runJob("xmipp_metadata_utilities",'-i %s --set join %s itemId'%(fnAligned,fnAllAngles),numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities",'-i %s --set join %s itemId'%(fnAligned,fnParticles),numberOfMpi=1)
        cleanPath(fnAllAngles)
        cleanPath(fnParticles)

    def getDigRes(self):
        maxRes = self.maxRes.get()
        if maxRes == -1:
            digRes = 0.5
        else:
            digRes = self.inputParticles.get().getSamplingRate() / self.maxRes.get()
        return digRes

    def reconstructHtbStep(self):
        """ Create the input file in STAR format as expected by Xmipp.
        If the input particles comes from Xmipp, just link the file. 
        """
        params =  '  -i %s' % self._getExtraPath("imagesAligned.xmd")
        params += '  -o %s' % self._getExtraPath('Htb.vol')
        params += ' --max_resolution %0.3f' %self.getDigRes()
        params += ' --padding %0.3f' % self.pad.get()
        params += ' --thr %d' % self.numberOfThreads.get()
        params += ' --sampling %f' % self.inputParticles.get().getSamplingRate()
        params += ' --iter 0'
        self.runJob('xmipp_reconstruct_fourier', params)
            
    def correctCTFStep(self):
        params =  '  -i %s --sampling_rate %f --correct_envelope' % (self._getExtraPath("imagesAligned.xmd"),self.inputParticles.get().getSamplingRate())
        if self.inputParticles.get().isPhaseFlipped():
            params += " --phase_flipped"
        self.runJob('xmipp_ctf_correct_wiener2d', params)

    def reconstructHtHStep(self):
        I = xmipp.Image()
        I.setDataType(xmipp.DT_DOUBLE)
        I.initLPFImpulseResponse(2*self.inputParticles.get().getDimensions()[0]-1,self.getDigRes())
        fnKernel = self._getExtraPath("kernel.xmp")
        I.write(fnKernel)

        fnKernelXmd = self._getExtraPath("kernel.xmd")
        params =  '  -i %s -o %s --fill image constant %s' % (self._getExtraPath("imagesAligned.xmd"),fnKernelXmd,fnKernel)
        self.runJob('xmipp_metadata_utilities', params, numberOfMpi=1)
        
        params =  '  -i %s' % fnKernelXmd
        params += '  -o %s' % self._getExtraPath('HtH.vol')
        params += ' --max_resolution %0.3f' %self.getDigRes()
        params += ' --padding %0.3f' % self.pad.get()
        params += ' --thr %d' % self.numberOfThreads.get()
        params += ' --sampling %f' % self.inputParticles.get().getSamplingRate()
        params += ' --iter 0'
        self.runJob('xmipp_reconstruct_fourier', params)
        
    def reconstructAdmmStep(self):
        """ Create the input file in STAR format as expected by Xmipp.
        If the input particles comes from Xmipp, just link the file. 
        """
        params =  ' -i %s' % self._getExtraPath("imagesAligned.xmd")
        params += ' --oroot %s' % self._getExtraPath('admm')
        params += ' --Htb %s'%self._getExtraPath('Htb.vol')
        params += ' --HtKH %s'%self._getExtraPath('HtH.vol')
        params += ' --positivity'
        #params += ' --firstVolume %s' %self._getExtraPath('Htb.vol')
        params += ' --cgiter 2'
        params += ' --admmiter 30'
        params += ' --lambda 1e4'
        params += ' --mu 1e4'
        self.runJob('xmipp_reconstruct_admm2', params,numberOfMpi=1)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputParticles, volume)

    def _validate(self):
        errors = []
        if not self.inputParticles.get().isPhaseFlipped():
            errors.append("This protocol requires the images to be phase flipped. Re-extract the particles with phase flipping")
        return errors