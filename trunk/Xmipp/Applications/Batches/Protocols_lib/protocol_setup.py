#!/usr/bin/env python
#------------------------------------------------------------------------
# General script to setup standard Xmipp protocols
#  - Preprocessing of micrographs
#  - Manual particle picking
#  - Preprocessing of extracted particles
#  - 2D image alignment and classification (by ml2d & kerdenSOM)
#  - 2D image analysis by classification of rotational spectra
#  - 3D classification by ml3D
#  - 3D refinement by standard projection matching
#  - 3D high resolution reconstruction
#
# Example use:
# ./setup_protocols.py
#
# Author: Sjors Scheres, March 2007
#------------------------------------------------------------------------
# Choose the protocol(s) you want to setup:
#------------------------------------------------------------------------
# {setup-pre} Preprocess micrographs
SetupPreProcessMicrographs=False
# {setup-pre} Manual particle selection
SetupParticlePick=False
# {setup-pre} Preprocess particles
SetupPreProcessParticles=False
# {setup-2d} ML2D classification
SetupML2D=False
# {setup-2d} kerdenSOM classification 
SetupKerdensom=False
# {setup-2d} Rotational spectra classification
SetupRotSpectra=False
# {setup-3d} Random Conical Tilt
SetupRCT=False
# {setup-3d} ML3D classification
SetupML3D=False
# {setup-3d} Projection matching refinement
SetupProjMatch=False
# {setup-3d} High Resolution reconstruction
SetupHighRes3D=False
#------------------------------------------------------------------------
# {section} Global Parameters
#------------------------------------------------------------------------
# Root directory name for this project:
ProjectDir="/home/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
# {expert} Directory name for preprocessing:
PreProcessDir="Preprocessing"
# {expert} Directory name for particle images:
ImagesDir="Images"
# {expert} Directory name for ML2D classification:
ML2DDir="ML2D"
# {expert} Directory name for rotational spectra classification:
RotSpectraDir="Rotspectra"
# {expert} Directory name for Random Conical Tilt reconstruction:
RCTDir="RCT"
# {expert} Directory name for ML3D classification:
ML3DDir="ML3D"
# {expert} Directory name for projection matching refinement:
ProjMatchDir="ProjMatching"
# {expert} Directory name for HighRes3D:
HighRes3DDir="HighRes3D"
#
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class setup_protocols_class:
       #init variables
        def __init__(self,
                     SetupPreProcessMicrographs,
                     SetupParticlePick,
                     SetupPreProcessParticles,
                     SetupML2D,
                     SetupKerdensom,
                     SetupRotSpectra,
                     SetupRCT,
                     SetupML3D,
                     SetupProjMatch,
                     SetupHighRes3D,
                     ProjectDir,
                     LogDir,
                     ImagesDir,
                     PreProcessDir,
                     ML2DDir,
                     RotSpectraDir,
                     RCTDir,
                     ML3DDir,
                     ProjMatchDir,
                     HighRes3DDir,
                     AutoLaunch):

            import os,sys

            self.SetupPreProcessMicrographs=SetupPreProcessMicrographs
            self.SetupParticlePick=SetupParticlePick
            self.SetupPreProcessParticles=SetupPreProcessParticles
            self.SetupML2D=SetupML2D
            self.SetupKerdensom=SetupKerdensom
            self.SetupRotSpectra=SetupRotSpectra
            self.SetupRCT=SetupRCT
            self.SetupML3D=SetupML3D
            self.SetupProjMatch=SetupProjMatch
            self.SetupHighRes3D=SetupHighRes3D

            self.ProjectDir=ProjectDir
            self.LogDir=LogDir
            self.ImagesDir=ImagesDir
            self.PreProcessDir=PreProcessDir
            self.ML2DDir=ML2DDir
            self.RotSpectraDir=RotSpectraDir
            self.RCTDir=RCTDir
            self.ML3DDir=ML3DDir
            self.ProjMatchDir=ProjMatchDir
            self.HighRes3DDir=HighRes3DDir

            self.AutoLaunch=AutoLaunch

            scriptdir=os.path.expanduser('~')+'/scripts/'
            sys.path.append(scriptdir) # add default search path
            self.SYSTEMSCRIPTDIR=scriptdir

            # Which scripts and which directories to use
            self.library={}
            self.library['SetupPreProcessMicrographs']=[self.SetupPreProcessMicrographs,
                                                self.PreProcessDir,
                                                ['protocol_preprocess_micrographs.py']]
            self.library['SetupParticlePick']=[self.SetupParticlePick,
                                                self.PreProcessDir,
                                                ['protocol_particle_pick.py']]
            self.library['SetupPreProcessParticles']=[self.SetupPreProcessParticles,
                                                self.PreProcessDir,
                                                ['protocol_preprocess_particles.py']]
            self.library['SetupML2D']=[self.SetupML2D,
                                         self.ML2DDir,
                                         ['protocol_ML2D.py','visualize_ML2D.py']]
            self.library['SetupKerdensom']=[self.SetupKerdensom,
                                              self.ML2DDir,
                                              ['protocol_kerdensom.py','visualize_kerdensom.py']]
            self.library['SetupRotSpectra']=[self.SetupRotSpectra,
                                               self.RotSpectraDir,
                                               ['protocol_rotspectra.py','visualize_rotspectra.py']]
            self.library['SetupRCT']=[self.SetupRCT,
                                        self.RCTDir,
                                        ['protocol_rct.py','visualize_rct.py']]
            self.library['SetupML3D']=[self.SetupML3D,
                                        self.ML3DDir,
                                        ['protocol_ML3D.py','visualize_ML3D.py']]
            self.library['SetupProjMatch']=[self.SetupProjMatch,
                                        self.ProjMatchDir,
                                        ['protocol_projmatch.py','visualize_projmatch.py']]
            self.library['SetupHighRes3D']=[self.SetupHighRes3D,
                                        self.HighRes3DDir,
                                        ['protocol_highres3D.py','visualize_highres3D.py']]

            # For automated editing of default directories in protocols
            self.DEFAULTDIRS={"ProjectDir":self.ProjectDir,
                              "LogDir":self.LogDir,
                              "PreProcessDir":self.PreProcessDir,
                              "ImagesDir":self.ImagesDir,
                              "ML2DDir":self.ML2DDir,
                              "RotSpectraDir":self.RotSpectraDir,
                              "RCTDir":self.RCTDir,
                              "ML3DDir":self.ML3DDir,
                              "ProjMatchDir":self.ProjMatchDir,
                              "HighRes3DDir":self.HighRes3DDir,
                              }
            

            # Perform the actual setup:
            if (self.AutoLaunch!=""):
                # A. Setup from GUI (Autolaunch)
                # This will copy the (modified) protocol script to the corresponding directory
                # and will automatically launch the GUI for this protocol
                self.setup_protocol(self.library[self.AutoLaunch][1],
                                    self.library[self.AutoLaunch][2])
            else:
                # B. Setup from this script:
                # This will only copy the (modified) protocol script to the corresponding directory
                for var in self.library:
                    if (self.library[var][0]):
                        self.setup_protocol(self.library[var][1],
                                            self.library[var][2])

        def modify_script_header(self,src):

            import sys
            # Read script header and body
            fh=open(src,'r')
            header_lines=[]
            body_lines=[]
            isheader=False
            while (isheader!=True):
                line=fh.readline()
                if line=="":
                    print "Error, this file does not have a {end-of-header} label"
                    sys.exit()
                header_lines.append(line)
                if "{end-of-header}" in line:
                    isheader=True

            body_lines=fh.readlines()
            fh.close()

            # Loop over all project-related directories and preset default directories
            for dir in self.DEFAULTDIRS:
                newheader_lines=[]
                for line in header_lines:
                    if ((not line[0]=="#" and
                         not line[0]==" " and
                         not line[0]=="\t" and
                         not line[0]=="\n" and
                         not line[0]=="\"") and (dir in line)):
                        args=line.split("=")
                        lineb=str(args[0])+'=\"'+self.DEFAULTDIRS[dir]+'\"\n'
                    else:
                        lineb=line
                    newheader_lines.append(lineb)
                header_lines=newheader_lines
            return header_lines+body_lines

        def setup_protocol(self,directory,scripts):
            import os
            import shutil

            os.chdir(self.ProjectDir)
            if not os.path.exists(directory):
                os.makedirs(directory)
                
            for script in scripts:
                src=str(self.SYSTEMSCRIPTDIR)+"/"+str(script)
                dst=str(directory)+"/"+str(script)
                if os.path.exists(dst):
                    src=dst
                    print "* File "+dst+" already existed (now updated)"

                print src,dst,directory
                text=self.modify_script_header(src)
                fh=open(dst,'w')
                fh.writelines(text)
                fh.close()

            # Only auto-launch the first script
            script=scripts[0]
            if (self.AutoLaunch!=""):
                os.chdir(directory)
                command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+str(script)+' &'
                os.system(command)
                os.chdir(self.ProjectDir)


#
# Main
#
if __name__ == '__main__':

    import sys
    if (len(sys.argv) < 2):
        AutoLaunch=""
    else:
        AutoLaunch=sys.argv[1]

    setup=setup_protocols_class(SetupPreProcessMicrographs,
                                SetupParticlePick,
                                SetupPreProcessParticles,
                                SetupML2D,
                                SetupKerdensom,
                                SetupRotSpectra,
                                SetupRCT,
                                SetupML3D,
                                SetupProjMatch,
                                SetupHighRes3D,
                                ProjectDir,
                                LogDir,
                                ImagesDir,
                                PreProcessDir,
                                ML2DDir,
                                RotSpectraDir,
                                RCTDir,
                                ML3DDir,
                                ProjMatchDir,
                                HighRes3DDir,
                                AutoLaunch)

