######################################################
#use -m64 sompiler option for all code
setenv COMPILER_BIT 64
setenv g2psoft /work/halla/g2p/disk1/centos62
setenv g4soft ${g2psoft}/geant
######################################################
###!set up cernlib

setenv CERN /usr/lib64/cernlib
setenv CERN_LEVEL 2006
setenv CERN_ROOT $CERN/$CERN_LEVEL
setenv PATH ${CERN_ROOT}/bin:${PATH}
if !($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH ${CERN_ROOT}/lib:/usr/lib64:/usr/local/lib64
else
    setenv LD_LIBRARY_PATH ${CERN_ROOT}/lib:${LD_LIBRARY_PATH}
endif
######################################################
###!set up root

setenv ROOTSYS ${g2psoft}/root-5.34.09-x86_64
setenv ROOTLIB $ROOTSYS/lib
setenv PATH ${ROOTSYS}/bin:${PATH}
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
######################################################
###!set up XERCES
setenv XERCESCROOT ${g2psoft}/xercesc
setenv PATH ${XERCESCROOT}/bin:${PATH}
setenv LD_LIBRARY_PATH ${XERCESCROOT}/lib:${LD_LIBRARY_PATH}
######################################################
###!set up CLHEP

#setenv  CLHEP_BASE_DIR     ${g4soft}/clhep/2.1.3.1
#setenv  CLHEP_INCLUDE_DIR  ${CLHEP_BASE_DIR}/include
#setenv  CLHEP_LIB_DIR      ${CLHEP_BASE_DIR}/lib
#setenv  CLHEP_LIB          CLHEP
#setenv  LD_LIBRARY_PATH ${CLHEP_LIB_DIR}:${LD_LIBRARY_PATH}
######################################################
###!set up Geant4

setenv  G4VERSION  9.6.3
setenv  G4SYSTEM   Linux-g++
setenv  G4INSTALL  ${g4soft}/share/Geant4-${G4VERSION}/geant4make
setenv  G4INCLUDE  ${g4soft}/include/Geant4
setenv  G4LIB      ${g4soft}/lib64/Geant4-${G4VERSION}

setenv G4DATADIR  ${g4soft}/share/Geant4-${G4VERSION}
setenv G4NEUTRONHPDATA   ${G4DATADIR}/data/G4NDL4.2
setenv G4LEDATA          ${G4DATADIR}/data/G4EMLOW6.32    
setenv G4LEVELGAMMADATA  ${G4DATADIR}/data/PhotonEvaporation2.3
setenv G4RADIOACTIVEDATA ${G4DATADIR}/data/RadioactiveDecay3.6
setenv G4NEUTRONXSDATA   ${G4DATADIR}/data/G4NEUTRONXS1.2       
setenv G4PIIDATA         ${G4DATADIR}/data/G4PII1.3                   
setenv G4REALSURFACEDATA ${G4DATADIR}/data/RealSurface1.0     
setenv G4SAIDXSDATA      ${G4DATADIR}/data/G4SAIDDATA1.1   

setenv G4UI_USE_TCSH  1

#setenv G4UI_BUILD_XAW_SESSION  1
#setenv G4UI_USE_XAW  1
#setenv XAWFLAGS  ""
#setenv XAWLIBS " -lXaw "

#setenv G4UI_BUILD_XM_SESSION  1
#setenv G4UI_USE_XM  1
#setenv XMFLAGS  ""
#setenv XMLIBS  " -lXm -lXpm "

setenv G4UI_BUILD_QT_SESSION  1
setenv G4UI_USE_QT  1

#setenv G4VIS_BUILD_DAWN_DRIVER  1
#setenv G4VIS_USE_DAWN  1

setenv G4VIS_BUILD_OPENGLX_DRIVER  1
setenv G4VIS_USE_OPENGLX  1

#setenv G4VIS_BUILD_OPENGLXM_DRIVER  1
#setenv G4VIS_USE_OPENGLXM  1

#setenv G4VIS_BUILD_RAYTRACERX_DRIVER  1
#setenv G4VIS_USE_RAYTRACERX  1

#setenv G4VIS_BUILD_VRML_DRIVER  1
#setenv G4VIS_USE_VRML  1

setenv G4VIS_BUILD_OPENGLQT_DRIVER  1
setenv G4VIS_USE_OPENGLQT  1

 
setenv G4LIB_BUILD_GDML  1
#setenv G4LIB_BUILD_ZLIB  1
#setenv G4LIB_BUILD_STATIC  1

#use the default qt of centos6.2
setenv QTDIR /usr
setenv QTLIBPATH ${QTDIR}/lib64
setenv QTMOC ${QTDIR}/bin/moc-qt4
#setenv QTDIR ${g4soft}/qt/4.8.5
#setenv QTLIBPATH ${QTDIR}/lib
#setenv QTMOC ${QTDIR}/bin/moc
setenv QTFLAGS "-I${QTDIR}/include -I${QTDIR}/include/QtCore -I${QTDIR}/include/QtGui  -I${QTDIR}/include/QtOpenGL"
setenv QTLIBS "-L${QTLIBPATH} -lQtCore_debug -lQtGui -lQtXml"
setenv GLQTLIBS "-L${QTLIBPATH} -lQtCore_debug -lQtGui -lQtOpenGL -lQtXml"

setenv  G4WORKDIR  /u/scratch/$USER/g4work_${G4VERSION}_${OS_NAME}
if !(-e  $G4WORKDIR) mkdir -p $G4WORKDIR
if !(-d  $G4WORKDIR) then 
    setenv G4WORKDIR  $HOME/.g4work_${G4VERSION}_${OS_NAME}
endif
#set G4BIN in home dir such that every pc can access, if set it in scratch
#some batch farm node will not see it since they do not share the same scratch disk 
#setenv G4BIN ${G4WORKDIR}/bin
setenv G4BIN $HOME/.g4work_${G4VERSION}_${OS_NAME}/bin

#SET G4TMP so one can recompile the lib with minimum rebuild
setenv G4TMP ${G4WORKDIR}/tmp

setenv PATH ${G4BIN}/${G4SYSTEM}:${g4soft}/bin:${PATH}
setenv LD_LIBRARY_PATH ${QTLIBPATH}:${G4LIB}/${G4SYSTEM}:${LD_LIBRARY_PATH}

