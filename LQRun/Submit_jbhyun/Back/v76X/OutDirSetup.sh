#!/bin/bash
# makes until catlep dir. ; Catout & CatVer & out|log & catlep

CATOUT='/data2/Users/jbhyun/cmssw/CATanalyzer/CATOut/'

if [[ -z ${CATVERSION} ]]; then echo "source setup.sh needed. exit"; exit 1; fi
if [[ ! -d $CATOUT ]]; then echo "mkdir $CATOUT"; mkdir $CATOut; fi
if [[ ! -d ${CATOUT}/${CATVERSION} ]]; then echo "mkdir ${CATOUT}/${CATVERSION}"; mkdir ${CATOUT}/${CATVERSION}; fi
CatLoc="${CATOUT}/${CATVERSION}"
echo "mkdir ${CatLoc}/output";              mkdir ${CatLoc}/output
echo "mkdir ${CatLoc}/logfiles";            mkdir ${CatLoc}/logfiles
echo "mkdir ${CatLoc}/output/Muon";         mkdir ${CatLoc}/output/Muon
echo "mkdir ${CatLoc}/output/Electron";     mkdir ${CatLoc}/output/Electron
echo "mkdir ${CatLoc}/output/ElectronMuon"; mkdir ${CatLoc}/output/ElectronMuon
