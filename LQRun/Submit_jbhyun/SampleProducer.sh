#!/bin/bash

#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -S SingleMuon -p G -m "Processing Missing sample(periodG)" -c v8-0-4 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleMuon -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S MuonEG -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq
nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleEG -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq

#List="TT"
#List="SignalMajor_1e2mu"
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -list ${List} -c ${CATVERSION} -m "Reproduction of 1e2mu channel for triggers"
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerDiLep -list ${List} -c ${CATVERSION} -m "First Time production of dilepskimmed powheg ttbar sample" -q longq
