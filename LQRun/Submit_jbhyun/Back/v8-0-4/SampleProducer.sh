#!/bin/bash

#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -S SingleMuon -p G -m "Processing Missing sample(periodG)" -c v8-0-4 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleMuon -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S MuonEG -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleEG -p ALL -m "First Production of TriLep Skimmed Datasetv804" -c v8-0-4 -q longq

#List="TT"
List="SignalMajor_All"
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -list ${List} -c ${CATVERSION} -m "First production of majot H+>WA sample in 806" -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerDiLep -list ${List} -c ${CATVERSION} -m "First production of major H+>WA sample in 806" -q longq
nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -list ${List} -c ${CATVERSION} -m "First production of major H+>WA sample in 806" -q longq
