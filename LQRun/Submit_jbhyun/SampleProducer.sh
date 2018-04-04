#!/bin/bash

#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -S SingleMuon -p G -m "Processing Missing sample(periodG)" -c v8-0-6 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleMuon -p ALL -m "First Production of TriLep Skimmed Datasetv806" -c v8-0-6 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S MuonEG -p ALL -m "First Production of TriLep Skimmed Datasetv806" -c v8-0-6 -q longq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -S DoubleEG -p ALL -m "First Production of TriLep Skimmed Datasetv806" -c v8-0-6 -q longq

#List="SignalMajor_All"
List="ListToProd"
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -list ${List} -c ${CATVERSION} -m "First production of new sample" -q fastq -SIG
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -list ${List} -c ${CATVERSION} -m "First production of new sample" -q fastq 
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerDiLep -list ${List} -c ${CATVERSION} -m "First production of new sample" -q fastq
#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerDiLep -list ${List} -c ${CATVERSION} -m "First production of new sample" -q fastq -SIG
nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMakerTriLep -list ${List} -c ${CATVERSION} -m "First production of new sample" -q fastq -SIG
