#!/bin/bash

#nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -S SingleMuon -p G -m "Processing Missing sample(periodG)" -c v8-0-4 -q longq
nohup bash ${LQANALYZER_BIN_PATH}/submitSKTree.sh -b True -M True -a SKTreeMaker -S SingleElectron -p E -m "Processing Missing sample(periodE)" -c v8-0-4 -q fastq
