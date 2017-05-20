#!/bin/bash
#Usage : ./shell <FileToCopy>
#     ex)./shell Aug.cc

if [[ -z $1 ]]; then echo "Write down the file to copy"; exit 1; fi;
FileToCp=$( basename $1 .cc )

cp ${FileToCp}.cc ${FileToCp}_D3lv.cc
find . -name "${FileToCp}_D3lv.cc" -type f -exec sed -i s/${FileToCp}/${FileToCp}_D3lv/g {} +
cp ${FileToCp}.cc ${FileToCp}_D2l2j.cc
find . -name "${FileToCp}_D2l2j.cc" -type f -exec sed -i s/${FileToCp}/${FileToCp}_D2l2j/g {} +

sed -i '/Basic\ Objects\ Distribution/ i \ \ \ //GenDecayCut///////////////////////////////////' ${FileToCp}_D3lv.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ Decay3lv=GenDecayInfo(truthColl,"3lv");' ${FileToCp}_D3lv.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ if(!Decay3lv) return;' ${FileToCp}_D3lv.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ ////////////////////////////////////////////////' ${FileToCp}_D3lv.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ ' ${FileToCp}_D3lv.cc


sed -i '/Basic\ Objects\ Distribution/ i \ \ \ //GenDecayCut///////////////////////////////////' ${FileToCp}_D2l2j.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ Decay3lv=GenDecayInfo(truthColl,"3lv");' ${FileToCp}_D2l2j.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ if(Decay3lv) return;' ${FileToCp}_D2l2j.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ ////////////////////////////////////////////////' ${FileToCp}_D2l2j.cc
sed -i '/Basic\ Objects\ Distribution/ i \ \ \ ' ${FileToCp}_D2l2j.cc
