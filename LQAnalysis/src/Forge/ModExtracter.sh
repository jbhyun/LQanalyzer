#!/bin/bash
# 23. Aug. 2016 Author: Jihwan Bhyun
# context: copy certain range of text from source file to destination file
#          For LQ merging stage, just from //My Modification ... to Endline seems enough
# ./shell <SourceFile> <DestFile>

SourceLoc=$1
DestLoc=$2

if [[ -z $CATTAG ]]; then echo "source setup needed, exit"; exit 1; fi;
if [[ ! -e ${SourceLoc} ]]; then echo "Source, Not Exists, exit"; exit 1; fi;
if [[ ! -e ${DestLoc} ]]; then echo "Dest, Not Exists, exit"; exit 1; fi
if [[ ${SourceLoc} = ${DestLoc} ]]; then echo "Source, Dest Same File, exit"; exit 1; fi
SourceName=$( readlink -f ${SourceLoc} | rev | cut -d '/' -f1 | rev )
DestName=$( readlink -f ${DestLoc} | rev | cut -d '/' -f1 | rev )
if [[ ${SourceName} != ${DestName} ]]; then echo "Source, Dest Diff File. But Not Same Name, exit"; exit 1; fi
if [[ ${SourceLoc} -nt ${DestLoc} ]];
then
  echo "Something weird, Source is Newer. Keep going? y,n?";
  read answer;
  if [[ $answer = 'n' ]]; then exit 1;
  elif [[ $answer = 'y' ]]; then echo "go";
  else echo "Answer Wrong. exit"; exit 1;
  fi
fi

LineN_Start=$( nl -ba ${SourceLoc} | egrep "Jihwan Bhyun Modification" | sed "s/^\ \ //g" | cut -d$'\t' -f1 )
LineN_End=$( wc -l ${SourceLoc} | cut -d ' ' -f1 )

#echo "${LineN_Start} ${LineN_End}"

sed -n ${LineN_Start},${LineN_End}p ${SourceLoc} >> ${DestLoc}
