#!/bin/bash
# Author : Jihwan Bhyun, Jan. 20. 2017
# Context : kill all the processes submitted with LQanalyzer or kill some specific process jobs.
# Usage : ./killer <LogPath>

KillAll="true"  ## true, True, false, False accepted.
#KillAll="false"
declare -a ToKillProcess=('TT_powheg')

if [[ ! -d $1 ]]; then echo "Log path wrong, exiting"; exit 1; fi
if [[ -e Assassin.sh ]]; then rm Assassin.sh; echo "Removed previous assassin code."; fi

if [[ $KillAll == "true" || $KillAll == "True" ]];
then
  echo "#!/bin/bash" > Assassin.sh
  grep -r "kill" $1 | sed "s/^.*source/source/g" >> Assassin.sh
elif [[ $KillAll == "false" || $KillAll == "False" ]];
then
  echo "#!/bin/bash" > Assassin.sh
  for Proc in "${ToKillProcess[@]}"; do
    grep -r "kill.*${Proc}" $1 | sed "s/^.*source/source/g" >> Assassin.sh
  done
else
  echo "Clarify what you want. Exiting"; exit 1;
fi

chmod 700 Assassin.sh
./Assassin.sh
