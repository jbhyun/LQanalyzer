#!/bin/bash

ps | egrep python | sed s/^[\ ]*//g | cut -d ' ' -f1 > ToKillList.txt
while read line
do
  kill -9 ${line}
done<ToKillList.txt

if [[ -e ToKillList.txt ]]; then rm ToKillList.txt; fi
