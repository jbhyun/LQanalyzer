export CATVERSION=v8-0-4
### If there is a small bug/new code then new subtag is made
export tag_numerator='.14'
if [[ $1 == 'branch' ]];
    then
    export CATTAG=
else
    export CATTAG=v8-0-4.14
fi

