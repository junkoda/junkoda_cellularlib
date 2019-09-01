#!/bin/bash

#
# sh push
# sh push -f   # Do push even if this branch is up to date
# sh push -a   # commit -a

# Increment patch version and git push

if [ "$1" = "-f" -o "$2" = "-f" ]; then
  # force update version
  :
else
  # check if local is ahead of origin
  if git status -uno | grep 'Your branch is ahead'; then
    :
  else
    echo "No need to push"
    exit
  fi
fi

if [ "$1" = "-a" -o "$2" = "-a" ]; then
  opt="-a"
fi


# Update py/ver
thisdir="$( cd "$(dirname "$0")" ; pwd -P )"
ver=`cat $thisdir/py/ver`
ver=`expr $ver + 1`
echo $ver > $thisdir/py/ver

# Commit py/ver
git add $thisdir/py/ver
git commit $opt -e -m "Push patch version $ver" || exit

git push || exit

echo "Pushed patch version $ver"
