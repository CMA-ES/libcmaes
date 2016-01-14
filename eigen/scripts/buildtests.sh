#!/bin/bash

if [[ $# != 1 || $1 == *help ]]
then
  echo "usage: ./check regexp"
  echo "  Builds tests matching the regexp."
  echo "  The EIGEN_MAKE_ARGS environment variable allows to pass args to 'make'."
  echo "    For example, to launch 5 concurrent builds, use EIGEN_MAKE_ARGS='-j5'"
  exit 0
fi

TESTSLIST=""
targets_to_make=`echo "$TESTSLIST" | egrep "$1" | xargs echo`

if [ -n "${EIGEN_MAKE_ARGS:+x}" ]
then
  /usr/local/bin/gmake $targets_to_make ${EIGEN_MAKE_ARGS}
else
  /usr/local/bin/gmake $targets_to_make 
fi
exit $?

