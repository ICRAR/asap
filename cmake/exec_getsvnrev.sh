#! /bin/bash

DIR=.
if [ $# -ge 1 ]; then
    DIR=$1
fi

pushd ${DIR}
${DIR}/getsvnrev.sh
pushd
