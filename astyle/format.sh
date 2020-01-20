#!/bin/bash
echo pwd:$PWD
astyle_executable=$1
astylerc=$2
do_format=$3
files_to_format="""
EKF/AlphaFilter.hpp
EKF/RingBuffer.h
"""
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

if [[ $do_format -eq 1 ]]
then
    echo formatting
    $astyle_executable ${files_to_format} --options=${astylerc} --preserve-date\
    echo -e ${GREEN}Formatting finished${NC}
else
    echo checking format...
    $astyle_executable --dry-run ${files_to_format} --options=${astylerc} --preserve-date | grep Formatted &>/dev/null
    if [[ $? -eq 0 ]]
    then
        echo -e ${RED}Error: need to format${NC}
        echo -e ${YELLOW}From cmake build directory run: 'make format'${NC}
        exit 1
    fi
    echo -e ${GREEN}no formatting needed${NC}
    exit 0
fi

# vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 :
