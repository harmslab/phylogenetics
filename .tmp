#!/bin/sh
#install.sh

# Installs *.py package by copying all .py files to install_dir. 

install_dir=${1}
if [ ${install_dir} ]; then
    current_dir=`pwd`
    cd ${install_dir}
    for script in `ls ${current_dir}/*.py`; do
        cp ${script} .
    done
    for script in `ls ${current_dir}/*.sh`; do
        if [ ${script} != "install.sh" ]; then
            cp ${script} .
        fi
    done
else
    echo "No install directory specified.  Exiting."
    echo "    install.sh install_dir"
    exit
fi

