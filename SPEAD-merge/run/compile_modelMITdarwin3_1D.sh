#!/bin/sh
##------------------------------------------
cd ../make
make Clean
##------------------------------------------
export ROOTDIR=${HOME}/LEGLAND/MITgcm
##------------------------------------------
# linux_amd64_gfortran modified to include more warnings (Le Gland, 17/03/2021)
${ROOTDIR}/tools/genmake2 \
    -mods=../code \
#   -of ../run/linux_amd64_gfortran
    -of ${ROOTDIR}/tools/build_options/linux_amd64_gfortran
##------------------------------------------
make depend
# I think the following step is necessary to avoid mixing with the next task
read -t 10 -p "Waiting 10 seconds before compiling"
make -j8
#read -t 10 -p "Waiting 10 seconds before moving mitgcmuv"
##------------------------------------------
rm ../run/mitgcmuv
mv mitgcmuv ../run
##------------------------------------------
