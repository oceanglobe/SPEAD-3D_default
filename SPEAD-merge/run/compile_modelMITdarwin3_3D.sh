#!/bin/sh
##------------------------------------------
cd ../make
make Clean
##------------------------------------------
export ROOTDIR=${HOME}/LEGLAND/MITgcm
export MPI_INC_DIR=/lib/x86_64-linux-gnu/openmpi/include
##------------------------------------------
${ROOTDIR}/tools/genmake2 -mpi \
    -of ${ROOTDIR}/tools/build_options/linux_amd64_gfortran -mods=../code
##------------------------------------------
make depend
# I think the following step is necessary to avoid mixing with the next task
#read -t 10 -p "Waiting 10 seconds before compiling"
make -j8
#read -t 10 -p "Waiting 10 seconds before moving mitgcmuv"
##------------------------------------------
mv mitgcmuv ../run
##------------------------------------------
