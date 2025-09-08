#!/bin/sh
##------------------------------------------
cp ../data/* .
./mitgcmuv > out.dat 
#read -t 10 -p "Waiting 10 seconds before moving files"
##------------------------------------------
rm ../output/meta-data/darwin_latest/*
#mv out.dat OUTPUT/META-DATA/darwin_latest/
mv *.data ../output/meta-data/darwin_latest/
mv *.meta ../output/meta-data/darwin_latest/
mv *.txt ../output/meta-data/darwin_latest/
##------------------------------------------
rm *data* 
##------------------------------------------

