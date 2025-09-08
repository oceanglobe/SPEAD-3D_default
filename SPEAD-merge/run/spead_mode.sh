#!/bin/sh
##----------------------------------------------------------

read -p "Number of spatial dimensions ? (1/3)" ndim
read -p "Discrete (D) / Continuous (C) approach ?" DC
read -p "First output year ?" nyear

#nday is 360*nyear if the output if the last year
#During test it is also possible to output only a few days
#nday=1
nday=$((360*nyear))
nstep=$((nday*8))

nb=7
if [ "$DC" = "D" ]
then
	nx=7
	ny=7
	nz=7
	ntrac=$((nx*ny*nz+1))
elif [ "$DC" = "C" ]
then
	ntrac=11
fi
	nptr=$((ntrac+nb))

##----------------------------------------------------------

#In SIZE.h, change domain size.
sed -i "/#define DIMENSIONS/c\\#define DIMENSIONS $ndim" ../code/SIZE.h 

#In DARWIN_OPTIONS.h, change trait diversity approach
if [ "$DC" = "D" ]
then
	sed -i "s/#define SPEAD_CONTINUOUS_TRAIT/#undef  SPEAD_CONTINUOUS_TRAIT/g" ../code/DARWIN_OPTIONS.h
	sed -i "s/#define SPEAD_CONTINUOUS_COVARIANCE/#undef  SPEAD_CONTINUOUS_COVARIANCE/g" ../code/DARWIN_OPTIONS.h
fi
if [ "$DC" = "C" ]
then
	sed -i "s/#undef  SPEAD_CONTINUOUS_TRAIT/#define SPEAD_CONTINUOUS_TRAIT/g" ../code/DARWIN_OPTIONS.h
	sed -i "s/#undef  SPEAD_CONTINUOUS_COVARIANCE/#define SPEAD_CONTINUOUS_COVARIANCE/g" ../code/DARWIN_OPTIONS.h
fi

#In data, data.off and data.darwin, change input files
if [ $ndim -eq 1 ]
then
        cp save/data_1D ../data/data
        cp save/data.off_1D ../data/data.off
	##--------------------------------------------------
	sed -i 's/ ironFile/#ironFile/g' ../data/data.darwin
	sed -i 's/ iceFile/#iceFile/g' ../data/data.darwin
	sed -i 's/ windFile/#windFile/g' ../data/data.darwin
	sed -i "0,/PARFile/ s/#PARFile/ PARFile/g" ../data/data.darwin
	sed -i "0,/PARFile/! {0,/PARFile/ s/ PARFile/#PARFile/g}" ../data/data.darwin
	##--------------------------------------------------
fi
if [ $ndim -eq 3 ]
then
        cp save/data_3D ../data/data
        cp save/data.off_3D ../data/data.off
	##--------------------------------------------------
	sed -i 's/#ironFile/ ironFile/g' ../data/data.darwin
	sed -i 's/#iceFile/ iceFile/g' ../data/data.darwin
	sed -i 's/#windFile/ windFile/g' ../data/data.darwin
	sed -i "0,/PARFile/ s/ PARFile/#PARFile/g" ../data/data.darwin
	sed -i "0,/PARFile/! {0,/PARFile/ s/#PARFile/ PARFile/g}" ../data/data.darwin
	##--------------------------------------------------
fi

#In data change number of time steps
echo $nstep
sed -i -E "s/( nTimeSteps=).*/\1$nstep,/g" ../data/data

#In data.diagnostics, change number of levels, output frequency and time phase
if [ $ndim -eq 1 ]
then
	# 1 day output frequency
	sed -i -E 's/(frequency\([0-9]+\)=).*/\1 86400.,/g' ../data/data.diagnostics
	# all depth levels 
	sed -i 's/ levels/#levels/g' ../data/data.diagnostics
fi
if [ $ndim -eq 3 ]
then
	# 1 month output frequency
	sed -i -E 's/(frequency\([0-9]+\)=).*/\1 2592000.,/g' ../data/data.diagnostics
	# only upper levels (number must be specified manually in data.diagnostics)
	sed -i 's/#levels/ levels/g' ../data/data.diagnostics
fi
sed -i -E "s/(timephase\([0-9]+\)=).*/\1 $((31104000*(nyear-1))).,/g" ../data/data.diagnostics

#In data.diagnostics, activate or deactivate TRAC outputs
#sed -i -E "s/( fields\().*(,1)/#fields\(\1,1/g" ../data/data.diagnostics
sed -i -E "s/ fields(\([0-9]+,1)/#fields\1/g" ../data/data.diagnostics
for i in `seq 1 $nptr`
do
	sed -i "s/#fields($i,1)/ fields($i,1)/g" ../data/data.diagnostics
done

##----------------------------------------------------------

#Build data.ptracers and data_ptracers_names.m from scratch
> data.ptracers

echo " &PTRACERS_PARM01" >> data.ptracers
echo " PTRACERS_numInUse= $nptr," >> data.ptracers
echo " PTRACERS_Iter0= 0," >> data.ptracers
echo "#" >> data.ptracers
echo " PTRACERS_advScheme(:)= $nptr*33," >> data.ptracers
echo " PTRACERS_diffKh(:)= $nptr*0.0," >> data.ptracers
echo " PTRACERS_diffKr(:)= $nptr*1D-05," >> data.ptracers
echo " PTRACERS_useGMRedi(:)= $nptr*.TRUE.," >> data.ptracers
echo " PTRACERS_useKPP(:)= $nptr*.TRUE.," >> data.ptracers
echo "#" >> data.ptracers
# The tracers should be renumbered if the compilation options are changed
#echo " ptracers_units(  )= 'mmol C/m^3'," >> data.ptracers
echo " ptracers_units( 1)= 'mmol N/m^3'," >> data.ptracers
echo " ptracers_units( 2)= 'mmol N/m^3'," >> data.ptracers
echo " ptracers_units( 3)= 'mmol N/m^3'," >> data.ptracers
echo " ptracers_units( 4)= 'mmol P/m^3'," >> data.ptracers
echo " ptracers_units( 5)= 'mmol Fe/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol Si/m^3'," >> data.ptracers
echo " ptracers_units( 6)= 'mmol C/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol N/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol P/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol Fe/m^3'," >> data.ptracers
echo " ptracers_units( 7)= 'mmol C/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol N/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol P/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol Fe/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol Si/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol C/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol eq/m^3'," >> data.ptracers
#echo " ptracers_units(  )= 'mmol O/m^3'," >> data.ptracers

if [ "$DC" = "D" ]
then
	echo " ptracers_units($((nb+1)):$nptr)= $ntrac*'mmol C/m^3'," >> data.ptracers
elif [ "$DC" = "C" ]
then
	echo " ptracers_units($((nb+1)))= 'mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+2)))= 'log(ESD) mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+3)))= 'degC mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+4)))= 'log(PAR) mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+5)))= 'log(ESD)^2 mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+6)))= '(degC)^2 mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+7)))= 'log(PAR)^2 mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+8)))= 'log(ESD) degC mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+9)))= 'log(ESD) log(PAR) mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+10)))= 'degC log(PAR) mmol C/m^3'," >> data.ptracers
	echo " ptracers_units($((nb+11)))= 'mmol C/m^3'," >> data.ptracers
fi
echo "#" >> data.ptracers

if [ $ndim -eq 1 ]
then
	echo " PTRACERS_initialFile( 4)= '../../SPEAD-1D/input/loc1_PO4_janprof_sargasso.bin'," >> data.ptracers
	echo "  ptracers_ref(:,1)= 80*1.6," >> data.ptracers
	echo "  ptracers_ref(:,2)= 80*0.1," >> data.ptracers
	echo "  ptracers_ref(:,3)= 80*0.1," >> data.ptracers
	echo "  ptracers_ref(:,6)= 80*0.0," >> data.ptracers
	echo "  ptracers_ref(:,7)= 80*0.0," >> data.ptracers
elif [ $ndim -eq 3 ]
then
	echo " "
        # echo " PTRACERS_initialFile( 1)='/home/legland/LEGLAND/SPEAD-3D/input/glodap_tco2_ann-3d.bin',"
        echo " PTRACERS_initialFile( 1)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_dimNO3_run82.bin'," >> data.ptracers
        echo " PTRACERS_initialFile( 3)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_dimNH4_run82.bin'," >> data.ptracers
        echo " PTRACERS_initialFile( 4)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_dimPO4_run82.bin'," >> data.ptracers
        echo " PTRACERS_initialFile( 5)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_dimFe2_run82.bin'," >> data.ptracers
        # echo " PTRACERS_initialFile( 7)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_dimSiO_run82.bin',"
        echo " PTRACERS_initialFile( 6)='/home/legland/LEGLAND/SPEAD-3D/input/doc_y2000.data'," >> data.ptracers
        # echo " PTRACERS_initialFile( 9)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_domNO3_run82.bin',"
        # echo " PTRACERS_initialFile(10)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_domPO4_run82.bin',"
        # echo " PTRACERS_initialFile(11)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_domFe2_run82.bin',"
        echo " PTRACERS_initialFile(7)='/home/legland/LEGLAND/SPEAD-3D/input/poc_y2000.data'," >> data.ptracers
        # echo " PTRACERS_initialFile(13)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_pomNO3_run82.bin',"
        # echo " PTRACERS_initialFile(14)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_pomPO4_run82.bin',"
        # echo " PTRACERS_initialFile(15)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_pomFe2_run82.bin',"
        # echo " PTRACERS_initialFile(16)='/home/legland/LEGLAND/SPEAD-3D/input/ICP/MITesm3D_pomSiO_run82.bin',"
        # echo " PTRACERS_initialFile(17)='/home/legland/LEGLAND/SPEAD-3D/input/pic_y2000.data',"
        # echo " PTRACERS_initialFile(18)='/home/legland/LEGLAND/SPEAD-3D/input/alk.bin',"
        # echo " PTRACERS_initialFile(19)='/home/legland/LEGLAND/SPEAD-3D/input/o2_molm3.bin',"
	echo "  ptracers_ref(:,3)= 23*0.0," >> data.ptracers
fi

echo "$resultat"
if [ "$DC" = "D" ]
then
	if [ $ndim -eq 1 ]
	then
		biomasse=0.75
		resultat=$(echo "scale=5; $biomasse / $ntrac" | bc)
		for i in `seq 1 $((ntrac-1))`
		do 
			echo "  ptracers_ref(:,$((nb+i)))= 80*$resultat," >> data.ptracers
		done
		echo "  ptracers_ref(:,$nptr)= 80*$biomasse," >> data.ptracers
	elif [ $ndim -eq 3 ]
	then
		for i in `seq 1 $((ntrac-1))`
		do 
			echo "  ptracers_ref(:,$((nb+i)))= 80*0.0001," >> data.ptracers
		done
		echo "  ptracers_ref(:,$nptr)= 80*0.01," >> data.ptracers
	fi
elif [ "$DC" = "C" ]
then
	if [ $ndim -eq 1 ]
	then
		echo "  ptracers_ref(:,$((nb+1)))= 80*0.75," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+2)))= 80*4.5," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+3)))= 80*18.0," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+4)))= 80*2.25," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+5)))= 80*27.75," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+6)))= 80*432.75," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+7)))= 80*6.825," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+8)))= 80*108.0," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+9)))= 80*13.5," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+10)))= 80*54.0," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+11)))= 80*0.75," >> data.ptracers
	elif [ $ndim -eq 3 ]
	then
		echo "  ptracers_ref(:,$((nb+1)))= 23*0.01," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+2)))= 23*0.03," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+3)))= 23*0.30," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+4)))= 23*0.03," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+5)))= 23*0.10," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+6)))= 23*9.10," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+7)))= 23*0.10," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+8)))= 23*0.90," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+9)))= 23*0.09," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+10)))= 23*0.90," >> data.ptracers
		echo "  ptracers_ref(:,$((nb+11)))= 23*0.01," >> data.ptracers
	fi
fi
echo "#" >> data.ptracers

echo " PTRACERS_names(1)= 'NO3'," >> data.ptracers
echo " PTRACERS_names(2)= 'NO2'," >> data.ptracers
echo " PTRACERS_names(3)= 'NH4'," >> data.ptracers
echo " PTRACERS_names(4)= 'PO4'," >> data.ptracers
echo " PTRACERS_names(5)= 'FeT'," >> data.ptracers
echo " PTRACERS_names(6)= 'DOC'," >> data.ptracers
echo " PTRACERS_names(7)= 'POC'," >> data.ptracers

if [ "$DC" = "D" ]
then
	for i in `seq 1 $ntrac`
	do 
		echo " PTRACERS_names($((i+7)))= 'c$i'," >> data.ptracers
	done
elif [ "$DC" = "C" ]
then
	echo " PTRACERS_names($((nb+1)))= 'c01'," >> data.ptracers
	echo " PTRACERS_names($((nb+2)))= 'xm01'," >> data.ptracers
	echo " PTRACERS_names($((nb+3)))= 'ym01'," >> data.ptracers
	echo " PTRACERS_names($((nb+4)))= 'zm01'," >> data.ptracers
	echo " PTRACERS_names($((nb+5)))= 'xxv01'," >> data.ptracers
	echo " PTRACERS_names($((nb+6)))= 'yyv01'," >> data.ptracers
	echo " PTRACERS_names($((nb+7)))= 'zzv01'," >> data.ptracers
	echo " PTRACERS_names($((nb+8)))= 'xyc01'," >> data.ptracers
	echo " PTRACERS_names($((nb+9)))= 'xzc01'," >> data.ptracers
	echo " PTRACERS_names($((nb+10)))= 'yzc01'," >> data.ptracers
	echo " PTRACERS_names($((nb+11)))= 'c02'," >> data.ptracers
fi
echo " /" >> data.ptracers

#---------------------------------

> data_ptracers_names.m

echo "function [PTRACERS_names] = data_ptracers_names()" >> data_ptracers_names.m
echo "" >> data_ptracers_names.m

echo " PTRACERS_names{1}= 'NO3';" >> data_ptracers_names.m
echo " PTRACERS_names{2}= 'NO2';" >> data_ptracers_names.m
echo " PTRACERS_names{3}= 'NH4';" >> data_ptracers_names.m
echo " PTRACERS_names{4}= 'PO4';" >> data_ptracers_names.m
echo " PTRACERS_names{5}= 'FeT';" >> data_ptracers_names.m
echo " PTRACERS_names{6}= 'DOC';" >> data_ptracers_names.m
echo " PTRACERS_names{7}= 'POC';" >> data_ptracers_names.m

if [ "$DC" = "D" ]
then
	for i in `seq 1 $ntrac`
	do 
		echo " PTRACERS_names{$((nb+i))}= 'c$i';" >> data_ptracers_names.m
	done
elif [ "$DC" = "C" ]
then
	echo " PTRACERS_names{$((nb+1))}= 'c01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+2))}= 'xm01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+3))}= 'ym01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+4))}= 'zm01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+5))}= 'xxv01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+6))}= 'yyv01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+7))}= 'zzv01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+8))}= 'xyc01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+9))}= 'xzc01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+10))}= 'yzc01';" >> data_ptracers_names.m
	echo " PTRACERS_names{$((nb+11))}= 'c02';" >> data_ptracers_names.m

fi

echo "" >> data_ptracers_names.m
echo " return" >> data_ptracers_names.m

mv data.ptracers ../data/data.ptracers
mv data_ptracers_names.m ../data/data_ptracers_names.m

##----------------------------------------------------------

#These last changes are only for the DISCRETE mode
if [ "$DC" = "D" ]
then
	#In DARWIN_SIZE.h and PTRACERS_SIZE.h, adapt the number of tracers (DISCRETE)
	sed -i "/#else \/\* NOT spead_continuous_trait \*\//{
	n; s/.*/      parameter\(nTrac=$ntrac\) /
	n; s/.*/      parameter\(nPlank=$ntrac\) /
	n; s/.*/      parameter\(nPhyp=$((ntrac-1))\) /
	n; s/.*/      parameter\(nPhoto=$((ntrac-1))\) /
        }" ../code/DARWIN_SIZE.h
	sed -i "/#else \/\* NOT spead_continuous_trait \*\//{
	n; s/.*/      PARAMETER\(PTRACERS_num = $nptr\) /
        }" ../code/PTRACERS_SIZE.h

	#In data.darwin change the number of trait values along each axis (DISCRETE)
	sed -i -E "s/( grp_nbiovol(:)=).*/\1         $nx, 1, ! Discrete/g" ../data/data.darwin
	sed -i -E "s/( grp_ntopt(:)=).*/\1           $ny, 1, ! Discrete/g" ../data/data.darwin
	sed -i -E "s/( grp_nparopt(:)=).*/\1         $nz, 1, ! Discrete/g" ../data/data.darwin
fi

##----------------------------------------------------------
##----------------------------------------------------------

#Compile SPEAD in 1D or 3D mode
read -p "Compile SPEAD ? (Y/N)" comp
if [ "$comp" = "Y" ]
then
        cd ../make
	make Clean
        ##--------------------------------------------------
        export ROOTDIR=${HOME}/LEGLAND/MITgcm
        ##--------------------------------------------------
	if [ $ndim -eq 1 ]
	then
                ${ROOTDIR}/tools/genmake2 \
                -of ${ROOTDIR}/tools/build_options/linux_amd64_gfortran -mods=../code
	fi
	if [ $ndim -eq 3 ] 
	then 
		export MPI_INC_DIR=/lib/x86_64-linux-gnu/openmpi/include
                ##------------------------------------------
                ${ROOTDIR}/tools/genmake2 -mpi \
                -of ${ROOTDIR}/tools/build_options/linux_amd64_gfortran -mods=../code
	fi
        ##--------------------------------------------------
        make depend
        make -j8
        ##--------------------------------------------------
        rm ../run/mitgcmuv
        mv mitgcmuv ../run
fi
##----------------------------------------------------------

