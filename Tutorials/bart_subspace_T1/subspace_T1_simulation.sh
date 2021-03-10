#!/bin/bash
# 
# Copyright 2021. Uecker Lab, University Medical Center Goettingen.
#
# Author: Xiaoqing Wang, 2020-2021
# xiaoqing.wang@med.uni-goettingen.de
#

set -e

export PATH=$TOOLBOX_PATH:$PATH

if [ ! -e $TOOLBOX_PATH/bart ] ; then
	echo "\$TOOLBOX_PATH is not set correctly!" >&2
	exit 1
fi

TR=0.0041
DIM=384
SPOKES=1
REP=1020
NC=8
NBR=$((DIM/2))

# create trajectory
bart traj -x $DIM -y $SPOKES -t $REP -c -r -G _traj
bart transpose 5 10 {_,}traj
bart scale 0.5 traj traj1

# create geometry basis functions
bart phantom -s$NC -T -k -b -t traj1 _basis_geom

# create simulation basis functions
bart signal -F -I -n$REP -r$TR  -1 3:3:1 -2 1:1:1 _basis_simu_water
bart signal -F -I -n$REP -r$TR  -1 0.2:2.2:10 -2 0.045:0.045:1 _basis_simu_tubes

bart scale 1. _basis_simu_tubes _basis_simu_sdim_tubes
bart join 6 _basis_simu{_water,_sdim_tubes,}

# create simulated dataset
bart fmac -s $(bart bitmask 6) _basis_geom _basis_simu phantom_ksp
bart phantom -x$NBR -T mask

# add noise to the simulated dataset 
bart noise -n500 phantom_ksp phantom_ksp_1

# create a reference T1 (noiseless) 
bart index 6 10 tmp_T1s
bart scale 0.2 tmp_T1s tmp_T1s_1
bart ones 7 1 1 1 1 1 1 10 tmp_ones_T1s
bart saxpy 0.2 tmp_ones_T1s tmp_T1s_1 tmp_T1s_2
bart ones 7 1 1 1 1 1 1 1 tmp_ones_T1s_1
bart scale 3.0 tmp_ones_T1s_1 tmp_ones_T1s_2
bart join 6 tmp_ones_T1s_2 tmp_T1s_2 ref_T1s
bart phantom -T -b -x$NBR phan_T1
bart fmac -s $(bart bitmask 6) phan_T1 ref_T1s tmp
bart transpose 0 1 tmp phan_ref_T1s
bart invert phan_ref_T1s phan_ref_R1s

# create a zero map 
dim0=$(bart show -d0 phan_ref_R1s)
dim1=$(bart show -d1 phan_ref_R1s)
bart zeros 2 $dim0 $dim1 zero_map

# Coil sensitivity estimation using NLINV

bart extract 5 $((REP-400)) $REP traj traj_state
bart extract 5 $((REP-400)) $REP phantom_ksp_1 phantom_ksp_state
bart transpose 2 5 traj_state traj_state1
bart transpose 2 5 phantom_ksp_state phantom_ksp_state1

ITER=8
DEBUG=4
bart nlinv -d$DEBUG -i$ITER -a1000 -t traj_state1 phantom_ksp_state1 reco sens_nlinv
bart resize -c 0 $DIM 1 $DIM sens_nlinv sens_nlinv1
bart scale 100. sens_nlinv1 sens

# data binning
nspokes_per_frame=10
bart reshape $(bart bitmask 4 5) $nspokes_per_frame $((REP/nspokes_per_frame)) traj traj1
bart transpose 4 2 traj1 traj

bart reshape $(bart bitmask 4 5) $nspokes_per_frame $((REP/nspokes_per_frame)) phantom_ksp_1 phantom_ksp_2
bart transpose 4 2 phantom_ksp_2 phantom_ksp_1

# time vector
TR1=4100
bart index 5 $((REP/nspokes_per_frame)) tmp1.ra
# use local index from newer bart with older bart
#./index 5 $num tmp1.coo
bart scale $(($nspokes_per_frame * $TR1)) tmp1.ra tmp2.ra
bart ones 6 1 1 1 1 1 $((REP/nspokes_per_frame)) tmp1.ra
bart saxpy $((($nspokes_per_frame / 2) * $TR1)) tmp1.ra tmp2.ra tmp3.ra
bart scale 0.000001 tmp3.ra TI

#-----------------------------------------------
#------------- Linear subspace reco ------------
#-----------------------------------------------

# Generate a dictionary for IR LL 
nR1s=1000
nMss=100
TR2=0.041
bart signal -F -I -1 5e-3:5:$nR1s -3 1e-2:1:$nMss -r$TR2 -n$((REP/nspokes_per_frame)) dicc

bart reshape $(bart bitmask 6 7) $((nR1s * nMss)) 1 dicc dicc1
bart squeeze dicc1 dicc2
bart svd -e dicc2 U S V

export ITER=100
export REG=0.0005

nCoe=4
bart extract 1 0 $nCoe U basis
bart transpose 1 6 basis basis1
bart transpose 0 5 basis1 basis_${nCoe}

bart pics -SeH -d5 -RW:$(bart bitmask 0 1):$(bart bitmask 6):$REG -i$ITER -t traj -B basis_${nCoe} phantom_ksp_1 sens subspace_reco
bart resize -c 0 $NBR 1 $NBR subspace_reco tmp
bart fmac mask tmp tmp_masked

# project "nCoe" coefficient maps to images and 
# perform pixel-wise fitting to obtain T1 map
bart fmac -s $(bart bitmask 6) basis_${nCoe} tmp_masked imgs
python3 mapping_piecewise.py imgs T1 TI maps

bart extract 2 0 3 maps tmp1
bart transpose 2 6 tmp1 tmp2

bart looklocker -t0.0 -D0.0 tmp2 tmp3
bart scale 0.5 tmp3 T1map
bart fmac mask T1map tmp4 
bart transpose 0 1 tmp4 subspace_T1map

bart saxpy -- -1.0 subspace_T1map phan_ref_T1s subspace_diff_T1map
bart fmac subspace_diff_T1map phan_ref_R1s subspace_T1_rela_diff
python3 save_maps.py subspace_T1map viridis 0 2.0 subspace_T1map.png
python3 save_maps.py subspace_T1_rela_diff viridis 0 0.4 subspace_T1_rela_diff.png