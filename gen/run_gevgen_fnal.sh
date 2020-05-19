#! /usr/bin/env bash

####################

HORN=$1
NPER=$2
FIRST=$3
TEST=$4
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

if [ "${NPER}" = "" ]; then
echo "Number of events per job not specified, using 1000"
NPER=1000
fi

if [ "${FIRST}" = "" ]; then
echo "First run number not specified, using 0"
FIRST=0
fi

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "In TEST mode, assuming interactive running"
#CP="cp"
PROCESS=0
fi

MODE="neutrino"
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
fi

echo "Running gevgen for ${NPER} events in ${HORN} mode"

RNDSEED=$((${PROCESS}+${FIRST}))
NEVENTS="-n ${NPER}"      # No. of events, -e XE16 for POT

GEOMETRY="lar_mpt"
TOPVOL="volGArTPC"

USERDIR="/pnfs/dune/persistent/users/marshalc/CAF"
OUTDIR="/pnfs/dune/persistent/users/sbjones/CAF/genie"

OUTFLAG="GAr"
RDIR=0$((${RNDSEED} / 1000))

####################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup dk2nu        v01_05_01b -q e15:prof
setup ifdhc

####################

## Copy stuff to the local node
${CP} ${USERDIR}/gevgen.tar.gz gevgen.tar.gz
tar -xzf gevgen.tar.gz

# tarball has contents in a folder to avoid tarbombing while testing
mv gevgen/* .

# Get flux files to local node
chmod +x copy_dune_ndtf_flux
./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60

####################

## Add the location of the GNuMIFlux.xml to the GENIE xml path

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

## Run GENIE and copy output file to dCache persistent

gevgen_fnal \
    -f local_flux_files/gsimple*.root,DUNEND \
    -g ${GEOMETRY}.gdml \
    -t ${TOPVOL} \
    -m ${GEOMETRY}.${TOPVOL}.maxpl.xml \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC


echo "${CP} ${MODE}.${RNDSEED}.ghep.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root"
${CP} ${MODE}.${RNDSEED}.ghep.root ${OUTDIR}/${OUTFLAG}/${HORN}/${RDIR}/${OUTFLAG}.${MODE}.${RNDSEED}.ghep.root
rm -f ${MODE}.${RNDSEED}.ghep.root
