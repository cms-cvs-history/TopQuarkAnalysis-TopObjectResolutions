#!/bin/sh

# This script allows to submit the ResolutionCreator jobs in parallel, and to merge and fit the 
# resolution (Et,Eta) bins for all objects afterwards (Tau still to be added!!!)
# author: J.Heyninck
# date: 30the of August 2007

# this should be your working dir
cd /scratch/delaer/top/step3/CMSSW_1_3_6/src/TopQuarkAnalysis/TopObjectResolutions/test

# protocol to open a ROOT-file + storage element
storageAndProtocol="file://"
#set storageAndProtocol 	= ""

# the area where the output of all the ResolutionCreator jobs arrived
inputDir="/scratch/delaer/top/step3/CMSSW_1_3_6/src/TopQuarkAnalysis/TopObjectResolutions/test/output"
#set inputDir   		= "."

# all (Et,Eta) bin content will be merged automatically
# with this flag also the fitting and the filling of the tree can be turned on 
dofit="true"

# suppose you have output files in different locations, then you want to do the merging in two steps: 
# first all files in one location, then the merged files from the different locations 
# this flag will be appended to the output name  
#set outputFlag 		= "tt0j"
outputFlag=""

# the target (electron, muon, ...
target="$1"

rm -f command.csh
touch command.csh
echo -n "hadd -f Resolutions_bins_${target}.root " >> command.csh
for f in `ls $inputDir | grep ${target}`; do
  echo -n "${storageAndProtocol}$inputDir/$f " >> command.csh
done
chmod u+x command.csh
./command.csh
rm command.csh
if [ "$dofit" = "true" ]; then
  root -l -b -q "ExtractFitCurvesFromMergedBins.C(\"Resolutions_bins_${target}.root\",\"Resolutions_${target}.root\")"
  rm -f Resolutions_bins_${target}.root
else 
  mv Resolutions_bins_${target}.root Resolutions_bins_${target}_${outputFlag}.root
fi
