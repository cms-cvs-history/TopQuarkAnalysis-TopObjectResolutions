#!/bin/tcsh

# This script allows to submit the ResolutionCreator jobs in parallel, and to merge and fit the 
# resolution (Et,Eta) bins for all objects afterwards (Tau still to be added!!!)
# author: J.Heyninck
# date: 30the of August 2007

# this should be your working dir
cd /beo5/heyninck/CMSSW/src/TopQuarkAnalysis/TopObjectResolutions/test/

# protocol to open a ROOT-file + storage element
set storageAndProtocol 		= "dcap://maite.iihe.ac.be"
#set storageAndProtocol 	= ""

# the area where the output of all the ResolutionCreator jobs arrived
set inputDir   			= "/pnfs/iihe/becms/heyninck/ResOutput_tt0j"
#set inputDir   		= "."

# all (Et,Eta) bin content will be merged automatically
# with this flag also the fitting and the filling of the tree can be turned on 
set dofit      			= "true"

# suppose you have output files in different locations, then you want to do the merging in two steps: 
# first all files in one location, then the merged files from the different locations 
# this flag will be appended to the output name  
#set outputFlag 		= "tt0j"
set outputFlag 			= ""



# electron
touch command.csh
echo -n "hadd Resolutions_bins_electron.root " >> command.csh
foreach f(`ls $inputDir|grep electron`)
  echo -n "${storageAndProtocol}$inputDir/$f " >> command.csh
end
chmod u+x command.csh
./command.csh
rm command.csh
if ($dofit == "true") then
  root -l -q 'ExtractFitCurvesFromMergedBins.C("Resolutions_bins_electron.root","Resolutions_electron.root")'
  rm -f Resolutions_bins_electron.root
else 
  mv Resolutions_bins_electron.root Resolutions_bins_electron_${outputFlag}.root
endif
