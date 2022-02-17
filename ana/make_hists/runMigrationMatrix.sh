#-----------------------------------------------------------------------------------------------
#MACROS HELP:
#
#        -./NeutronSelectionHists  Name_and_path_to_Output_file Playlist Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files 
#
#        -Name_and_path_to_Output_file    =       Name of and path to the Output ROOT file that will be created 
#        -Playlist        =       Name of the playlist you want to run over e.g. minerva1 
#        -Selection_Signal        =       Can be: Signal, MichelSideBand, SideBand (this is:non-vtx vs Q2 sideband) 
#        -Make_Flux_Constraint_Histo      =       If TRUE Enter 1; If FALSE Enter 0 
#        -Multiplicity    =       Enter 0, 1 or 2. Here 0: 1 or 2 tracks; 1: 1-track events only; 2: 2-track events only 
#        -Number_MC_files         =       Number of MonteCarlo files. To use all files, set this to: -1 
#        -Number_DATA_files       =       Number of Data files. To use all files, set this to: -1
#-----------------------------------------------------------------------------------------------



outputfolder=/minerva/data/users/tejinc/CCQENu/antinu/hists/
filename=neutronselhists.root
cmd="./MuonSelection $outputfolder$filename minervame5A Signal 0 1 -1 -1"
echo $cmd

#export DOPIONTUNE=1
#export DO2P2HTUNE=1
#export DOCCQERPA=1
#export DOMINERVALOWQ2=1
#export DOELASTICFSIQE=1
#export DOQ2KFREWEIGHTSF=1

#export DONEUTRONCVRW=1



playlist=6H
fluxHist=1
applyFlux=1
multiplicity=1
sample="hydrogen"
npMode=0
nFiles=1
cmd="./MigrationMatrixHists ./test_Mig.root minervame${playlist} $fluxHist $applyFlux $multiplicity $sample ${nFiles}  $npMode Q2"
echo $cmd
$cmd
