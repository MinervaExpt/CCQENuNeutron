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
cmd="./MuonSelection $outputfolder$filename minervame1D Signal 0 1 -1 -1"
echo $cmd
playlist=6D

export DOPIONTUNE=1
export DO2P2HTUNE=1
export DOCCQERPA=1
export DOMINERVALOWQ2=1
export DOELASTICFSIQE=1
export DOQ2KFREWEIGHTSF=1

for playlist in 6G #6G #5A 6A 6B #6C 6D 6E 6G
do
  nohup ./RecoilEnergyStudy recoil_${playlist}.root minervame${playlist} Signal 0 1 -1 -1 0  &
done
  #nohup ./MuonSelectionHistsMap /minerva/data/users/tejinc/CCQENu/antinu/hists/test_neutronSel_minervame5A_$mode.root minervame5A $mode 0 1 -1 -1 2>&1 >log/neutron_minervame5A_$mode.log &
  #nohup ./MuonSelectionHistsMap /minerva/data/users/tejinc/CCQENu/antinu/hists/test_neutronSel_minervame6A_$mode.root minervame6A $mode 0 1 -1 -1 2>&1 >log/neutron_minervame6A_$mode.log &
  #nohup ./MuonSelectionHistsMap /minerva/data/users/tejinc/CCQENu/antinu/hists/test_neutronSel_minervame6B_$mode.root minervame6B $mode 0 1 -1 -1 2>&1 >log/neutron_minervame6B_$mode.log &
#eval $cmd
