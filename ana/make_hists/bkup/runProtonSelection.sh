#-----------------------------------------------------------------------------------------------
#MACROS HELP:
#
#        -./ProtonSelectionHists  Name_and_path_to_Output_file Playlist Selection_Signal  Make_Flux_Constraint_Histo  Multiplicity  Number_MC_files  Number_DATA_files 
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
cmd="./ProtonSelectionHistsMap2 /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode.root minervame1A Signal 0 2 -1 -1 2>&1 >log/proton_minervame1A.log &"
echo $cmd
#nohup ./ProtonSelectionHistsMap2 /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode_wrongSign.root minervame1A Signal 0 2 -1 -1 2>&1 >log/proton_minervame1A.log &
for mode in MichelSideBand BlobSideBand MicBlobSideBand
#for mode in Signal
do
  cmd="nohup ./ProtonSelectionHistsMap2 /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode_$mode.root minervame1A $mode 0 2 -1 -1 2>&1 >log/proton_minervame1A_$mode.log &"
  echo $cmd
  eval $cmd
done
#nohup ./ProtonSelectionHistsMap2 /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode_MichelSideBand.root minervame1A MichelSideband 0 2 -1 -1 2>&1 >log/proton_minervame1A.log &
#nohup ./ProtonSelectionHistsMap2 /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode_MichelSideBand.root minervame1A MichelSideband 0 2 -1 -1 2>&1 >log/proton_minervame1A.log &
#eval $cmd
