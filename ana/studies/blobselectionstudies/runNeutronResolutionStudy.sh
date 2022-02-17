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
cmd="./NeutronResolutionStudies /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_minervame5A.root minervame5A SignalTest 0 1 -1 -1 0 2>&1 >log/proton_minervame5A.log &"
echo $cmd

nohup ./NeutronResolutionStudy2 /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_minervame6A.root minervame6A SignalTest 0 1 -1 -1 0 2>&1 >log/minervame6A.log &
nohup ./NeutronResolutionStudy2 /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_minervame6B.root minervame6B SignalTest 0 1 -1 -1 0 2>&1 >log/minervame6B.log &
./NeutronResolutionStudy2 /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_minervame5A.root minervame5A SignalTest 0 1 -1 -1 0 
#eval $cmd
