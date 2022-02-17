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
cmd="./RecoilEnergyStudies /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_minervame5A.root minervame5A Signal 0 1 -1 -1 2>&1 >log/proton_minervame5A.log &"
echo $cmd
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6A.root minervame6A Signal 0 1 -1 -1 0 0  2>&1 >log/minervame6A.log &
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6D.root minervame6D Signal 0 1 -1 -1 0 0  2>&1 >log/minervame6D.log &
./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame5A.root minervame5A Signal 0 1 -1 -1 0 0  
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6B.root minervame6B Signal 0 1 -1 -1 0 0  2>&1 >log/minervame6B.log &
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6C.root minervame6C Signal 0 1 -1 -1 0 0 2>&1 >log/minervame6C.log &
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6E.root minervame6E Signal 0 1 -1 -1 0 0 2>&1 >log/minervame6E.log &
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6F.root minervame6F Signal 0 1 -1 -1 0 0 2>&1 >log/minervame6F.log &
#nohup ./NeutronPurityStudy /minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6G.root minervame6G Signal 0 1 -1 -1 0 0 2>&1 >log/minervame6G.log &
#eval $cmd
