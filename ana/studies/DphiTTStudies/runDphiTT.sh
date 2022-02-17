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



export DOPIONTUNE=1
export DO2P2HTUNE=1
export DOCCQERPA=1
export DOMINERVALOWQ2=1
export DOELASTICFSIQE=1
export DOQ2KFREWEIGHTSF=1
#export DONEUTRONCVRW=1

#export DONEUTRONCVRW=1
#export DONEUTRONCVRW2=1
#./MuonSelectionHists2DNeutron test_muon_ela_SF_2D.root minervame${playlist} Signal 0 1 -1 -1 0 

#for pl in 6A 6B 6C 6D 6E 6F 6G
output="/minerva/data/users/tejinc/CCQENu/test/hadronRW/"
#playlists="5A 6A 6B 6C 6D 6E 6F 6G"
#playlists="6A 6B 6C 6D 6E 6F"
#playlists="5A 6A 6B 6C 6D 6E 6F 6G"
#playlists="5A"
#playlists="6B 6C 6D 6E 6F 6G"
#playlists="6G"
#playlists="1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P"
playlists="1A 1B 1C 1D 1E 1F 1G 1L 1M 1N 1O 1P"
modMin=0
modMax=0

rwType=1

#modMin=0
#modMax=40
#
#
#modMin=40
#modMax=70
#
#
#modMin=70
#modMax=100
for pl in $playlists
do
 
  mult=2
  flux=0
  fname="pndphitt_minervame${pl}.root"
  cmd="./DphiTTStudy ${output}/${fname} minervame${pl} Signal $flux $mult -1 -1 10 $modMin $modMax 0 $rwType"
  if [ "$pl" = "1E" ] || [ "$pl" = "1M" ] 
  then
    echo $cmd
    $cmd
  else
    echo "$pl to nohup"
    nohup $cmd > $fname.log 2>&1 &
  fi
  #$cmd
  ln -s ${output}/${fname} root_files/
  #./MuonSelectionHists2DNeutron test_muon_ela_SF_neutronW.root minervame${pl} Signal 0 1 -1 0 0

done


#export DOELASTICFSIQE=1
#export DOQ2KFREWEIGHTSF=0
#./MuonSelectionHists test_muon_ela.root minervame${playlist} Signal 0 1 5 10 0 
#
#export DOELASTICFSIQE=0
#export DOQ2KFREWEIGHTSF=0
#./MuonSelectionHists test_muon.root minervame${playlist} Signal 0 1 5 10 0 

