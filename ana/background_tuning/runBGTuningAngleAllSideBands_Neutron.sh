model=""

mult=1
neutrino="minervameAntiNu"
tag=""

source ../include/scripts_source.sh

#for tune in  CV2 CV2Elastic CVElasticNuwroSF CVElasticNuwroSFn CVElasticNuwroSFn2
#for tune in  CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 #pionrpaElasticNuwroSF
for tune in  CV2ElasticNuwroSFn 
#for tune in  CV2ElasticNuwroSF CV2 CV2Elastic 
do
  #playlist=All

  folder=../rootfiles/${neutrino}/cuts${cuts}/${tune}${tag}/
  folder2=../rootfiles/${neutrino}/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/${tune}${tag}/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-Signal-Angles_CombinedPlaylists${model}.root
  blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-BlobSideBand-Angles_CombinedPlaylists${model}.root
  michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MichelSideBand-Angles_CombinedPlaylists${model}.root
  micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MicBlobSideBand-Angles_CombinedPlaylists${model}.root
  
  #signalweight=15
  #lambda="50"
  #signalweight="30"
  #lambda="6"
  #lambda="4900000"

  signalweight="1"
  #lambda="3085"
  #lambda="450"
  lambda="550"
  lambda="230"
  lambda="320"
  lambda="500"
  lambda="400"
  lambda="550"
  lambda="390"
  lambda="200"
  #lambda="20"
  #lambda="200"
  #lambda="15"

  #output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_qe_2p2h_snp_scp_mp.root
  output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_qe_res_2p2h_snp_scp_3.root

  savefolder=$folder
  #mkdir -p $savefolder
  #ln -sf $folder/*root ${savefolder}
  #cmd="./SideBandFitAngleAllSamples_Neutron $signal $blob $michel $micblob  $output ${signalweight} ${lambda}"
  cmd="./SideBandFitAngleAllSamples_Neutron $signal $blob $michel $micblob  $output ${signalweight} ${lambda}"
  echo $cmd
  $cmd
  #nohup $cmd > nohup_${tune}$cuts.log &
done


#signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root
#blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_minervame1D.root
#michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_minervame1D.root
#micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_minervame1D.root
##nohup ./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root > log02 &
#./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
