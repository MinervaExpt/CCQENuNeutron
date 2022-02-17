
mult=2
neutrino="minervameNu"
source ../include/scripts_source_nu.sh

#for tune in CV2 CV2Elastic CV2ElasticNuwroSF CV2AllElasticNuwroSF pionrpaElastic pionrpaAllElastic
for tune in CV2Elastic 
do
  folder=../rootfiles/${neutrino}/cuts${cuts}/${tune}${tag}/
    signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-Signal_CombinedPlaylists.root
      blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-BlobSideBand_CombinedPlaylists.root
    michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MichelSideBand_CombinedPlaylists.root
   micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MicBlobSideBand_CombinedPlaylists.root
  mkdir -p $folder

  signalweight="1"
  lambda="0"
  output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_snp_scp_mp.root
  cmd="./SideBandFitAngleAllSamples_PN $signal $blob $michel $micblob  $output ${signalweight} ${lambda}"
  echo $cmd
  $cmd
done


#signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root
#blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_minervame1D.root
#michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_minervame1D.root
#micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_minervame1D.root
##nohup ./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root > log02 &
#./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
