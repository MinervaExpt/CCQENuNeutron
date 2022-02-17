#signal=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_Signal.root
#blob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_BlobSideBand.root
#michel=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MichelSideBand.root
#micblob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MicBlobSideBand.root

#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV/
#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV_Elastic_FSI/
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_optimized_bins_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Processing_R02_20200323/
basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Processing_20200408/
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Neutron_test_6_optimizing_bins_


#for tag in CV CV_Absorption_FSI CV_All_FSI CV_Elastic_FSI
#for tag in pion_rpa
model=""
#model="_nuwroSF"

mult=1
neutrino="minervameNu"

#for tag in CV2 CV2Elastic #CV2ElasticNuwroLFG CV2ElasticNuwroSF # #CV2 CV2ElasticCai CV2Elastic #CV CV2 CV2Elastic CV2ElasticCai pion_rpa default
#for tag in CV2 CV2Elastic CV2ElasticNuwroSF CV2ElasticNuwroLFG 
for tag in CV2ElasticNuwroSFn 
do
  #playlist=All

  folder=../rootfiles/${neutrino}/${tag}/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-Signal-Angles_CombinedPlaylists${model}.root
  blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-BlobSideBand-Angles_CombinedPlaylists${model}.root
  michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MichelSideBand-Angles_CombinedPlaylists${model}.root
  micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-MicBlobSideBand-Angles_CombinedPlaylists${model}.root
  #signalweight=15
  #lambda="50"
  signalweight="20"
  lambda="550"

  #output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_noRES.root
  output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_noRES_SignalBlob.root
  #output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}${model}.root

  savefolder=$folder
  #mkdir -p $savefolder
  #ln -sf $folder/*root ${savefolder}
  cmd="./SideBandFitAngleAllSamples_Neutron $signal $blob $michel $micblob  $output ${signalweight} ${lambda}"
  echo $cmd
  nohup $cmd  &
  #$cmd
  #./SideBandFit $signal $blob $michel $micblob ./bkgd_weights_minervame1D_ptmu_flux_$tag.root
done


#signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root
#blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_minervame1D.root
#michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_minervame1D.root
#micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_minervame1D.root
##nohup ./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root > log02 &
#./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
