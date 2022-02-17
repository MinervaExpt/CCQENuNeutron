#!/bin/sh

source ../include/scripts_source_nu.sh

playlist=minervameNu
saveplaylist=NP
func="./plot_pn_eventrate_regions_neutron_proton_func"
multiplicity=2

#for tune in  CV2 CV2Elastic CV2ElasticNuwroSF CV2AllElasticNuwroSF pionrpaElastic pionrpaAllElastic #CV2 CV2ElasticNuwroSF #CV2Elastic CV2ElasticNuwroSF CV2ElasticNuwroLFG #CV2 CV2Elastic #CV CV2 CV2Elastic CV2ElasticCai pion_rpa default;

for cuts in _np_NoEB_openAngle_newCat5_blob_2_angle_1 #_np_NoEB_openAngle_newCat5 _np_NoEB_openAngle_newCat5_pass10Sel _np_NoEB_openAngle_newCat5_pass10SelNew
do
for tune in  CV2Elastic #CV2ElasticNuwroSF CV2AllElasticNuwroSF pionrpaElastic pionrpaAllElastic #CV2 CV2ElasticNuwroSF #CV2Elastic CV2ElasticNuwroSF CV2ElasticNuwroLFG #CV2 CV2Elastic #CV CV2 CV2Elastic CV2ElasticCai pion_rpa default;
do

  w=1
  for lam in 0
  do
    for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
    do
      inputdir="../rootfiles/${playlist}/cuts${cuts}/${tune}/"
      savefolder=plots/${saveplaylist}/cuts${cuts}/${tune}/${sample}_${w}_${lam}/
      echo save in $savefolder

      input="${inputdir}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_CombinedPlaylists.root"
      weight="${inputdir}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_snp_scp_mp.root"
      #weight="${inputdirWeight}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_1_0_snp_scp_mp.root"
      #weight="${inputdir}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp.root"
      cmd="${func} ${input}  ${weight} "
      echo $cmd
      $cmd
      mkdir -p $savefolder
      mv plots/*pdf $savefolder
    done
  done
done
done
