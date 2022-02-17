#!/bin/sh

tag=CV2 
SampleType=Neutron 
#processTag=_Recoil
processTag=""
Use2D=""
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0"
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_newBlobCut-1"

source ../include/scripts_source.sh

for tag in CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 #CVElasticNuwroSF pionrpaElasticNuwroSF  CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2 CV2Elastic
do
  oldtag=$tag
  tag=${tag}${Use2D}
  playlist=minervameAntiNu
  neutrino=0
  multiplicity=1
  nPlaylist=7


   w=1
  for lam in 550 #2500 #2300 #10 4900000
  do
    inputdir="../rootfiles/minervameAntiNu/cuts${cuts}_2D/${tag}/"
    #signalweight="../rootfiles/minervameAntiNu/cuts${cuts}/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_2p2h_snp_scp_mp.root"
    signalweight="../rootfiles/minervameAntiNu/cuts${cuts}/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp.root"

    for sample in Signal #BlobSideBand #MichelSideBand MicBlobSideBand
    do
      #func="./plot_nu_eventrate_regions_neutron_proton_reweights"
      inputorig="${inputdir}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_CombinedPlaylists.root"
      inputregion="${inputdir}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_CombinedPlaylists.root"

      func="./plot_nu_eventrate_regions_neutron_proton_2D"
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      echo $cmd
      $cmd

      savefolder=plots/${playlist}${processTag}/cuts${cuts2D}/${tag}/${sample}_${w}_${lam}/
      mkdir -p $savefolder
      mv plots/*pdf $savefolder
    done
  done
done
