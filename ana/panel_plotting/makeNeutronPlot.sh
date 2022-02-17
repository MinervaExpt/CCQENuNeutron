#!/bin/sh

#processTag=_Recoil
processTag=""
Use2D=""

playlist=minervameNu
multiplicity=2

source ../include/scripts_source_nu.sh

playlist=minervameAntiNu
multiplicity=1

source ../include/scripts_source.sh

for tag in CV2ElasticNuwroSFn #CV2ElasticNuwroSFz CV2ElasticNuwroSFn CV2ElasticNuwroSFnz #CV2ElasticNuwroSF CV2ElasticNuwroSFz #CV2ElasticNuwroSFz CV2 CV2Elastic #pionrpaElasticNuwroSF  CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2 CV2Elastic
#for tag in CV2ElasticNuwroSFn
do
  oldtag=$tag
  tag=${tag}${Use2D}
  neutrino=0
  nPlaylist=7


   w=1
  #for lam in 550 #2500 #2300 #10 4900000
  for lam in 200  #390 #550 #400 #500 #400 #320 #230 #550 #2500 #2300 #10 4900000
  do
    inputdir="../rootfiles/${playlist}/cuts${cuts}/${tag}/"
    signalweight="../rootfiles/${playlist}/cuts${cuts}/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp_3.root"
    #signalweight2="../rootfiles/minervameAntiNu/cuts${cuts}/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp_3.root"

    #for sample in Signal #BlobSideBand #MichelSideBand MicBlobSideBand
    for sample in Signal BlobSideBand #MichelSideBand MicBlobSideBand
    do
      #func="./plot_nu_eventrate_regions_neutron_proton_reweights"
      inputorig="${inputdir}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_CombinedPlaylists.root"
      inputregion="${inputdir}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_CombinedPlaylists.root"

      func="./plot_nu_eventrate1" #drawfit, drawchunk3d
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      ###$cmd

      func="./plot_nu_eventrate2" #regions/systematics
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      $cmd

      func="./plot_nu_eventrate3" #diagnostic
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      ###$cmd

      func="./plot_nu_eventrate4" #recoil
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      ###$cmd

      func="./plot_nu_eventrate5" #1D distro
      cmd="${func} ${inputorig} ${inputregion} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
      $cmd


      savefolder=plots/${playlist}${processTag}/cuts${cuts}/${tag}/${sample}_${w}_${lam}/
      savefolder=plots/${playlist}${processTag}/cuts${cuts}/${tag}/redraw_${sample}_${w}_${lam}/
      mkdir -p $savefolder
      mv plots/*pdf $savefolder
    done
  done
done
