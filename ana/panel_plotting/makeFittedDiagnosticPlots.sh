#!/bin/sh

tag=CV2 
SampleType=Neutron 
#processTag=_Recoil
processTag=""
Use2D=""
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0"
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_newBlobCut-1"

source ../include/scripts_source.sh

#for tag in CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 #CVElasticNuwroSF pionrpaElasticNuwroSF  CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2 CV2Elastic
for tag in CV2ElasticNuwroSFn #CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 #CVElasticNuwroSF pionrpaElasticNuwroSF  CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2 CV2Elastic
do
  oldtag=$tag
  tag=${tag}${Use2D}
  playlist=minervameAntiNu
  neutrino=0
  multiplicity=1
  nPlaylist=7


   w=1
  #for lam in 550 #2500 #2300 #10 4900000
  for lam in 230 #2500 #2300 #10 4900000
  do
    signal="../studies/HadronReweightStudies/root_files/dbs_rw-0_CombinedPlaylists.root"
    weight="../rootfiles/minervameAntiNu/cuts${cuts}/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp.root"

    #func="./plot_nu_eventrate_regions_neutron_proton_reweights"

    func="./plot_nu_diagnostic_plots" #drawfit, drawchunk3d
    cmd="${func} ${signal} ${weight} "
    $cmd

    savefolder=plots/${playlist}${processTag}/studies_rw0/
    mkdir -p $savefolder
    mv plots/*pdf $savefolder
  done
done
