neutrino=minervameNu
multiplicity=2
#cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-1"
#cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_newBlobCut-1"
#cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_newBlobCut-1"

source ../include/scripts_source_nu.sh

#for tag in  CV2ElasticNuwroLFG  CV2 CV2Elastic CV2ElasticCai pion_rpa default
#for tag in  CV2ElasticNuwroSFn #CV2ElasticNuwroLFG   CV2ElasticNuwroSF CV2 CV2Elastic CV2ElasticNuwroSFz 
#for tag in  CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 # CV2ElasticNuwroSFz CV2ElasticNuwroSFn2 #CV2ElasticNuwroLFG   CV2ElasticNuwroSF CV2 CV2Elastic CV2ElasticNuwroSFz 
#for tag in pionrpaElasticNuwroSF CVElasticNuwroSF CVElasticNuwroSFn CVElasticNuwroSFn2 CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2Elastic CV2
#for tag in CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2Elastic CV2
#for tag in CV2ElasticNuwroSFn CV2ElasticNuwroSF CV2ElasticNuwroSFnz CV2ElasticNuwroSFz
for tag in CV2ElasticNuwroSF CV2Elastic CV2
do
  #for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
  #for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
  for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
  do
    for pl in CombinedPlaylists
    do 
      folder=../rootfiles/${neutrino}/cuts${cuts}/${tag}/
      signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_$pl.root
      output=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_$pl.root
      nplaylist=12

      mkdir -p $folder
      #ln -sf $folder/*root ${savefolder}
      cmd="./GenerateSideBandHists $signal $output $nplaylist 1"
      echo $cmd
      $cmd
      #nohup $cmd  > nohup_${tag}_${sample}.log &2>1 &
    done
  done
    #echo $cmd
  #./SideBandFit $signal $blob $michel $micblob ./bkgd_weights_minervame1D_ptmu_flux_$tag.root
done


#signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root
#blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_minervame1D.root
#michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_minervame1D.root
#micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_minervame1D.root
##nohup ./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root > log02 &
#./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
