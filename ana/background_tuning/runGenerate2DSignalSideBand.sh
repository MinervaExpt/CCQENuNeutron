#signal=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_Signal.root
#blob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_BlobSideBand.root
#michel=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MichelSideBand.root
ProcessFolder=Neutron
neutrino=minervameAntiNu
multiplicity=1

source ../include/scripts_source.sh

export MALLOC_CHECK_=0

#for tag in pionrpaElasticNuwroSF CVElasticNuwroSF CVElasticNuwroSFn CVElasticNuwroSFn2 CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2Elastic CV2
for tag in CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2
do
  for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
  do 
    folder=../rootfiles/${neutrino}/cuts${cuts}_2D/${tag}/
    signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_CombinedPlaylists.root
    output=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_CombinedPlaylists${model}.root
    nplaylist=12

    mkdir -p $folder
    #ln -sf $folder/*root ${savefolder}
    cmd="./GenerateSideBandHists2D $signal $output $nplaylist 1"
    echo $cmd
    $cmd 
  done
  #echo $cmd
  #./SideBandFit $signal $blob $michel $micblob ./bkgd_weights_minervame1D_ptmu_flux_$tag.root
done

