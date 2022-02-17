mult=1
neutrino="minervameAntiNu"
tag=""

source ../include/scripts_source.sh

for tune in  CV2ElasticNuwroSF #CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 pionrpaElasticNuwroSF
do
  folder=../rootfiles/${neutrino}/cuts${cuts}/${tune}${tag}/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-Signal_CombinedPlaylists${model}.root
  signal_angle=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-Signal-Angles_CombinedPlaylists${model}.root
  

  signalweight="1"
  lambda="490"

  wtag=qe_res_2p2h_snp_scp
  weight=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${signalweight}_${lambda}_${wtag}.root
  output=${folder}/MuonEventSelection_MakeFlux-1_Multiplicity-${mult}_Sample-EvtRate_CombinedPlaylists_${signalweight}_${lambda}_${wtag}.root

  savefolder=$folder
  #mkdir -p $savefolder
  #ln -sf $folder/*root ${savefolder}
  #cmd="./SideBandFitAngleAllSamples_Neutron $signal $blob $michel $micblob  $output ${signalweight} ${lambda}"
  cmd="./EventRateFitter $signal $signal_angle $weight $output"
  echo $cmd
  $cmd
  #nohup $cmd > nohup_${tune}$cuts.log &
done


