#MACROS HELP:
#      ./CrossSectionHists  Output_file  Input_file  Input_file_BKG  Input_file_MIGRATION  Input_file_EFFICIENCY FluxConstraint Num_iterations  Run_Type DoENU 
#      -output_file = Cross Section output file
#      -input_file = Muon Event Selection(Signal) for Neutron samples 
#      -input_file_angles = Muon Event Selection(Signal) for Neutron Region samples 
#      -input_file_bkg = Background weights filename
#      -input_file_migration = Migration Histogram
#      -input_file_efficiency = Efficiency Histograms
#      -UseFluxConstraint = UseFluxConstraint 
#      -Number of Iterations = Number of Iters  
#      -Run Type, unused = can be any integer  
#      -Do cross section with flux normalization =  0 or 1  
#      ********************************************************************************************** 
#       Please see : MuonSelectionHists.cxx, BackgroundWeights.cxx, MigrationMatrixHists.cxx, EffPurityHists.cxx and Ana/Flux/python/compute_flux.py for getting the necessary input files std::endl; 
#


neutrinoMode=minervameAntiNu
multiplicity=1
w=1
lam=550
lam=200

flux=0
runType=-1
DoENU=0

source ../include/scripts_source.sh

#for tune in  CV2ElasticNuwroSF   CV2ElasticNuwroSFn2   CV2ElasticNuwroSFn   
for lam in 200 #390
do
  for tune in    CV2ElasticNuwroSFn   #CV2ElasticNuwroSF   CV2Elastic CV2
  do
    for iterations in  4 #5 #6 7 8 9 10
    #for iterations in 15 20  #6 7 8 9 10
    do
      for pl in CombinedPlaylists
      do
      #pl=CombinedPlaylists
      #pl=minervame6D
      folder=../rootfiles/${neutrinoMode}/cuts${cuts}/${tune}/
      signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-Signal-Angles_${pl}.root
      bcg=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_snp_scp_3.root
      mig=${folder}/Migration_$pl.root
      eff=${folder}/EffPurity_$pl.root

      output=${folder}/CrossSectionHist_It${iterations}_w_lambda_${w}_${lam}_${pl}
      cmd="./CrossSectionHists_Tejin $output $signal $bcg $mig $eff  $flux $iterations $runType $DoENU"
      $cmd


      eff=${folder}/EffPurity_FullKine_$pl.root
      output=${folder}/CrossSectionHist_FullKine_It${iterations}_w_lambda_${w}_${lam}_${pl}
      cmd="./CrossSectionHists_Tejin $output $signal $bcg $mig $eff  $flux $iterations $runType $DoENU"
      $cmd

      output=${folder}/CrossSectionHist2D_It${iterations}_w_lambda_${w}_${lam}_${pl}
      cmd="./CrossSectionHists_w2D_Tejin $output $signal $bcg $mig $eff  $flux $iterations $runType $DoENU"
      echo $cmd
      #$cmd
      done
    done
  done
done
