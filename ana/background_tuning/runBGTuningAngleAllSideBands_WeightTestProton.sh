
#for tag in CV CV_Absorption_FSI CV_All_FSI CV_Elastic_FSI
#for tag in pion_rpa



source ../include/scripts_source_nu.sh

for tune in CV2ElasticNuwroSF #CiV CV2 CV2Elastic pion_rpa default
do
  #playlist=All
  folder=../rootfiles/minervameNu/cuts${cuts}/${tune}/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal-Angles_CombinedPlaylists.root
  blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand-Angles_CombinedPlaylists.root
  michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand-Angles_CombinedPlaylists.root
  micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand-Angles_CombinedPlaylists.root

  savefolder=$folder
  #mkdir -p $savefolder
  #ln -sf $folder/*root ${savefolder}
  #for w in 1 3 5 10 15 20 50 100 200 500 1000 
  #for w in 0.5 1 2 3 4 5 10 14 15 16 17 18 19 20 21 25 30 40 100
  for w in 1 
  #for w in 1
  do
    #for lam in 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1  2  3 4 5 6 7 8 9 10 11 12 13 14 15  20  30  40 50 100 150 170 200 230 250 300 350 400 450 500 600 700 800 900 1000 1500 1800 1900 2000 2100 2150 2200 2250 2300 2400 2500 3000 4000 5000 6000 7000 8000 9000 10000 15000 20000
    for lam in 0  1  2  3 4 5 6 7 8 9 10 11 12 13 14 15  20  30  40 50 100 150 170 180 190 195 200 205 210 220 230 250 275 300 325 350 375 390 395 400 405 410 415 420 425 450 455 460 465 470 475 480 485 490 495 500 525 550 575 600 700 750 790 800 810 850 900 1000 1100 1200 1300 1400 1500 1650 1700 1800 1850 1890 1900 2000
    do
      mkdir -p ${folder}/test/
      #output=${folder}/test/test_signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_2p2h_scp_snp_noNorm_blob99.root
      output=${folder}/test/test_signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_${w}_${lam}_qe_res_2p2h_scp_snp_noNorm_blob99.root
      #rm $output
      if test -f "$output"
      then
        echo "File exists"
        rm $output
        continue
      fi
      cmd="./SideBandFitAngleAllSamples_test_proton $signal $blob $michel $micblob  $output $w $lam"
      echo $cmd
      $cmd 
    done
  done
done

