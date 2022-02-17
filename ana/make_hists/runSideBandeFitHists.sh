neutrinoMode=minervameAntiNu
multiplicity=1
flux=0

basefolder=../rootfiles/${neutrinoMode}/
for tag in CV2ElasticNuwroSF CV2ElasticNuwroLFG CV2Elastic CV2
do

  folder=$basefolder/$tag/
  signal=MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root
  bckweight=signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_20_1000_noRES.root
  out=Hydrogen_Data_Bcksub.root
  ./SideBandFitHists_Tejin ${folder}/$out ${folder}/${signal} ${folder}/${bckweight} ${flux} ${flux}
done
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
