#signal=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_Signal.root
#blob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_BlobSideBand.root
#michel=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MichelSideBand.root
#micblob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MicBlobSideBand.root

#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV/
#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV_Elastic_FSI/
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_optimized_bins_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Processing_R02_20200323/
basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Processing_20200408/
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Neutron_test_6_optimizing_bins_


#for tag in CV CV_Absorption_FSI CV_All_FSI CV_Elastic_FSI
#for tag in pion_rpa
ProcessFolder=Proton
for tag in CV #CV2 CV2Elastic pion_rpa default
do
  #playlist=All
  folder=../rootfiles/minervameNu/${tag}/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal-Angles_CombinedPlaylists.root
  bgweight=${folder}/bkgd_weights_minervameNu_q2qe_flux_2D_optimized_Type.root
  output=${folder}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type.root

  savefolder=$folder
  #mkdir -p $savefolder
  #ln -sf $folder/*root ${savefolder}
  cmd="./SideBandFitAngle2 $signal $bgweight  $output"
  echo $cmd
  $cmd
  #./SideBandFit $signal $blob $michel $micblob ./bkgd_weights_minervame1D_ptmu_flux_$tag.root
done


#signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root
#blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_minervame1D.root
#michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_minervame1D.root
#micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_minervame1D.root
##nohup ./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root > log02 &
#./SideBandFit_Tejin $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_elastic.root
#nohup ./SideBandFit $signal $blob $michel $micblob test00.root > log00 &
