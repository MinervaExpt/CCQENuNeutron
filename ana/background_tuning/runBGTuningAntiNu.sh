#signal=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_Signal.root
#blob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_BlobSideBand.root
#michel=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MichelSideBand.root
#micblob=/minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_new_MicBlobSideBand.root

#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV/
#folder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_flux_CV_Elastic_FSI/
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_2_optimized_bins_
#basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Proton_test_3_reweight_optimized_bins2_20200215_
basefolder=/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Neutron_test_6_optimizing_bins_
basefolder=../rootfiles/minervameAntiNu/


#for tag in CV CV_Absorption_FSI CV_All_FSI CV_Elastic_FSI
#for tag in pion_rpa
for tag in CV
do
  playlist=All
  folder=$basefolder$tag/
  signal=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_${playlist}.root
  blob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand_${playlist}.root
  michel=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-MichelSideBand_${playlist}.root
  micblob=$folder/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-MicBlobSideBand_${playlist}.root

  #nohup ./SideBandFit $signal $blob $michel $micblob bkgd_weights_minervame1D_q2qe_flux_$tag.root > log_$tag &
  savefolder=../rootfiles/minervameAntiNu/${tag}/
  mkdir -p $savefolder
  #ln -s $folder/*root $savefolder
  ./SideBandFit2 $signal $blob $michel $micblob ${savefolder}/bkgd_weights_${playlist}_ptmu_flux_${tag}_2D_optimized.root
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
