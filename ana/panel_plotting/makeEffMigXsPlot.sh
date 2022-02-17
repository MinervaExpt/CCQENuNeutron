#!/bin/sh

tag=CV2 
SampleType=Neutron 
#processTag=_Recoil
processTag=""
Use2D=""
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0"
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_newBlobCut-1"

source ../include/scripts_source.sh


w=1
lam=550
sample="Signal"
neutrino="minervameAntiNu"
var="hydrogen"
for tag in CV2ElasticNuwroSFn #CV2ElasticNuwroSFz CV2Elastic
do
  for lam in 200 #390
  do
    #for iteration in 3 4 5 6 7 8 9 10
    for iteration in 4 
    do
      inputdir="../rootfiles/minervameAntiNu/cuts${cuts}/${tag}/"

      #func="./plot_nu_eventrate_regions_neutron_proton_reweights"
      input_eff="${inputdir}/EffPurity_CombinedPlaylists.root"
      input_eff_fk="${inputdir}/EffPurity_FullKine_CombinedPlaylists.root"
      #input_eff="/pnfs/minerva/scratch/users/tejinc/CCQENu/CCQENuNeutron/Processing_20200904/NoMuonPCut_NoForwardCut_NoVtxCut_NoN2DCut_NewBlobCut2_ExclusiveInnerAngle/Neutron_CV2ElasticNuwroSFn2/EffPurity_MakeFlux-1_minervame5A.root"
      input_mig="${inputdir}/Migration_CombinedPlaylists.root"
      input_xs="${inputdir}/CrossSectionHist_It${iteration}_w_lambda_${w}_${lam}_CombinedPlaylists_${var}.root"
      input_xs_fk="${inputdir}/CrossSectionHist_FullKine_It${iteration}_w_lambda_${w}_${lam}_CombinedPlaylists_${var}.root"
      input_xs_2D="${inputdir}/CrossSectionHist2D_It${iteration}_w_lambda_${w}_${lam}_CombinedPlaylists_${var}.root"

      func="./plot_nu_eventrate_regions_neutron_proton_eff_mig_xs"
      cmd="${func} ${input_eff} ${input_eff_fk} ${input_mig} ${input_xs} ${input_xs_fk} ${input_xs_2D} ${iteration}"
      echo $cmd
      $cmd

      savefolder=plots/${neutrino}/cuts${cuts}/${tag}/${sample}_${w}_${lam}/
      mkdir -p $savefolder
      mv plots/*pdf $savefolder
      mv plots/FaFitHisto*.root $savefolder
      mv plots/AMatrix*.root $savefolder
      mv plots/*it*.root $savefolder
    done
  done
done
