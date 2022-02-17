#!/bin/sh

tag=CV2
SampleType=Proton
#for tag in CV2 CV pion_rpa;
#for tag in CV2 CV pion_rpa CV2Elastic;
for tag in CV CV2 CV2Elastic pion_rpa default;
do
  #playlist=minervameAntiNu
  #neutrino=0
  #multiplicity=1
  #nPlaylist=7
  playlist=minervameNu
  neutrino=1
  multiplicity=2
  nPlaylist=12 #12 for Nu, 7 for antinu
  #for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand;
  for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand;
  do
    mkdir -p plots/${playlist}/${tag}/${sample}/
    mkdir -p plots/${playlist}/${tag}/${sample}/
    func="./plot_nu_eventrate_2D_neutron_proton_reweights"
    inputorig="../rootfiles/${playlist}/${tag}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_CombinedPlaylists.root"
    inputregion="../rootfiles/${playlist}/${tag}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_CombinedPlaylists.root"
    signalweight="../rootfiles/minervameNu/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_noRes_chi2.root"
    #cmd="./plot_nu_eventrate_regions_neutron_proton_reweights ../rootfiles/${playlist}/${tag}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}-Angles_CombinedPlaylists.root ../rootfiles/minervameNu/${tag}/bkgd_weights_minervameNu_q2qe_flux_2D_optimized_Type.root ../rootfiles/minervameNu/${tag}/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_noRes_chi2.root    ${sample} ${neutrino} ${nPlaylist} "
    cmd="${func} ${inputorig} ${signalweight}  ${sample} ${neutrino} ${nPlaylist} "
    echo $cmd
    $cmd
    mv plots/*pdf plots/${playlist}/${tag}/${sample}/;
    #mv plots/*png plots/${playlist}/${tag}/${sample}/;
  done
done
