
#part="Neutron"
#for mode in Signal SideBand
#do
#  nohup python plot_maker.py $part $mode /minerva/data/users/tejinc/CCQENu/antinu/hists/test_neutronSel_minervame_antinu_$mode.root neutron_plots_3/$mode/ > logs/plot_neutron_$mode.log &
#done

#part="Proton"
#for mode in Signal MichelSideBand BlobSideBand MicBlobSideBand SideBand
##for mode in Signal
#do
#  nohup python plot_maker.py $part $mode /minerva/data/users/tejinc/CCQENu/antinu/hists/test_protonSel_minervame1A_neutMode_$mode.root proton_plots_3/$mode/ > logs/plot_proton_$mode.log & 
#done



part="ProtonBK"
for mode in Signal
#for mode in Signal
do
  #nohup python plot_maker.py $part $mode /minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/out.root   proton_sub_plots_3/$mode/ > logs/plot_proton_sub_$mode.log & 
  cmd="python plot_maker.py $part $mode /minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervame1D/bkg_sub_hists_ptmu_optimized_bins.root   proton_sub_plots_4/$mode/"
  echo $cmd
  #python plot_maker.py $part $mode /minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervame1D/bkg_sub_hists_ptmu_optimized_bins.root   proton_sub_plots_4/$mode/# > logs/plot_proton_sub_$mode.log & 

done

