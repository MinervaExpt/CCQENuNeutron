doMuonSel=1
doMigration=0
doEffPurity=0
playlists="minervame5A minervame6A minervame6B minervame6C minervame6D minervame6E minervame6F minervame6G " 
fbase0="./"
#for tag in 0_20 0_40 #40_70 70_100;
#do
#  fbase="root_files/inelastic_graph_${tag}_"
#  cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
#  #$cmd
#  fbase="root_files/inelastic_graph_abf_${tag}_"
#  cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
#  #$cmd
#  fbase="root_files/inelastic_graph_abf_inela_${tag}_"
#  cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
#  $cmd
#done

#for tag in 0_45 
#do
#  fbase="root_files/inelastic_graph_drop_0d3_${tag}_"
#  cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
#  $cmd
#done
mode=1
#                  dbs_rw-1_neutroncvrw
fbase="root_files/dbs_rw-${mode}_neutroncvrw_"
cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
$cmd
