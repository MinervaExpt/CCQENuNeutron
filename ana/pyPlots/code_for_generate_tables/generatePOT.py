import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G","minervame6H", "minervame6I", "minervame6J"]
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/CV2/"
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat_SigFunc2/CV2ElasticNuwroSFn/"
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat2_SigFunc2/CV2ElasticNuwroSFn/"

pot_data,pot_mc=0,0
for pl in playlists:
  fpath = folder+"MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_"+pl+".root"
  f = ROOT.TFile( fpath, "read")
  data_pot = f.pot.X()
  mc_pot = f.pot.Y()
  pot_data+=data_pot
  pot_mc+=mc_pot
  f.Close()

  print pl, data_pot, mc_pot

print "total: ", pot_data, pot_mc
