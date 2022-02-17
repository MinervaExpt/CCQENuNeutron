import ROOT, PlotUtils
ROOT.gROOT.SetBatch()

path="../../rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_SigFunc/CV2ElasticNuwroSFn/"
#../rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_SigFunc/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root

for pl in ["5A","6A","6B", "6C", "6D","6E","6F","6G"]:
  mig="Migration_minervame%s.root"%pl
  evt="MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_minervame%s.root"%pl

  fmig=ROOT.TFile(path+mig)
  fevt=ROOT.TFile(path+evt)

  print pl, fmig.h_q2qe_angle_10_reco_qelike_qe_h.GetEntries(), fevt.h_enu_q2qe_qelike_qe_h.GetEntries()
