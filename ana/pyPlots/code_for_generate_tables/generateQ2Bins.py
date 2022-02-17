import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervame5A"]
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/CV2/"

for pl in playlists:
  fpath = folder+"MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_"+pl+".root"
  f = ROOT.TFile( fpath, "read")
  h = ROOT.h_q2qe_mc


  for i in range(h.GetNbinsX() ):
    ibin = i+1
    lowEdge = h.GetBinLowEdge(ibin)
    width = h.GetBinWidth(ibin)
    highEdge = lowEdge+width

    center = h.GetBinCenter(ibin)

    print "%d & %s & %s & %s \\\\\hline" % (ibin, lowEdge, highEdge, center)


  f.Close()


