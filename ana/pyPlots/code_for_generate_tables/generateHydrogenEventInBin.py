import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervame5A"]
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat_SigFunc2/CV2ElasticNuwroSFn/"

for pl in playlists:
  fpath = folder+"MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root"
  f = ROOT.TFile( fpath, "read")
  h = ROOT.h_q2qe_region_00_qelike_qe_h


  for i in range(h.GetNbinsX() ):
    ibin = i+1
    lowEdge = h.GetBinLowEdge(ibin)
    width = h.GetBinWidth(ibin)
    highEdge = lowEdge+width

    center = h.GetBinCenter(ibin)
    value = h.GetBinContent(ibin)
    err = h.GetBinError(ibin)
    ratio = err/value*100 if value > 0 else 100


    print "%d & %s & %s & %s & %.2f & %.2f &%.2f\\\\\hline" % (ibin, lowEdge, highEdge, center, value, err,ratio)


  f.Close()


