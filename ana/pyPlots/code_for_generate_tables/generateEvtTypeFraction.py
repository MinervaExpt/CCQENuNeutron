import PlotUtils
import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervameAntiNu"]
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat_SigFunc2/CV2ElasticNuwroSFn/"

evtTypes=["qelike_qe_h", "qelike_qe_oth", "qelike_2p2h", "qelike_res","qelikenot_singleneutralpion","qelikenot_singlechargedpion"]

line="region & N Events &"
for tp in evtTypes:
  line+=" & %s"%tp
print line
for pl in playlists:
  fpath = folder+"MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root"
  f = ROOT.TFile( fpath, "read")
  #f.ls("h_q2qe_*_mc")
  for regions in [0,1,2,3,4,5,99]:
    hmc = f.Get("h_q2qe_region_%02d_mc"%regions)
    integral=hmc.Integral()
    hists =  [f.Get("h_q2qe_region_%02d_%s"%(regions, tp)) for tp in evtTypes ]
    #integral = hists[1].Integral()
    fracs = [ h.Integral()/integral for h in hists]
    lfracs = " & %.1f"%integral
    for frac in fracs: lfracs+="& %.2f  "%frac
    print regions, lfracs

for pl in playlists:
  fpath = folder+"MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand-Angles_CombinedPlaylists.root"
  f = ROOT.TFile( fpath, "read")
  #f.ls("h_q2qe_*_mc")
  for regions in [0,1,2,3,4,5,99]:
    hmc = f.Get("h_q2qe_region_%02d_mc"%regions)
    integral=hmc.Integral()
    hists =  [f.Get("h_q2qe_region_%02d_%s"%(regions, tp)) for tp in evtTypes ]
    #integral = hists[1].Integral()
    fracs = [ h.Integral()/integral for h in hists]
    lfracs = " & %.1f"%integral
    for frac in fracs: lfracs+="& %.2f  "%frac
    print regions, lfracs


    #print "%d & %s & %s & %s & %.2f & %.2f &%.2f\\\\\hline" % (ibin, lowEdge, highEdge, center, value, err,ratio)


  f.Close()


