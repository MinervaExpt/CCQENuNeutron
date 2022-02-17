import PlotUtils
import ROOT

ROOT.gROOT.SetBatch()

path1 = "../test_MuonEventSelection_flux.root"
path2 = "../test_MuonEventSelection_flux_2.root"

f1 = ROOT.TFile(path1,"read")
f2 = ROOT.TFile(path2,"read")

histoname = "vtx_energy"
#histoname = "pnAngle_ptmu"
evtType="qe"

h1=f1.Get("h_%s_qelike_%s_notpn"%(histoname, evtType))
h2=f2.Get("h_%s_qelike_%s_notpn"%(histoname, evtType))
h22=f2.Get("h_%s_qelike_%s_pn"%(histoname, evtType))

#h1.ProjectionX("h1p").Draw()
#h2.ProjectionX("h2p").Draw("histsame")
h22.Draw("")
h2.Draw("histsame")

ROOT.c1.Print("test.pdf")
