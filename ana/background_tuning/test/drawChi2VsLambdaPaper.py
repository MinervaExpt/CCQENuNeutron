import PlotUtils
import ROOT
ROOT.gROOT.SetBatch()

f = ROOT.TFile("wlam.root")
graph = f.Get("lamvswGraph1D")

graph.GetXaxis().SetTitle("#lambda")
graph.GetXaxis().SetTitleOffset(1.3)
graph.GetXaxis().CenterTitle()
graph.GetXaxis().SetTitleSize(0.04)

graph.GetYaxis().SetTitle("Weighted #chi^{2}/DOF (DOF=9)")
graph.GetYaxis().SetTitleOffset(1.3)
graph.GetYaxis().SetTitleSize(0.04)
graph.GetYaxis().CenterTitle()


graph.GetYaxis().SetRangeUser(1.2,2)

c1=ROOT.TCanvas("c1","c1")
graph.SetLineWidth(2)
graph.Draw("al")
c1.SetLogx()

c1.SetGrid()
c1.Print("lcurve.pdf")


