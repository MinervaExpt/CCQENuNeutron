import PlotUtils
from PlotUtils import MnvH1D
import ROOT

flux_file = ROOT.TFile("../../../../PlotUtils/data/flux/flux-gen2thin-pdg-14-minervame6A.root")
print flux_file.ls()

flux_orig = flux_file.flux_E_cvweighted

flux_new = ROOT.TH1D("flux_rebinned","flux_rebinned", 1000,0,100)

for i in range(1,flux_new.GetNbinsX()+1):
  E = flux_new.GetBinCenter(i)
  flux = flux_orig.GetBinContent( flux_orig.FindBin(E) )
  flux_new.SetBinContent( i, flux )

f = ROOT.TFile("flux_rhc_numubar.root", "recreate")
f.cd()
flux_new.Write("flux_rhc_numubar_rebinned")
f.Close()
