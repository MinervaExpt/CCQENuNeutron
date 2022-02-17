import PlotUtils
import ROOT
import array

ROOT.gROOT.SetBatch()
print "Set Batch"

plotter=PlotUtils.MnvPlotter()

plotter.SetROOT6Palette(55)
plotter.axis_maximum=0.3
plotter.axis_minimum=0.0

c = ROOT.TCanvas("c","c")

print "RHC"
path="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg-14-minervame6A.root"
f=ROOT.TFile(path, 'read')
flux = f.Get("flux_E_cvweighted")
flux.PopVertErrorBand("Flux")
flux.GetXaxis().SetRangeUser(0,22)
flux.GetXaxis().SetTitle("Energy (GeV)")
flux.GetYaxis().SetRangeUser(0,0.3)
plotter.DrawErrorSummary(flux,"TR", False)
c.Print("plots/flux_rhc_numubar_sys.pdf")



print "FHC"
path="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg14-minervame1D1M1NWeightedAve.root"
f=ROOT.TFile(path, 'read')
flux = f.Get("flux_E_cvweighted")
flux.PopVertErrorBand("Flux")
flux.GetXaxis().SetRangeUser(0,22)
flux.GetXaxis().SetTitle("Energy (GeV)")
flux.GetYaxis().SetRangeUser(0,0.3)
plotter.DrawErrorSummary(flux,"TR", False)
c.Print("plots/flux_fhc_numu_sys.pdf")


