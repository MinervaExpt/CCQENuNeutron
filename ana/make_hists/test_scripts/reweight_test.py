import ROOT
from ROOT import TFile, TCanvas, TLegend
import PlotUtils 
from PlotUtils import MnvH1D, MnvH2D, MnvPlotter

f = TFile("/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy/CV2ElasticNuwroSFn/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_1_550_qe_res_2p2h_snp_scp.root", "read")

f2 = TFile("/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root")

plotter = MnvPlotter()
plotter.ApplyStyle( 6 )
plotter.legend_n_columns=2
#f.ls()
fitNames=["qelike_qe_oth", "qelike_2p2h", "qelike_res", "qelikenot_snp", "qelikenot_scp"]
errorgroups = plotter.error_summary_group_map["Others"]
for name in fitNames:
  w = f.Get("hs_weights_yvar_bgType_%s"%name)
  c = TCanvas("c1","c1",1200,900 )

  latNames = list( w.GetLatErrorBandNames() )
  vertNames = list( w.GetVertErrorBandNames() )

  asfrac = False
  plotter.DrawErrorSummary( w , "TL", True, True, 0.00001, False, "Others", asfrac);
  c.SetLogx()
  c.Print("fit_err_%s.pdf"%name)



leg = TLegend(.8,.5,1,1)

name="qelike_2p2h"
h = f2.Get("h_q2qe_region_00_%s"%name)
w = f.Get("hs_weights_yvar_bgType_%s"%name)

h.Draw("hist")
leg.AddEntry(h,"cv")

for i,err in enumerate(errorgroups):
  if not h.HasVertErrorBand(err): continue
  for idx in [0,1]:
    hists = h.GetVertErrorBand(err).GetHists()
    hists[idx].SetLineColor( i*2)
    hists[idx].SetLineStyle(1+idx)
    leg.AddEntry(hists[idx], "%s %d"%(err,idx) )
    hists[idx].Draw("histsame")

leg.Draw()
c.Print("h_err.pdf")

h.Multiply(h, w)
for i,err in enumerate(errorgroups):
  if not h.HasVertErrorBand(err): continue
  for idx in [0,1]:
    hists = h.GetVertErrorBand(err).GetHists()
    hists[idx].SetLineColor( i*2)
    hists[idx].SetLineStyle(1+idx)
    hists[idx].Draw("histsame")

leg.Draw()
c.Print("h_fit.pdf")

