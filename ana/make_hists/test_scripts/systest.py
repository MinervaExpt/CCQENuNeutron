from ROOT import TFile, TCanvas, TLegend
import PlotUtils as pu
from PlotUtils import MnvH2D

f = TFile("../test_muon_ela_SF_neutronW.root", "read")

h_enu_q2qe_mc = f.Get("h_enu_q2qe_mc") 

print list( h_enu_q2qe_mc.GetLatErrorBandNames () )
print list( h_enu_q2qe_mc.GetVertErrorBandNames () )


c = TCanvas("c1","c1")
lat = h_enu_q2qe_mc.GetLatErrorBand("Muon_Energy_MINOS")
hists = lat.GetHists()

leg = TLegend(.8,.8,1,1)
leg.AddEntry( hists[0].ProjectionY(), "cv")
leg.AddEntry( hists[1].ProjectionY(), "err")

hists[0].ProjectionY().Draw("hist")
hists[1].ProjectionY().Draw("histsame")
c.SetLogx()
leg.Draw()
c.Print("test_q2_1.pdf")

hists[0].ProjectionX().Draw("hist")
hists[1].ProjectionX().Draw("histsame")
c.SetLogx()
leg.Draw()
c.Print("test_enu_1.pdf")


vert = h_enu_q2qe_mc.GetVertErrorBand("Reweight_Neutron")
hists = vert.GetHists()

leg = TLegend(.8,.8,1,1)
leg.AddEntry( hists[0].ProjectionY(), "cv")
leg.AddEntry( hists[1].ProjectionY(), "err")

hists[0].ProjectionY().Draw("hist")
hists[1].ProjectionY().Draw("histsame")
c.SetLogx()
leg.Draw()
c.Print("test_neutron_q2.pdf")

hists[0].ProjectionX().Draw("hist")
hists[1].ProjectionX().Draw("histsame")
c.SetLogx()
leg.Draw()
c.Print("test_neutron_enu.pdf")
