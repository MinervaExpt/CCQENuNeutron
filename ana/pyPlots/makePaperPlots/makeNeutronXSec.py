import ROOT, PlotUtils
from ROOT import TFile, TH2D, TH1D,TGraph
from PlotUtils import MinervaStyle

MinervaStyle

path="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/hadronReweight/"
hydrogen = TFile(path+"cross_section_neutron_hydrogen_hd.root")
carbon_new = TFile(path+"new_geant4_xsec.root")
carbon = TFile(path+"cross_section_neutron_carbon_hd.root")

hydrogen.ls()
carbon_new.ls()
carbon.ls()

def FormatXs(total, ela,inel,color):
  for g in [total,ela,inel]:
    g.SetLineColor(color)
    g.SetLineWidth(4)
  total.SetLineStyle(1)
  ela.SetLineStyle(2)
  inel.SetLineStyle(3)

c=ROOT.TCanvas("c","c")

c.SetLogy()
c.SetLogx()
c.SetGrid()

h_total, h_elastic, h_inelastic = hydrogen.totalXS, hydrogen.elXS, hydrogen.inelXS
c_total, c_elastic, c_inelastic = carbon_new.total_xsec, carbon_new.elastic_xsec, carbon_new.inelastic_xsec



FormatXs( h_total, h_elastic, h_inelastic, ROOT.kRed )
FormatXs( c_total, c_elastic, c_inelastic, ROOT.kBlue)

#hydrogen.totalXS.SetLineColor(ROOT.kRed)
#hydrogen.inelXS.SetLineColor(ROOT.kRed)
#hydrogen.elXS.SetLineColor(ROOT.kRed)
#
#carbon_new.total_xsec.SetLineColor(ROOT.kBlue)
#carbon_new.inelastic_xsec.SetLineColor(ROOT.kBlue)
#carbon_new.elastic_xsec.SetLineColor(ROOT.kBlue)

h_total.GetXaxis().SetRangeUser(0.1,1000)
h_total.GetYaxis().SetRangeUser(5,100000)

h_total.Draw("Al")
h_elastic.Draw("l")
h_inelastic.Draw("l")

c_total.Draw("l")
c_elastic.Draw("l")
c_inelastic.Draw("l")

leg = ROOT.TLegend(0.5,0.6,0.84,0.9)
leg.SetBorderSize(0)
leg.AddEntry( h_total, "Hydrogen Total XS", "l")
leg.AddEntry( h_elastic, "Hydrogen Elastic XS", "l")
leg.AddEntry( h_inelastic, "Hydrogen Inelastic XS", "l")

leg.AddEntry( c_total, "Carbon Total XS", "l")
leg.AddEntry( c_elastic, "Carbon Elastic XS", "l")
leg.AddEntry( c_inelastic, "Carbon Inelastic XS", "l")

leg.Draw()
c.Print("Neutron_XSec.pdf")

