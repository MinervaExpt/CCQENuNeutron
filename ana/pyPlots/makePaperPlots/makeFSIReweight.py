import ROOT, PlotUtils
from ROOT import TFile, TH2D, TH1D,TGraph
from PlotUtils import MinervaStyle

MinervaStyle

fsi=PlotUtils.weight_fsi()

def generateRange(x0,x1,n):
  step = (x1-x0)/n
  return [ x0+step*i for i in range(n) ]

def FormatGraph( g, color, width, style, xtitle, ytitle ):
  g.SetLineColor(color)
  g.SetLineWidth(width)
  g.SetLineStyle(style)
  g.GetXaxis().SetTitle(xtitle)
  g.GetYaxis().SetTitle(ytitle)
  g.GetXaxis().SetTitleOffset(1.3)
  g.GetYaxis().SetTitleOffset(1.3)

  g.GetXaxis().SetTitleFont(42)
  g.GetYaxis().SetTitleFont(42)
  g.GetXaxis().SetTitleSize(0.04)
  g.GetYaxis().SetTitleSize(0.04)


x0,x1,n = 25, 500., 1000
KE = generateRange(x0,x1,n)

#hp1 = TH1D("proton_toNoFSI","Proton to No FSI",n,x0,x1)
#hn1 = TH1D("neutron_toNoFSI","Neutron to No FSI",n,x0,x1)
#hp2 = TH1D("proton_toOtherFSI","Proton to Other FSI",n,x0,x1)
#hn2 = TH1D("neutron_toOtherFSI","Neutron to Other FSI",n,x0,x1)


hp1 = TGraph()
hn1 = TGraph()
hp2 = TGraph()
hn2 = TGraph()

for i, keMev in enumerate(KE): 
  ke = keMev/1000.
  proton_w = list( fsi.GetQEElasticFSIWeightFromParam( ke, 12, 2212 ) )
  neutron_w = list( fsi.GetQEElasticFSIWeightFromParam( ke, 12, 2112 ) )

  hp1.SetPoint( i,  keMev, proton_w[0] )
  hn1.SetPoint( i,  keMev, neutron_w[0] )
  hp2.SetPoint( i,  keMev, proton_w[1] )
  hn2.SetPoint( i,  keMev, neutron_w[1] )


hn1.Draw("al")
FormatGraph(hn1, ROOT.kRed, 2,1, "Kinetic Energy (MeV)", "Weight")
hn1.SetMinimum(1)
hp1.Draw("l")
FormatGraph(hp1, ROOT.kBlue, 2,1, "Kinetic Energy (MeV)", "Weight")

leg = ROOT.TLegend(0.65,0.75,0.85,0.9)
leg.AddEntry( hp1, "Proton", "l")
leg.AddEntry( hn1, "Neutron", "l")
leg.SetBorderSize(0)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw()

#hn2.Draw("l")
#FormatGraph(hn2, ROOT.kRed, 2,2, "Kinetic Energy (MeV)", "Weight")
#hp2.Draw("l")
#FormatGraph(hp2, ROOT.kBlue, 2,2, "Kinetic Energy (MeV)", "Weight")


ROOT.c1.Update()
ROOT.c1.SetLogx()
ROOT.c1.SetLogy()
ROOT.c1.SetGrid()
ROOT.c1.Print("plots/FSIReweight.pdf")
