import ROOT, PlotUtils
print "Loaded PlotUtils"

from ROOT import TFile, TH2D, TH1D
ROOT.gROOT.SetBatch()

from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

fpath ="/minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame6D.root" 
fpath ="/minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_combined2.root" 
f = TFile(fpath, "read")

plotUtils = CCQENuPlotUtils()

ROOT.Cintex.Enable()
plotter = PlotUtils.MnvPlotter
plotter.SetROOT6Palette(55)

histoMap2D = dict()


histoMap2D["h_q2qe_recoil"] =                 plotUtils.BookHistos(f, "h_q2qe_recoil", 1, 1 )
histoMap2D["h_q2qe_recoil_SF"] =              plotUtils.BookHistos(f, "h_q2qe_recoil_SF", 1, 1 )
histoMap2D["h_q2qe_recoil_SF_blob"] =         plotUtils.BookHistos(f, "h_q2qe_recoil_SF_blob", 1, 1 )
histoMap2D["h_q2qe_recoil_SF_neutron"] =      plotUtils.BookHistos(f, "h_q2qe_recoil_SF_neutron", 1, 1 )
histoMap2D["h_q2qe_recoil_SF_blob_neutron"] = plotUtils.BookHistos(f, "h_q2qe_recoil_SF_blob_neutron", 1, 1 )
histoMap2D["h_q2qe_recoil_blob_neutron"] = plotUtils.BookHistos(f, "h_q2qe_recoil_blob_neutron", 1, 1 )

categories = [kHistos.kQELike_QE_H, kHistos.kQELike_QE_OTH, kHistos.kQELike_QENot, kHistos.kQELikeNot]
categories_names = ["CCQE Hydrogen", "QELike QE Carbon", "QELike Not QE", "Not QELike"]

categories_signal = [kHistos.kQELike_QE_H, kHistos.kQELike_QE_OTH, kHistos.kQELike_2p2h, kHistos.kQELike_RES, kHistos.kQELike_DIS]
categories_background = [kHistos.kQELikeNot_SingleNeutralPion, kHistos.kQELikeNot_SingleChargedPion, kHistos.kQELikeNot_MultiPion, kHistos.kQELikeNot_NoPions]

def FormatHistos( hlist ):
  for cat in categories:
    hlist[cat].SetLineColor( colors[cat] )
    hlist[cat].SetMarkerColor( colors[cat] )
    hlist[cat].SetMarkerSize( 0.1 )
    hlist[cat].SetFillColor( colors[cat] )
    hlist[cat].GetXaxis().SetTitleOffset(1.3)
    hlist[cat].GetYaxis().SetTitleOffset(1.3)
    hlist[cat].GetXaxis().SetTitleSize(0.04)
    hlist[cat].GetYaxis().SetTitleSize(0.04)
  for cat in categories_signal:
    hlist[cat].SetFillStyle( 1001 )
  for cat in categories_background:
    hlist[cat].SetFillStyle( 3001 )




def CutFunc(x):
  offset = 0.05
  ret = 0
  if(x<0.175):
    ret=0.08
  elif(x<1.4):
    ret=0.03+0.3*x
  else:
    ret = 0.45
  return ret+offset

def CutFuncTest(x):
  offset = 0.05
  ret = 0
  #if (x<0.05):
  #  ret = 0.06
  if(x<0.3):
    ret = 0.04+0.43*x
  elif(x<1.4):
    ret=0.03+0.3*x+offset
  else:
    ret = 0.45+offset
  return ret


def CutFuncOpt(x):
  offset = 0.00
  ret = 0
  if(x<0.035):
    ret=0.08
  elif(x<1.3):
    ret=0.0737+0.18*x
  else:
    ret = 0.30
  return ret+offset


cut_graph = ROOT.TGraph()
cut_graph_test = ROOT.TGraph()
for i in range(0,1600):
  q2 = i/1000.
  cut_graph.SetPoint(i,q2, CutFunc(q2) )
  cut_graph_test.SetPoint(i,q2, CutFuncTest(q2) )
cut_graph_test.SetLineStyle(1)

def DrawCutFunction():
  pass

def Plot2D( hlist, name ):
  FormatHistos(hlist)
  c = ROOT.TCanvas("c","c")
  opt = ""
  leg = ROOT.TLegend(.60,.6,.84,.9)
  leg.SetBorderSize(0)
  #leg.SetFillStyle(0)

  for cat in categories[::-1]:
    hlist[cat].GetXaxis().SetRangeUser(0,1.6)
    hlist[cat].GetYaxis().SetRangeUser(0,.6)
    hlist[cat].GetXaxis().SetTitle("Reconstructed Q^{2} (GeV^{2})")
    hlist[cat].GetYaxis().SetTitle("Recoil Energy (GeV)")

    if cat == kHistos.kQELike_QE_H:
      hlist[cat].SetMarkerSize(0.1)


    hlist[cat].Draw(opt)
    if "same" not in opt:
      opt+="same"

  for i, cat in enumerate(categories):
    leg.AddEntry( hlist[cat], categories_names[i] )

  #cut_graph.Draw("l")
  cut_graph_test.Draw("l")

  leg.Draw()

  c.SaveAs("%s.eps"%name)

  h1=hlist[kHistos.kQELike_QE_H].Clone("ratio")
  h2=hlist[kHistos.kMC].Clone("mc1")
  h1.Rebin2D(10,1)
  h2.Rebin2D(10,1)

  h1.Divide(h1, h2)
  h1.GetZaxis().SetRangeUser(0,0.3)
  h1.Draw("colz")
  cut_graph.Draw("l")
  cut_graph_test.SetLineStyle(2)
  cut_graph_test.Draw("l")
  c.Print("%s_ratio.pdf"%name)


Plot2D(histoMap2D["h_q2qe_recoil"] ,                   "./plots/recoilQ2_all")
Plot2D(histoMap2D["h_q2qe_recoil_SF"] ,                "./plots/recoilQ2_all_SF")
Plot2D(histoMap2D["h_q2qe_recoil_SF_blob"] ,           "./plots/recoilQ2_all_SF_blob")
Plot2D(histoMap2D["h_q2qe_recoil_SF_blob_neutron"] ,   "./plots/recoilQ2_all_SF_blob_neutron")
Plot2D(histoMap2D["h_q2qe_recoil_SF_neutron"] ,        "./plots/recoilQ2_all_SF_neutron")
Plot2D(histoMap2D["h_q2qe_recoil_blob_neutron"] ,   "./plots/recoilQ2_all_blob_neutron")


def GetEfficiency( hlist, func, name ):
  FormatHistos(hlist)
  c = ROOT.TCanvas("c","c")
  opt = ""
  leg = ROOT.TLegend(.60,.6,.84,.9)
  leg.SetBorderSize(0)
  #leg.SetFillStyle(0)

  signal = hlist[kHistos.kQELike_QE_H].Clone("signal")
  mc = hlist[kHistos.kMC].Clone("mc")
  signal.Rebin2D(10,1)
  mc.Rebin2D(10,1)

  hEff = PlotUtils.MnvH1D( signal.ProjectionX("hEfficiency") )
  hEff.Reset()

  hPurity = hEff.Clone("hPurity")

  print "Calculating"

  for xbin in range(1, hEff.GetNbinsX()+1):
    Q2 = hEff.GetBinCenter(xbin)
    recoil_lim = func(Q2)
    ybin_lim = signal.GetYaxis().FindBin( recoil_lim )
    print Q2, recoil_lim, ybin_lim

    signal_cut = 0
    signal_sum = 0
    mc_cut = 0


    signal_sum = signal.Integral(xbin,xbin,1, signal.GetNbinsY() )
    signal_cut = signal.Integral(xbin,xbin,1, ybin_lim )
    mc_cut = mc.Integral(xbin,xbin,1, ybin_lim )

    pur = signal_cut/mc_cut
    eff = signal_cut/signal_sum
    hEff.SetBinContent( xbin, eff )
    hPurity.SetBinContent(xbin, pur )

  hEff.GetXaxis().SetTitle("Reconstructed Q^{2}_{QE} (GeV)")
  hEff.GetYaxis().SetTitle("Eff/Purity")

  hEff.SetLineColor(ROOT.kRed)
  hPurity.SetLineColor(ROOT.kBlue)
  hProduct = hPurity.Clone("product")
  hProduct.Multiply( hEff, hPurity)
  hProduct.SetLineColor(ROOT.kBlack)

  hEff.SetFillStyle(0)
  hProduct.SetFillStyle(0)
  hPurity.SetFillStyle(0)

  hEff.GetYaxis().SetRangeUser(0.01,1.01)
  hEff.GetXaxis().SetRangeUser(0.001,1.6)
  c = ROOT.TCanvas("c","c")
  hEff.Draw("hist")
  hPurity.Draw("histsame")
  hProduct.Draw("histsame")

  leg = ROOT.TLegend(.6,.6,.85,.8)
  leg.SetBorderSize(0)
  leg.AddEntry( hEff, "Efficiency","l")
  leg.AddEntry( hPurity, "Purity","l")
  leg.AddEntry( hProduct, "Efficiency x Purity","l")
  leg.Draw()

  c.SetLogx()
  c.SetLogy()
  c.SetGrid()
  c.Print(name+".pdf")
  return (hEff,hPurity,hProduct)


def GetEfficiencyFull( hlist, name, rbx,rby ):
  FormatHistos(hlist)
  c = ROOT.TCanvas("c","c")
  opt = ""
  leg = ROOT.TLegend(.60,.6,.84,.9)
  leg.SetBorderSize(0)
  #leg.SetFillStyle(0)

  signal = hlist[kHistos.kQELike_QE_H].Clone("signal")
  mc = hlist[kHistos.kQELike].Clone("mc")
  signal.Rebin2D(rbx,rby)
  mc.Rebin2D(rbx,rby)

  hEff = signal.Clone("hEff")
  hEff.Reset()
  hPurity = signal.Clone("hPurity")

  nbinsX = hEff.GetNbinsX()
  nbinsY = hEff.GetNbinsY()

  for xbin in range(1, nbinsX+1):
    Q2 = hEff.GetXaxis().GetBinCenter(xbin)

    signal_sum = signal.Integral(xbin,xbin,1, nbinsY)
    signal_cut = 0  #total number of events within the recoil cut bin
    mc_cut = 0
    for ybin in range(1, nbinsY+1 ):
      signal_cut= signal.Integral(xbin,xbin,1,ybin)
      mc_cut = mc.Integral(xbin,xbin,1,ybin)
      signal_content = signal.GetBinContent(xbin,ybin)

      eff_bin = signal_cut/signal_sum if signal_sum != 0 else 0.0000
      pur_bin = signal_cut/mc_cut if mc_cut != 0 else 0.0000

      offset1 = 0 if pur_bin>0 else 0.00001
      offset2 = 0 if eff_bin>0 else 0.00001
      hEff.SetBinContent(xbin, ybin, eff_bin*(signal_content>0))
      hPurity.SetBinContent(xbin, ybin, pur_bin*(signal_content>0))
      if Q2> 1 and ybin < 10:
        print Q2,ybin, hPurity.GetBinContent(xbin,ybin)

  hEff.GetXaxis().SetTitle("Reconstructed Q^{2}_{QE} (GeV)")
  hEff.GetYaxis().SetTitle("Eff")
  hPurity.GetXaxis().SetTitle("Reconstructed Q^{2}_{QE} (GeV)")
  hPurity.GetYaxis().SetTitle("Purity")
  hProduct = hEff.Clone("Product")
  hProduct.Reset()
  hProduct.Multiply( hEff, hPurity )

  hEff.GetZaxis().SetRangeUser(0,1)
  hPurity.GetZaxis().SetRangeUser(0,1)
  hProduct.GetZaxis().SetRangeUser(0,1)


  hEff.GetYaxis().SetRangeUser(0,0.8)
  hPurity.GetYaxis().SetRangeUser(0,0.8)
  hProduct.GetYaxis().SetRangeUser(0,0.8)

  c = ROOT.TCanvas("c","c")

  hEff.Draw("colz")
  cut_graph.Draw("l")
  cut_graph_test.Draw("l")
  c.Print(name+"_2D_Eff.pdf")

  hPurity.Draw("colz")
  cut_graph.Draw("l")
  cut_graph_test.Draw("l")
  c.SetLogz()
  c.Print(name+"_2D_Purity.pdf")

  hProduct.Draw("colz")
  cut_graph.Draw("l")
  cut_graph_test.Draw("l")
  c.Print(name+"_2D_Product.pdf")


  return (hEff, hPurity, hProduct)

GetEfficiencyFull(histoMap2D["h_q2qe_recoil_blob_neutron"], "plots/RecoilSearch", 10,1)

effpur_test=GetEfficiency( histoMap2D["h_q2qe_recoil_blob_neutron"], CutFuncTest, "plots/EffPur_Test")
effpur_curr=GetEfficiency( histoMap2D["h_q2qe_recoil_blob_neutron"], CutFunc, "plots/EffPur_Current")


effpur_test[2].Divide( effpur_test[2], effpur_curr[2] )
effpur_test[2].Draw("hist")
ROOT.c1.Print("plots/EffPur_Ratio.pdf")
    
