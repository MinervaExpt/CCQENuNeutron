import ROOT, PlotUtils
print "Loaded PlotUtils"
ROOT.gROOT.SetBatch()

from ROOT import TFile, TH2D, TH1D

from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

angleLimit=60

ROOT.Cintex.Enable()

plotter = PlotUtils.MnvPlotter

plotter.SetROOT6Palette(55)

fpath ="/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root"
f = TFile(fpath, "read")

plotUtils = CCQENuPlotUtils()

histoMap2D = dict()


histoMap2D["h_dthetaP_dthetaR"] =                 plotUtils.BookHistos(f, "h_dthetaP_dthetaR", 1, 1 )

categories = [kHistos.kQELike_QE_H, kHistos.kQELike_QE_OTH, kHistos.kQELike_QENot, kHistos.kQELikeNot]
categories_names = ["CCQE Hydrogen", "QELike QE Carbon", "QELike Not QE", "Not QELike"]

categories_signal = [kHistos.kQELike_QE_H, kHistos.kQELike_QE_OTH, kHistos.kQELike_2p2h, kHistos.kQELike_RES, kHistos.kQELike_DIS]
categories_background = [kHistos.kQELikeNot_SingleNeutralPion, kHistos.kQELikeNot_SingleChargedPion, kHistos.kQELikeNot_MultiPion, kHistos.kQELikeNot_NoPions]

categories_other = [kHistos.kMC, kHistos.kData]

def FormatHistos( hlist ):
  for cat in categories+categories_signal+categories_background+categories_other:
    hlist[cat].SetLineColor( colors[cat] )
    hlist[cat].SetMarkerColor( colors[cat] )
    hlist[cat].SetMarkerSize( 0.1 )
    hlist[cat].SetFillColor( colors[cat] )
    hlist[cat].GetXaxis().SetTitleOffset(1.3)
    hlist[cat].GetXaxis().SetTitleSize(0.04)

    hlist[cat].GetYaxis().SetTitleOffset(1.3)
    hlist[cat].GetYaxis().SetTitleSize(0.04)

    hlist[cat].GetZaxis().SetTitleOffset(1.15)
    hlist[cat].GetZaxis().SetTitleSize(0.04)

    hlist[cat].GetXaxis().SetRangeUser(-angleLimit,angleLimit)
    hlist[cat].GetYaxis().SetRangeUser(-angleLimit,angleLimit)
    hlist[cat].GetXaxis().SetTitle("#delta#theta_{P} (degree)")
    hlist[cat].GetYaxis().SetTitle("#delta#theta_{R} (degree)")
  for cat in categories_signal:
    hlist[cat].SetFillStyle( 1001 )
  for cat in categories_background:
    hlist[cat].SetFillStyle( 3001 )


def FormatText(t,align,font,size,color):
  t.SetTextAlign(align)
  t.SetTextFont(font)
  t.SetTextColor(color)
  t.SetTextSize(size)

def DrawRegions(c):
  c.cd()

  linewidth=3


  #C = [ROOT.kRed, ROOT.kRed+1, ROOT.kRed+2, ROOT.kGray+1, ROOT.kGray, ROOT.kBlack]
  C = [ROOT.kGray]*6
  b0 = ROOT.TBox(-10,-10,10,10 )
  b0.SetLineColor(C[0])
  b0.SetLineWidth(linewidth)
  b0.SetFillStyle(0)

  b1 = ROOT.TBox(-20,-20,20,20 )
  b1.SetLineColor(C[1])
  b1.SetLineWidth(linewidth)
  b1.SetFillStyle(0)

  b2 = ROOT.TBox(-30,-30,30,30 )
  b2.SetLineColor(C[2])
  b2.SetLineWidth(linewidth)
  b2.SetFillStyle(0)

  b3 = ROOT.TBox(-20,-40,20,-30 )
  b3.SetLineColor(C[3])
  b3.SetLineWidth(linewidth)
  b3.SetFillStyle(0)

  b4 = ROOT.TBox(-55,-40,55,-30 )
  b4.SetLineColor(C[4])
  b4.SetLineWidth(linewidth)
  b4.SetFillStyle(0)

  b5 = ROOT.TBox(-55,-55,55,-40 )
  b5.SetLineColor(C[5])
  b5.SetLineWidth(linewidth)
  b5.SetFillStyle(0)
  b5.Draw()
  b4.Draw()
  b3.Draw()
  b2.Draw()
  b1.Draw()
  b0.Draw()

  t0=ROOT.TText(5,0,"0")
  FormatText(t0,22,42,0.06,C[0])
  t0.Draw()
  t1=ROOT.TText(15,-12,"1")
  FormatText(t1,22,42,0.06,C[1])
  t1.Draw()
  t2=ROOT.TText(25,-22,"2")
  FormatText(t2,22,42,0.06,C[2])
  t2.Draw()
  t3=ROOT.TText(15,-34,"3")
  FormatText(t3,22,42,0.05,C[3])
  t3.Draw()
  t41=ROOT.TText(40,-34,"4")
  FormatText(t41,22,42,0.05,C[4])
  t41.Draw()
  t42=ROOT.TText(-40,-34,"4")
  FormatText(t42,22,42,0.05,C[4])
  t42.Draw()

  t5=ROOT.TText(40,-47,"5")
  FormatText(t5,22,42,0.06,C[5])
  t5.Draw()

  t99=ROOT.TText(40,30,"99")
  FormatText(t99,22,42,0.06,ROOT.kWhite)
  t99.Draw()


  return [b0,b1,b2,b3,b5,b4,t0,t1,t2,t3,t41,t42,t5,t99]



def Plot2D( hlist, name, cat, title="" ):
  FormatHistos(hlist)
  c = ROOT.TCanvas("c","c")
  z1=ROOT.TPad("z1","z1", 0.1, 0.1, 0.9, 0.9)
  z2=ROOT.TPad("z2","z2", 0.1, 0.1, 0.9, 0.9)
  z1.Draw()
  z2.Draw()
  z2.SetFillStyle(4000)
  opt = ""
  leg = ROOT.TLegend(.60,.6,.84,.9)
  leg.SetBorderSize(0)
  #leg.SetFillStyle(0)

  total = hlist[kHistos.kMC]
  
  comp = hlist[cat].Clone("ratio")
  comp.Divide(comp,total)
  comp.GetZaxis().SetRangeUser(0,.7)
  comp.GetZaxis().SetTitle("Fraction of Events")

  #z1.cd()
  #comp.Draw("cont4z")
  comp.Draw("colz")

  #z2.cd()
  ##a=comp.Clone("axis")
  ##a.Reset()
  ##a.SetFillStyle(4000)
  ##a.Draw("axis")
  #a=126
  #b=-10
  #z2.Range(-a,-a+b,a,a+b)

  boxes=DrawRegions(c)

  t=ROOT.TLatex()
  t.SetTextColor(ROOT.kWhite)
  t.SetTextSize(0.05)
  t.SetTextAlign(22)
  t.DrawLatexNDC(0.5,0.85,title)
  

  c.Print( name )



#Plot2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_H.pdf", kHistos.kQELike_QE_H,  "CCQE Hydrogen")
#Plot2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_C.pdf", kHistos.kQELike_QE_OTH, "QELike QE Carbon")
#Plot2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_2p2h.pdf", kHistos.kQELike_2p2h, "QELike 2p2h")
#Plot2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_RES.pdf", kHistos.kQELike_RES, "QELike RES")
#Plot2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_QELikeNot.pdf", kHistos.kQELikeNot, "Pions and Others")




def Print2D( hlist, name, cat, title="" ):

  total = hlist[kHistos.kMC]
  
  comp = hlist[cat].Clone("ratio")
  comp.Divide(comp,total)
  comp.GetZaxis().SetRangeUser(0,.7)
  comp.GetZaxis().SetTitle("Fraction of Events")

  r0 = comp.GetXaxis().FindBin(-60);
  r1 = comp.GetXaxis().FindBin(60);

  print title

  label=",".join( [ str(comp.GetYaxis().GetBinCenter(r)) for r in range(r0,r1+1)] )
  print label

  print "{"
  for ybin in range(r0,r1+1):
    xvalues = [ str(comp.GetBinContent( xbin, ybin )) for xbin in range(r0,r1+1) ]
    line=",".join(xvalues)
    print "{ %s },"%line
  print "}"

  print "\n"

Print2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_H.pdf", kHistos.kQELike_QE_H,  "CCQE Hydrogen")
Print2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_C.pdf", kHistos.kQELike_QE_OTH,  "CCQE Carbon")
Print2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_C.pdf", kHistos.kQELike_QENot,  "QELike NotQE")
Print2D(histoMap2D["h_dthetaP_dthetaR"] ,    "plots/angular_ratio_C.pdf", kHistos.kQELikeNot,  "QELikeNot")
