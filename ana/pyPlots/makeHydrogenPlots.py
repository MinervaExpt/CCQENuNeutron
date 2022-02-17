from ROOT import TFile
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

savefolder="./plots/"

f =  TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/studies/DphiTTStudies/root_files/pndphitt_Combined.root")

plotUtils = CCQENuPlotUtils()

histoMap2D = dict()
histoMapDphiTT = dict()
histoMapEnu = dict()

histoMap2D["h_enuTrue_dphitt"] = plotUtils.BookHistos(f, "h_enuTrue_dphitt",2,2)
histoMap2D["h_enuCCQE_dphitt"] = plotUtils.BookHistos(f, "h_enuCCQE_dphitt",2,2)
histoMap2D["h_enuFit_dphitt"] = plotUtils.BookHistos(f, "h_enuFit_dphitt",2,2)
histoMap2D["h_enuNoFit_dphitt"] = plotUtils.BookHistos(f, "h_enuNoFit_dphitt",2,2)
histoMap2D["h_denuCCQE_dphitt"] = plotUtils.BookHistos(f, "h_denuCCQE_dphitt",2,2)
histoMap2D["h_denuFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFit_dphitt",2,2)
histoMap2D["h_denuNoFit_dphitt"] = plotUtils.BookHistos(f, "h_denuNoFit_dphitt",2,2)
histoMap2D["h_denuFracCCQE_dphitt"] = plotUtils.BookHistos(f, "h_denuFracCCQE_dphitt",2,2)
histoMap2D["h_denuFracFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFracFit_dphitt",2,2)
histoMap2D["h_denuFracNoFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFracNoFit_dphitt",2,2)


for name, histos in histoMap2D.items():
  xmin,xmax=1,histos[0].GetNbinsX()
  ymin,ymax=1,histos[0].GetNbinsY()

  histoMapDphiTT[name] = plotUtils.Project2DHistos(histos, "y", xmin,xmax)
  histoMapEnu[name] = plotUtils.Project2DHistos(histos, "x",ymin,ymax)

  print histos[kHistos.kMC].Integral()
  c = ROOT.TCanvas("c1","c1")
  plotUtils.DrawStacked1D(c, histoMapEnu[name],None,opt="hist",xtitle="Reco E_{#nu} (GeV)",ytitle="Evt Rate /Bin", ymax=8000)
  c.Print("%s/%s_enu.pdf"%(savefolder,name) )

  plotUtils.DrawStacked1D(c, histoMapDphiTT[name],None,opt="hist",xtitle="Reco #delta #phi_{TT} (degree)",ytitle="Evt Rate /Bin")
  c.Print("%s/%s_dphitt.pdf"%(savefolder,name) )

for ih in [kHistos.kQELike_QE_OTH, kHistos.kQELike_2p2h, kHistos.kQELike_RES, kHistos.kQELike]:
  histos = []
  histos2= []
  hnames = ["EnuTrue", "EnuCCQE", "EnuFit", "EnuMuProton"]
  htitles= ["h_enuTrue_dphitt", "h_enuCCQE_dphitt", "h_enuFit_dphitt", "h_enuNoFit_dphitt" ]
  c = ROOT.TCanvas("c1","c1")
  for hname in htitles:
    histos.append( histoMapEnu[hname][ih] )
    histos2.append( histoMapDphiTT[hname][ih] )

  plotUtils.FormatHistosLine( histos )
  plotUtils.FormatHistosLine( histos2 )
  plotUtils.DrawHistos1D(c,histos,hnames,"E_{#nu} (GeV)","Evt Rate", "hist")
  c.SetGrid()
  c.Print("%s/EnuType_%s.pdf"%(savefolder, names[ih]) )

  plotUtils.DrawHistos1D(c,histos2,hnames,"#delta #phi_{TT} (degree)","Evt Rate", "hist")
  c.SetGrid()
  c.Print("%s/DphiTTType_%s.pdf"%(savefolder, names[ih]) )

  #do Ratio
  print histos
  h0 = histos[0].Clone("true")
  print h0
  histos2 = [ h.Clone( h.GetName()+"_clone") for h in histos ]
  for i in range(len(histos2)):
    histos2[i].Divide(histos2[i],h0)
  print histos2
  plotUtils.DrawHistos1D(c,histos2,hnames,"E_{#nu} (GeV)","Ratio to True", "hist")
  c.SetGrid()
  c.Print("%s/EnuTypeRatio_%s.pdf"%(savefolder, names[ih]) )

  

  hnames = ["d"+name for name in hnames]
  htitles= [ title.replace("enu","denu") for title in htitles ]
  c = ROOT.TCanvas("c1","c1")
  histos = [None]
  histosProfile = [None]
  for hname in htitles[1:]:
    histos.append( histoMapEnu[hname][ih] )
    h = histoMap2D[hname][ih].Clone( histoMap2D[hname][ih].GetName()+"_clone")
    h.Rebin2D(1,2)
    histosProfile.append( h.ProfileY() )

  plotUtils.FormatHistosLine( histos )
  plotUtils.FormatHistosLine( histosProfile )

  plotUtils.DrawHistos1D(c,histos[1:],hnames[1:],"#delta E_{#nu} (GeV)","Evt Rate", "hist", doStat=True)
  c.SetGrid()
  c.Print("%s/DEnuType_%s.pdf"%(savefolder, names[ih]) )

  plotUtils.DrawHistos1D(c,histosProfile[1:],hnames[1:],"#delta #phi_{TT} (degree)","E (GeV)", "hist", doStat=False)
  c.SetGrid()
  c.Print("%s/DEnuTypeProfile_%s.pdf"%(savefolder, names[ih]) )

c = ROOT.TCanvas("c1","c1")
histoMap2D["h_denuCCQE_dphitt"][kHistos.kQELike].Draw("colz")
c.Print("plots/2D.pdf")
