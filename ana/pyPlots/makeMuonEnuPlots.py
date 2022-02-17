from ROOT import TFile, TCanvas, TLegend
import PlotUtils
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

savefolder="./plots/Muon_Scale/"
folder_input = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_3/CV2ElasticNuwroSF/"

#f =  TFile("%s/CrossSectionHist_It4_w_lambda_1_550_hydrogen.root"%folder_input,"read")
f =  TFile("%s/CrossSectionHist2D_It10_w_lambda_1_550_hydrogen.root"%folder_input,"read")
f2=  TFile("%s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root"%folder_input,"read")


plotUtils = CCQENuPlotUtils()

f.ls()


h_enu_q2qe_data_nobck = f.Get("h_enu_q2qe_data_nobck_hydrogen")
h_enu_q2qe_mc_nobck = f.Get("h_enu_q2qe_qelike_qe_h_nobck_hydrogen")

h_data_nobck = h_enu_q2qe_data_nobck.ProjectionX("h_enu_data_nobck")
h_mc_nobck = h_enu_q2qe_mc_nobck.ProjectionX("h_enu_mc_nobck")

h_enu_q2qe_qelike_qe_h = f2.Get("h_enu_q2qe_qelike_qe_h")
#h_enu_qelike_qe_h = h_enu_q2qe_qelike_qe_h.ProjectionX("h_enu_mc_nobck")
h_mc_nobck = h_enu_q2qe_qelike_qe_h.ProjectionX("h_enu_mc_nobck")

print list( h_mc_nobck.GetVertErrorBandNames() )
print list( h_mc_nobck.GetLatErrorBandNames() )

def DrawErrorband( h, name, tag, doVert=True, doRatio = False ):
  errband = h.GetVertErrorBand( name ) if doVert else h.GetLatErrorBand( name )
  errband=PlotUtils.MnvVertErrorBand( errband ) if doVert else PlotUtils.MnvLatErrorBand( errband )
  if doRatio:
    cv = ROOT.TH1D( errband )
    errband.DivideSingle( errband, cv )

  c = TCanvas("c","c")
  errband.DrawAll("hist", True )

  savename = ("vert_" if doVert else "lat_") + name +"_"+tag+"_ratio-%d.pdf"%doRatio

  latex = ROOT.TLatex()
  latex.SetTextAlign(22)
  latex.DrawLatexNDC(0.5,0.95, name )


  c.SetGrid()
  c.SetLogx()
  c.Print(savefolder+savename)
  del c
  #leg = TLegend(0.7,0.7,0.85,0.85)
  #errband.Draw("histl")
  #leg.AddEntry(
  #for
def DrawErrorbandCV( h, name, tag, doVert=True, doRatio = False ):
  errband = h.GetVertErrorBand( name ) if doVert else h.GetLatErrorBand( name )
  errband=PlotUtils.MnvVertErrorBand( errband ) if doVert else PlotUtils.MnvLatErrorBand( errband )
  if doRatio:
    cv = ROOT.TH1D( errband )
    errband.DivideSingle( errband, cv )
    errband.GetYaxis().SetRangeUser(0.5,1.5)
    errband.GetYaxis().SetTitle( "Ratio to MnvGENIEv1" )

  c = TCanvas("c","c")
  errband.SetLineWidth(2)
  errband.SetLineColor(ROOT.kBlack)
  errband.Draw("hist")
  for h in errband.GetHists():
    h.Draw("histsame")

  savename = ("vert_cv_" if doVert else "lat_cv_") + name +"_"+tag+"_ratio-%d.pdf"%doRatio

  latex = ROOT.TLatex()
  latex.SetTextAlign(22)
  latex.DrawLatexNDC(0.5,0.95, name )


  c.SetGrid()
  c.SetLogx()
  c.Print(savefolder+savename)
  del c


def DrawErrorbandCV2( h,hdata, name, tag, doVert=True, doRatio = False ):
  errband = h.GetVertErrorBand( name ) if doVert else h.GetLatErrorBand( name )
  errband=PlotUtils.MnvVertErrorBand( errband ) if doVert else PlotUtils.MnvLatErrorBand( errband )
  hdata = hdata.Clone("clone")
  if doRatio:
    cv = ROOT.TH1D( errband )
    errband.DivideSingle( errband, cv )
    errband.GetYaxis().SetRangeUser(0,2)
    errband.GetYaxis().SetTitle( "Ratio to MnvGENIEv1" )


    hdata.Divide(hdata,h)
  else:
    pass
    errband.Scale(1.1)
    errband.SetMaximum( max(errband.GetMaximum(), hdata.GetMaximum() )*1.2)
    #errband.Scale(0.006,"width")
    #hdata.Scale(0.006,"width")

  errband.SetLineWidth(2)
  errband.SetLineColor(ROOT.kBlack)

  c = TCanvas("c","c")
  errband.Draw("hist")
  for h in errband.GetHists():
    h.Draw("histsame")

  hdata.Draw("same");
  savename = ("vert_cv_" if doVert else "lat_cv_") + name +"_"+tag+"_ratio-%d.pdf"%doRatio

  latex = ROOT.TLatex()
  latex.SetTextAlign(22)
  latex.DrawLatexNDC(0.5,0.95, name )


  c.SetGrid()
  c.SetLogx()
  c.Print(savefolder+savename)
  del hdata
  del errband
  del c



vertNames = ["MINOS_Reconstruction_Efficiency", "Reweight_Neutron"]
latNames = ['Muon_Energy_MINERvA', 'Muon_Energy_MINOS', 'Muon_Energy_Resolution']

for name in vertNames:
  DrawErrorbandCV(h_data_nobck, name, "data", True)
  DrawErrorbandCV(h_data_nobck, name, "data", True, True)

  DrawErrorbandCV2(h_mc_nobck, h_data_nobck, name, "comp",True )
  DrawErrorbandCV2(h_mc_nobck, h_data_nobck, name, "comp",True , True)

for name in latNames:
  #DrawErrorband(h_data_nobck, name, "data", False)
  #DrawErrorband(h_data_nobck, name, "data", False, True)
  DrawErrorbandCV(h_data_nobck, name, "data", False)
  DrawErrorbandCV(h_data_nobck, name, "data", False, True)

  DrawErrorbandCV2(h_mc_nobck, h_data_nobck, name, "comp", False)
  DrawErrorbandCV2(h_mc_nobck, h_data_nobck, name, "comp", False, True)

#histoMap2D = dict()
#histoMapDphiTT = dict()
#histoMapEnu = dict()
#
#histoMap2D["h_enuTrue_dphitt"] = plotUtils.BookHistos(f, "h_enuTrue_dphitt",2,2)
#histoMap2D["h_enuCCQE_dphitt"] = plotUtils.BookHistos(f, "h_enuCCQE_dphitt",2,2)
#histoMap2D["h_enuFit_dphitt"] = plotUtils.BookHistos(f, "h_enuFit_dphitt",2,2)
#histoMap2D["h_enuNoFit_dphitt"] = plotUtils.BookHistos(f, "h_enuNoFit_dphitt",2,2)
#histoMap2D["h_denuCCQE_dphitt"] = plotUtils.BookHistos(f, "h_denuCCQE_dphitt",2,2)
#histoMap2D["h_denuFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFit_dphitt",2,2)
#histoMap2D["h_denuNoFit_dphitt"] = plotUtils.BookHistos(f, "h_denuNoFit_dphitt",2,2)
#histoMap2D["h_denuFracCCQE_dphitt"] = plotUtils.BookHistos(f, "h_denuFracCCQE_dphitt",2,2)
#histoMap2D["h_denuFracFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFracFit_dphitt",2,2)
#histoMap2D["h_denuFracNoFit_dphitt"] = plotUtils.BookHistos(f, "h_denuFracNoFit_dphitt",2,2)
#
#
#for name, histos in histoMap2D.items():
#  xmin,xmax=1,histos[0].GetNbinsX()
#  ymin,ymax=1,histos[0].GetNbinsY()
#
#  histoMapDphiTT[name] = plotUtils.Project2DHistos(histos, "y", xmin,xmax)
#  histoMapEnu[name] = plotUtils.Project2DHistos(histos, "x",ymin,ymax)
#
#  print histos[kHistos.kMC].Integral()
#  c = ROOT.TCanvas("c1","c1")
#  plotUtils.DrawStacked1D(c, histoMapEnu[name],None,opt="hist",xtitle="Reco E_{#nu} (GeV)",ytitle="Evt Rate /Bin", ymax=8000)
#  c.Print("%s/%s_enu.pdf"%(savefolder,name) )
#
#  plotUtils.DrawStacked1D(c, histoMapDphiTT[name],None,opt="hist",xtitle="Reco #delta #phi_{TT} (degree)",ytitle="Evt Rate /Bin")
#  c.Print("%s/%s_dphitt.pdf"%(savefolder,name) )
#
#for ih in [kHistos.kQELike_QE_OTH, kHistos.kQELike_2p2h, kHistos.kQELike_RES, kHistos.kQELike]:
#  histos = []
#  histos2= []
#  hnames = ["EnuTrue", "EnuCCQE", "EnuFit", "EnuMuProton"]
#  htitles= ["h_enuTrue_dphitt", "h_enuCCQE_dphitt", "h_enuFit_dphitt", "h_enuNoFit_dphitt" ]
#  c = ROOT.TCanvas("c1","c1")
#  for hname in htitles:
#    histos.append( histoMapEnu[hname][ih] )
#    histos2.append( histoMapDphiTT[hname][ih] )
#
#  plotUtils.FormatHistosLine( histos )
#  plotUtils.FormatHistosLine( histos2 )
#  plotUtils.DrawHistos1D(c,histos,hnames,"E_{#nu} (GeV)","Evt Rate", "hist")
#  c.SetGrid()
#  c.Print("%s/EnuType_%s.pdf"%(savefolder, names[ih]) )
#
#  plotUtils.DrawHistos1D(c,histos2,hnames,"#delta #phi_{TT} (degree)","Evt Rate", "hist")
#  c.SetGrid()
#  c.Print("%s/DphiTTType_%s.pdf"%(savefolder, names[ih]) )
#
#  #do Ratio
#  print histos
#  h0 = histos[0].Clone("true")
#  print h0
#  histos2 = [ h.Clone( h.GetName()+"_clone") for h in histos ]
#  for i in range(len(histos2)):
#    histos2[i].Divide(histos2[i],h0)
#  print histos2
#  plotUtils.DrawHistos1D(c,histos2,hnames,"E_{#nu} (GeV)","Ratio to True", "hist")
#  c.SetGrid()
#  c.Print("%s/EnuTypeRatio_%s.pdf"%(savefolder, names[ih]) )
#
#  
#
#  hnames = ["d"+name for name in hnames]
#  htitles= [ title.replace("enu","denu") for title in htitles ]
#  c = ROOT.TCanvas("c1","c1")
#  histos = [None]
#  histosProfile = [None]
#  for hname in htitles[1:]:
#    histos.append( histoMapEnu[hname][ih] )
#    h = histoMap2D[hname][ih].Clone( histoMap2D[hname][ih].GetName()+"_clone")
#    h.Rebin2D(1,2)
#    histosProfile.append( h.ProfileY() )
#
#  plotUtils.FormatHistosLine( histos )
#  plotUtils.FormatHistosLine( histosProfile )
#
#  plotUtils.DrawHistos1D(c,histos[1:],hnames[1:],"#delta E_{#nu} (GeV)","Evt Rate", "hist", doStat=True)
#  c.SetGrid()
#  c.Print("%s/DEnuType_%s.pdf"%(savefolder, names[ih]) )
#
#  plotUtils.DrawHistos1D(c,histosProfile[1:],hnames[1:],"#delta #phi_{TT} (degree)","E (GeV)", "hist", doStat=False)
#  c.SetGrid()
#  c.Print("%s/DEnuTypeProfile_%s.pdf"%(savefolder, names[ih]) )
#
#c = ROOT.TCanvas("c1","c1")
#histoMap2D["h_denuCCQE_dphitt"][kHistos.kQELike].Draw("colz")
#c.Print("plots/2D.pdf")
