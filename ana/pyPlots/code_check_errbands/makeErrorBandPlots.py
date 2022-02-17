from ROOT import TFile
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

ROOT.gROOT.SetBatch()

savefolder="./plots/"

tune="CV2ElasticNuwroSFn/"
cut0="cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/"
cut1="cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_NoSigFunc/"
path_base="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/"
fname="MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand-Angles_CombinedPlaylists.root"

f0 =  TFile(path_base+cut0+tune+fname)
f1 =  TFile(path_base+cut1+tune+fname)

f0.ls()

plotUtils = CCQENuPlotUtils()

histos=dict()
histos["h_q2qe_non99"] = ( plotUtils.BookHistos(f0, "h_q2qe_non99",1,1),plotUtils.BookHistos(f1, "h_q2qe_non99",1,1) )
histos["h_q2qe_region_99"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_99",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_99",1,1) )

histos["h_q2qe_region_00"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_00",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_00",1,1) )
histos["h_q2qe_region_01"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_01",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_01",1,1) )
histos["h_q2qe_region_02"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_02",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_02",1,1) )
histos["h_q2qe_region_03"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_03",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_03",1,1) )
histos["h_q2qe_region_04"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_04",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_04",1,1) )
histos["h_q2qe_region_05"] = ( plotUtils.BookHistos(f0, "h_q2qe_region_05",1,1),plotUtils.BookHistos(f1, "h_q2qe_region_05",1,1) )

def compareSys(h0,h1, name, canvas, doVert=True, savepath=""):
  sys0 = h0.GetVertErrorBand(name) if doVert else h0.GetLatErrorBand(name)
  sys1 = h1.GetVertErrorBand(name) if doVert else h1.GetLatErrorBand(name)
  nUniv = sys0.GetNHists()

  canvas.cd()

  sys0.DivideSingle(sys0,h0)
  sys1.DivideSingle(sys1,h1)

  sys0.GetYaxis().SetRangeUser(-2,2)
  sys0.Draw("histl")

  sys1.SetLineWidth(3)

  sys1.Draw("histlsame")

  nUniv = 0
  for eb in [sys0,sys1]:
    hists = eb.GetHists()
    nUniv=len(hists)
    for h in hists:
      h.Draw("histlsame")

  latex = ROOT.TLatex()
  latex.DrawLatexNDC(0.1,0.95, name+" nHist: %d"%nUniv )
  canvas.SetLogx()
  canvas.Print(savepath)
  

for hid in plotUtils.QELikeNotPionHistos+plotUtils.QELikeHistos:
  c = ROOT.TCanvas("c","c")
  savename="region_99/test_99_%s.pdf"%names[hid]

  c.Print(savename+"[")
  h0=histos["h_q2qe_region_99"][0][hid]
  h1=histos["h_q2qe_region_99"][1][hid]

  vertNames = h0.GetVertErrorBandNames()
  for name in vertNames:
    compareSys( h0, h1, name, c, True, savename )
  c.Print(savename+"]")

for hid in plotUtils.QELikeNotPionHistos+plotUtils.QELikeHistos:
  c = ROOT.TCanvas("c","c")
  savename="region_non99/test_non99_%s.pdf"%names[hid]

  c.Print(savename+"[")
  h0=histos["h_q2qe_non99"][0][hid]
  h1=histos["h_q2qe_non99"][1][hid]

  vertNames = h0.GetVertErrorBandNames()
  for name in vertNames:
    compareSys( h0, h1, name, c, True, savename )
  c.Print(savename+"]")



for hid in plotUtils.QELikeNotPionHistos+plotUtils.QELikeHistos:
  c = ROOT.TCanvas("c","c")
  savename="region_non99sum/test_non99_sum_%s.pdf"%names[hid]

  c.Print(savename+"[")
  h0=histos["h_q2qe_region_00"][0][hid]
  h1=histos["h_q2qe_region_00"][1][hid]
  for region in ["01","02","03","04","05"]:
    h0.Add(histos["h_q2qe_region_%s"%region][0][hid] )
    h1.Add(histos["h_q2qe_region_%s"%region][1][hid] )



  vertNames = h0.GetVertErrorBandNames()
  for name in vertNames:
    compareSys( h0, h1, name, c, True, savename )
  c.Print(savename+"]")
