from ROOT import TFile
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

ROOT.gROOT.SetBatch()

savefolder="./plots/"

tune="CV2ElasticNuwroSFn/"
cut0="cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/"
cut0="cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_NoSigFunc/"
cut0="cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/"
cut0="cuts_baseCut_fullStat_SigFunc2/"
path_base="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/"
fname="CrossSectionHist_It5_w_lambda_1_550_CombinedPlaylists_hydrogen.root"

f0 =  TFile(path_base+cut0+tune+fname)

num= f0.h_q2qe_angle_10_qelike_qe_h
denom=f0.h_truth_q2qe_qelike_qe_h


def compareSys(h0, name, canvas, doVert=True, savepath=""):
  sys0 = h0.GetVertErrorBand(name) if doVert else h0.GetLatErrorBand(name)
  nUniv = sys0.GetNHists()

  canvas.cd()

  sys0.DivideSingle(sys0,h0)

  sys0.GetYaxis().SetRangeUser(.5,1.5)
  sys0.Draw("histl")


  nUniv = 0
  for eb in [sys0]:
    hists = eb.GetHists()
    nUniv=len(hists)
    for h in hists:
      h.Draw("histlsame")

  latex = ROOT.TLatex()
  latex.DrawLatexNDC(0.1,0.95, name+" nHist: %d"%nUniv )
  canvas.SetLogx()
  canvas.Print(savepath)
  

c = ROOT.TCanvas("c","c")

savename="eff_num.pdf"
c.Print(savename+"[")
vertNames = num.GetVertErrorBandNames()
for name in vertNames:
  compareSys( num, name, c, True, savename )
c.Print(savename+"]")

savename="eff_denom.pdf"
c.Print(savename+"[")
vertNames = denom.GetVertErrorBandNames()
for name in vertNames:
  compareSys( denom, name, c, True, savename )
c.Print(savename+"]")

savename="data_nobck.pdf"
data=f0.h_q2qe_region_00_data_nobck_hydrogen
c.Print(savename+"[")
vertNames = data.GetVertErrorBandNames()
for name in vertNames:
  compareSys( data, name, c, True, savename )
c.Print(savename+"]")

savename="data_nobck_unfold.pdf"
data=f0.h_q2qe_region_00_data_nobck_hydrogen_unfold
c.Print(savename+"[")
vertNames = data.GetVertErrorBandNames()
for name in vertNames:
  compareSys( data, name, c, True, savename )
c.Print(savename+"]")

savename="data_nobck_unfold_effcor.pdf"
data=f0.h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor
c.Print(savename+"[")
vertNames = data.GetVertErrorBandNames()
for name in vertNames:
  compareSys( data, name, c, True, savename )
c.Print(savename+"]")


