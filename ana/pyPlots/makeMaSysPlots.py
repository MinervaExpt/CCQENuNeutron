from ROOT import TFile, TCanvas, TLegend
import PlotUtils
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

savefolder="./plots/MaSystematics/"
folder_input = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_3/CV2ElasticNuwroSF/"

f_fullKine =  TFile("%s/CrossSectionHist_FullKine_It4_w_lambda_1_550_hydrogen.root"%folder_input,"read")
f =  TFile("%s/CrossSectionHist_It4_w_lambda_1_550_hydrogen.root"%folder_input,"read")

plotUtils = CCQENuPlotUtils()

f.ls()

h_data_nobck_fullKine = f_fullKine.Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor")
h_mc_nobck_fullKine = f_fullKine.Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor")

h_data_nobck = f.Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor")
h_mc_nobck = f.Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor")

print list( h_data_nobck.GetVertErrorBandNames() )
print list( h_data_nobck.GetLatErrorBandNames() )


def GetMaNorm(h):
  errband = h.GetVertErrorBand( "GENIE_NormCCQE" ) 
  norms = []
  for h in errband.GetHists():
    norms.append( h.Integral()/errband.Integral() )
  return norms

def DrawErrorbandCV( h,hdata, name, tag, doVert=True, doRatio = False ):
  errband = h.GetVertErrorBand( name ) if doVert else h.GetLatErrorBand( name )
  errband=PlotUtils.MnvVertErrorBand( errband ) if doVert else PlotUtils.MnvLatErrorBand( errband )
  hdata = hdata.Clone("clone")
  norms = GetMaNorm( h )

  errband.GetHist(0).Scale(norms[0])
  errband.GetHist(1).Scale(norms[1])
  if doRatio:
    cv = ROOT.TH1D( errband )
    errband.DivideSingle( errband, cv )
    errband.GetYaxis().SetRangeUser(0,2)
    errband.GetYaxis().SetTitle( "Ratio to MnvGENIEv1" )
    hdata.Divide(hdata,h)
  else:
    errband.Scale(0.006,"width")
    hdata.Scale(0.006,"width")

  errband.SetLineWidth(2)
  errband.SetLineColor(ROOT.kBlack)

  c = TCanvas("c","c")
  errband.Draw("hist")
  for i, h in enumerate(errband.GetHists()):
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

vertNames = ["GENIE_MaCCQEshape"]
latNames = []

for name in vertNames:
  DrawErrorbandCV(h_mc_nobck, h_data_nobck, name, "comp", True)
  DrawErrorbandCV(h_mc_nobck, h_data_nobck, name, "comp", True, True)

for name in latNames:
  #DrawErrorband(h_data_nobck, name, "data", False)
  #DrawErrorband(h_data_nobck, name, "data", False, True)
  DrawErrorbandCV(h_data_nobck, name, "data", False)
  DrawErrorbandCV(h_data_nobck, name, "data", False, True)



