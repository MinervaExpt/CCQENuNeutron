import ROOT, PlotUtils
from ROOT import TFile
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPalette( 53 )

savefolder="/minerva/data/users/tejinc/Documents/Neutrons/figures/mhrw/"

f =  TFile("/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist2D_It4_w_lambda_1_550_CombinedPlaylists_hydrogen.root")
f_mhrw =  TFile("/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSFn/CrossSectionHist2D_It4_w_lambda_1_550_CombinedPlaylists_hydrogen.root")


h_eff_denom = f.h_truth_q2qe_qelike_qe_h
h_eff_denom_mhrw = f_mhrw.h_truth_q2qe_qelike_qe_h

h_eff_num = f.h_q2qe_angle_10_qelike_qe_h
h_eff_num_mhrw = f_mhrw.h_q2qe_angle_10_qelike_qe_h

h_eff = h_eff_num.Clone("h_eff")
h_eff_mhrw = h_eff_num.Clone("h_eff_mhrw")

h_eff.Divide(h_eff_num, h_eff_denom)
h_eff_mhrw.Divide(h_eff_num_mhrw, h_eff_denom_mhrw)

def FormatHisto(h,color=None,linewidth=None, scale=False):
  if color:
    h.SetMarkerColor(color)
    h.SetLineColor(color)
  if linewidth:
    h.SetLineWidth(linewidth)
  if scale:
    print "scaling"
    h.Scale(1,"width")

FormatHisto(h_eff_denom,ROOT.kBlack,2, True)
FormatHisto(h_eff_denom_mhrw,ROOT.kRed,2, True)

FormatHisto(h_eff_num,ROOT.kBlack,2, True)
FormatHisto(h_eff_num_mhrw,ROOT.kRed,2, True)

FormatHisto(h_eff,ROOT.kBlack,2, False)
FormatHisto(h_eff_mhrw,ROOT.kRed,2, False)

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


def DrawRatio( h,h_mhrw, tag):
  c=ROOT.TCanvas("c1","c1")
  c.cd()
  h.Draw("histl")
  h_mhrw.Draw("histlsame")
  leg = ROOT.TLegend(.6,.6,.85,.75)
  leg.AddEntry(h, "CV2ElasticNuwroSF")
  leg.AddEntry(h_mhrw, "CV2ElasticNuwroSF+NeutronCV")
  leg.Draw()
  c.SetLogx()
  c.Print("%s/compare_eff_%s.pdf"%(savefolder,tag))


  hratio = h.Clone("ratio")
  hratio.GetYaxis().SetTitle("ratio")

  hratio.Divide(h_mhrw, h )
  leg2 = ROOT.TLegend(.6,.6,.85,.75)
  leg2.AddEntry(hratio, "Ratio without/with neutronCV")
  hratio.Draw("hist")
  leg2.Draw()
  c.Print("%s/compare_eff_%s_ratio.pdf"%(savefolder,tag))

DrawRatio(h_eff_denom, h_eff_denom_mhrw, "denom")
DrawRatio(h_eff_num, h_eff_num_mhrw, "num")
DrawRatio(h_eff, h_eff_mhrw, "eff")


