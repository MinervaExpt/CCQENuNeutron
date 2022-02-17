import ROOT, PlotUtils
from ROOT import TFile,TCanvas
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPalette( 53 )

savefolder="/minerva/data/users/tejinc/Documents/Neutrons/figures/mhrw/"

f =  TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/make_hists/test_Mig.root")
f_mhrw =  TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/make_hists/test_Mig_neutronCVRW.root")


h_q2qe_reco = f.h_q2qe_angle_10_reco_qelike_qe_h
h_q2qe_reco_mhrw = f_mhrw.h_q2qe_angle_10_reco_qelike_qe_h

h_q2qe_truth = f.h_q2qe_angle_10_truth_qelike_qe_h
h_q2qe_truth_mhrw = f_mhrw.h_q2qe_angle_10_truth_qelike_qe_h

def FormatHisto(h,color=None,linewidth=None, scale=False):
  if color:
    h.SetMarkerColor(color)
    h.SetLineColor(color)
  if linewidth:
    h.SetLineWidth(linewidth)
  if scale:
    print "scaling"
    h.Scale(1,"width")

FormatHisto(h_q2qe_reco,ROOT.kBlack,2, False)
FormatHisto(h_q2qe_reco_mhrw,ROOT.kRed,2, False)
FormatHisto(h_q2qe_truth,ROOT.kBlack,2, False)
FormatHisto(h_q2qe_truth_mhrw,ROOT.kRed,2, False)

def DrawErrorbandCV( h, name, tag, doVert=True, doRatio = False ):
  errband = h.GetVertErrorBand( name ) if doVert else h.GetLatErrorBand( name )
  errband=PlotUtils.MnvVertErrorBand( errband ) if doVert else PlotUtils.MnvLatErrorBand( errband )
  if doRatio:
    cv = ROOT.TH1D( errband )
    errband.DivideSingle( errband, cv )
    errband.GetYaxis().SetRangeUser(0.5,1.5)
    errband.GetYaxis().SetTitle( "Ratio to CV" )

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

DrawRatio(h_q2qe_reco, h_q2qe_reco_mhrw, "mig_reco")
DrawRatio(h_q2qe_truth, h_q2qe_truth_mhrw, "mig_truth")


DrawErrorbandCV( h_q2qe_reco, "Reweight_Neutron", "noMHRW", doVert=True, doRatio = False )
DrawErrorbandCV( h_q2qe_reco, "Reweight_Neutron", "noMHRW", doVert=True, doRatio = True )

DrawErrorbandCV( h_q2qe_reco_mhrw, "Reweight_Neutron", "MHRW", doVert=True, doRatio = False )
DrawErrorbandCV( h_q2qe_reco_mhrw, "Reweight_Neutron", "MHRW", doVert=True, doRatio = True )

