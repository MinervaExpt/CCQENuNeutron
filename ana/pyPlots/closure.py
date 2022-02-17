from ROOT import TFile, TH2D, TH1D

from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

import PlotUtils
import UnfoldUtils

folder="../rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_3/CV2ElasticNuwroSF/"

pl="CombinedPlaylists"
pl="minervame6D"

f_migration = TFile(folder+"Migration_%s.root"%pl)
f_eff = TFile(folder+"EffPurity_%s.root"%pl)
f_xs = TFile(folder+"CrossSectionHist_It4_w_lambda_1_550_%s_hydrogen.root"%pl)
f_genie= TFile("/minerva/data/users/tejinc/Documents/Neutrons/figures/cross_sections/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_%s.root"%pl)


h_mc_xs = f_xs.h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section.GetCVHistoWithStatError()
h_genie_xs = f_genie.ds_dq2_xsec.GetCVHistoWithStatError()

h_mc_xs.SetLineColor(ROOT.kBlack)
h_genie_xs.SetLineColor(ROOT.kBlue)

c=ROOT.TCanvas("c","c")
c.SetLogx()
c.SetGrid()

h_mc_xs.Scale(1,"width")

h_mc_xs.GetYaxis().SetRangeUser(0,25e-39)
h_mc_xs.Draw("p")
h_genie_xs.Draw("histsame")
c.Print("plots/closure/xs_compare_ratio-0.pdf")

h_mc_xs.Divide( h_mc_xs, h_genie_xs)
h_mc_xs.GetYaxis().SetRangeUser(0.9,1.1)
h_mc_xs.GetYaxis().SetTitle("MC/MnvGENIEv1")
h_mc_xs.Draw("hist")
c.Print("plots/closure/xs_compare_ratio-1.pdf")

print h_mc_xs


#h_mc_unfolded=None
#dummyCovMatrix=None
#h_mc_nobck = f_eff.Get("h_q2qe_angle_10_qelike_qe_h")
#h_mc_truth = f_eff.Get("h_truth_q2qe_qelike_qe_h")
#
#h_q2qe_angle_10_migration = f_migration.Get("h_q2qe_angle_10_migration_qelike_qe_h") #reco vs true
#h_q2qe_angle_10_reco  = f_migration.Get("h_q2qe_angle_10_reco_qelike_qe_h") #reco vs true
#h_q2qe_angle_10_truth = f_migration.Get("h_q2qe_angle_10_truth_qelike_qe_h") #reco vs true
#h_q2qe_angle_10_proj_reco = h_q2qe_angle_10_migration.ProjectionX() #reco 
#h_q2qe_angle_10_proj_truth = h_q2qe_angle_10_migration.ProjectionY() #truth
#
#h_diff = h_q2qe_angle_10_truth.Clone("h_diff")
#h_diff.Add( h_mc_nobck, -1 )
##h_diff.Add( h_q2qe_angle_10_proj_reco, -1 )
#
#
#
#text = ROOT.TLatex()
#text.SetTextSize(0.04)
#c=ROOT.TCanvas("c1","c1")
#
##===========1 test if there is difference between migration and efficiency numerator===
#h_diff.Draw("hist")
#text.DrawLatexNDC(0.3,0.95,"migration truth - eff num ")
#c.SetLogx()
#c.Print("plots/Hydrogen/h_diff.pdf")
#
#h_diff.Divide( h_diff, h_mc_nobck )
#h_diff.Draw("hist")
#text.DrawLatexNDC(0.3,0.95,"(migration truth - eff num)/eff num ")
#c.Print("plots/Hydrogen/h_diff_frac.pdf")
#
#h_mc_nobck.SetLineColor(ROOT.kBlue)
#h_q2qe_angle_10_truth.SetLineColor(ROOT.kBlack)
#leg1 = ROOT.TLegend(.7,.8,.9,.9)
#leg1.AddEntry(h_mc_nobck,"eff num QE-H")
#leg1.AddEntry(h_q2qe_angle_10_truth,"migration truth QE-H")
#h_mc_nobck.Draw("hist")
#h_q2qe_angle_10_truth.Draw("histsame")
#leg1.Draw()
#c.Print("plots/Hydrogen/compare_hists.pdf")
#
##=========2 test if there is closure ==================================================
#
#c.SetLogx(False)
#c.SetLogy(False)
#
#unfold=UnfoldUtils.MnvUnfold()
#h_unfolded_1 = h_q2qe_angle_10_reco.Clone("h_unfolded_1")
#h_unfolded_1.Reset()
#h_unfolded_4 = h_unfolded_1.Clone("h_unfolded_4")
#dummyCovMatrix = ROOT.TMatrixD()
#
#h_reco = h_q2qe_angle_10_reco
#h_truth = h_q2qe_angle_10_truth
#h_migration=h_q2qe_angle_10_migration
#print(type(h_migration))
#
#unfold.UnfoldHisto( h_unfolded_1, dummyCovMatrix, h_migration, h_reco, 1, 1, False, False )
#unfold.UnfoldHisto( h_unfolded_4, dummyCovMatrix, h_migration, h_reco, 1, 4, False, False )
#
#
#h_reco.SetLineColor(ROOT.kBlack)
#h_truth.SetLineColor(ROOT.kBlue)
#h_unfolded_1.SetLineColor(ROOT.kBlack)
#h_unfolded_1.SetLineStyle(2)
#h_unfolded_4.SetLineColor(ROOT.kRed)
#h_unfolded_4.SetLineStyle(4)
#
#h_truth.Draw("hist")
#h_reco.Draw("histsame")
#h_unfolded_1.Draw("histsame")
#h_unfolded_4.Draw("histsame")
#
#
#
#h_mc_nobck.SetLineColor(ROOT.kGreen)
#h_mc_nobck.SetLineWidth(1)
#h_mc_nobck.Draw("histsame")
#
#c.SetLogx()
#
#leg = ROOT.TLegend(.7,.7,.9,.9)
#leg.AddEntry(h_truth,"truth","l")
#leg.AddEntry(h_reco,"reco","l")
#leg.AddEntry(h_unfolded_1,"unfolded it1","l")
#leg.AddEntry(h_unfolded_4,"unfolded it4","l")
#leg.AddEntry(h_unfolded_4,"Eff Numerator","l")
#leg.Draw()
#c.Print("plots/Hydrogen/compare_unfold.pdf")
