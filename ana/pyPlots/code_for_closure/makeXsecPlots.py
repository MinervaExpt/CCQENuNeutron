import ROOT, PlotUtils
f_xsec = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist_It5_w_lambda_1_550_minervame6D_hydrogen.root")

f_gxe = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/GENIEXSecExtract/cmt/GENIEXSECEXTRACT_CCQENuNeutron_Q2QE_minervame6D.root")

pot = f_xsec.pot
scale = pot.Y()/pot.X()
h_mc_effcor =   f_xsec.h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor
h_mc_effcor.Scale(scale)
h_data_effcor = f_xsec.h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor
h_eff_denum =   f_xsec.h_truth_q2qe_qelike_qe_h
h_gxe_evt =     f_gxe.ds_dq2_xsec_evRate


print "h_gxe_evt", h_gxe_evt.GetEntries(),h_gxe_evt.Integral()
print "h_eff_denum", h_eff_denum.GetEntries(),h_eff_denum.Integral()

h_mc_effcor.SetLineColor( ROOT.kBlue )
h_eff_denum.SetLineColor( ROOT.kRed )
h_gxe_evt.SetLineColor( ROOT.kBlack )

c1=ROOT.TCanvas("c1","c1")
c1.SetLogx()
h_gxe_evt.Draw("hist")
h_eff_denum.Draw("same")
h_mc_effcor.Draw("same")

leg = ROOT.TLegend(.6,.5,.85,.85 )
leg.AddEntry( h_gxe_evt, "GXE Evt Rate")
leg.AddEntry( h_eff_denum, "MC Eff Denom")
leg.AddEntry( h_mc_effcor, "MC Eff Corrected")
leg.Draw()
c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/minervame6D_effcorr_gxe_evtrate.pdf")

h_mc_effcor.Divide( h_mc_effcor, h_eff_denum )
h_gxe_evt.DivideSingle( h_gxe_evt, h_eff_denum )
h_mc_effcor.GetXaxis().SetRangeUser(0.05,2)
h_mc_effcor.GetYaxis().SetRangeUser(0.95,1.05)
h_mc_effcor.Draw()
h_gxe_evt.Draw("same")

leg2 = ROOT.TLegend(.6,.65,.9,.85 )
leg2.AddEntry( h_gxe_evt, "GXE/EffDenom","l")
leg2.AddEntry( h_mc_effcor, "MC/EffDenom","l")
leg2.Draw()

c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/minervame6D_effcorr_gxe_ratio.pdf")




#=============================

#f_xsec = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist_It5_w_lambda_1_550_minervame6D_hydrogen.root")
#f_xsec = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist_It5_w_lambda_1_550_CombinedPlaylists_hydrogen.root")
h_xs_mc = f_xsec.h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section
h_xs_eff= f_xsec.h_effhist_demH_cross_section
h_xs_gxe= f_gxe.ds_dq2_xsec

h_xs_mc.Scale(1,"width")
h_xs_eff.Scale(1/scale, "width")

print h_xs_mc.Integral()
print h_xs_eff.Integral()
print h_xs_gxe.Integral()

h_xs_mc.SetLineColor( ROOT.kBlue )
h_xs_eff.SetLineColor( ROOT.kRed )
h_xs_gxe.SetLineColor( ROOT.kBlack )

c1=ROOT.TCanvas("c1","c1")
c1.SetLogx()
h_xs_gxe.Draw("hist")
h_xs_eff.Draw("histlsame")
h_xs_mc.Draw("same")

leg = ROOT.TLegend(.6,.5,.85,.85 )
leg.AddEntry( h_xs_gxe, "GXE")
leg.AddEntry( h_xs_eff, "Eff Denom XS")
leg.AddEntry( h_xs_mc, "MC XS")
leg.Draw()
c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/minervame6D_xs_evtrate.pdf")

h_xs_mc.Divide( h_xs_mc, h_xs_gxe )
h_xs_gxe.DivideSingle( h_xs_gxe, h_xs_eff )
h_xs_mc.GetXaxis().SetRangeUser(0.05,2)
h_xs_mc.GetYaxis().SetRangeUser(0.95,1.05)
h_xs_mc.Draw("hist")
h_xs_gxe.Draw("histsame")

leg2 = ROOT.TLegend(.6,.65,.9,.85 )
leg2.AddEntry( h_xs_gxe, "GXE/EffDenom","l")
leg2.AddEntry( h_xs_mc, "MC/GXE","l")
leg2.Draw()

c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/minervame6D_xs_ratio.pdf")

#============ Testing Covariance Matrix ================

h_xs_data = f_xsec.h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section
print list(h_xs_data.GetSysErrorMatricesNames() )

