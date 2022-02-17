import ROOT, PlotUtils
ROOT.gROOT.SetBatch()

path="/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_NoSigFunc"
#path="/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW"

f_evt0=ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/%s/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root"%path)
#h_dthetaPdthetaR_q2qe_fine

f_evt = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/%s/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root"%path)

f_eff = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/%s/CV2ElasticNuwroSFn/EffPurity_CombinedPlaylists.root"%path)

f_mig = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/%s/CV2ElasticNuwroSFn/Migration_CombinedPlaylists.root"%path)

print f_evt, f_eff, f_mig

h_evt0     = f_evt0.h_dthetaPdthetaR_q2qe_fine_qelike_qe_h.ProjectionY("h_evt0")
h_evt      = f_evt.h_q2qe_region_00_qelike_qe_h
h_evt = h_evt0
h_mig_reco = f_mig.h_q2qe_angle_10_reco_qelike_qe_h
h_mig_true = f_mig.h_q2qe_angle_10_truth_qelike_qe_h
h_eff_num  = f_eff.h_q2qe_angle_10_qelike_qe_h


h_eff_num.SetLineStyle(2)
h_eff_num.SetLineColor(ROOT.kBlue)
h_mig_true.SetLineStyle(2)
h_mig_true.SetLineColor(ROOT.kRed)

h_mig_reco.SetLineColor(ROOT.kRed)
h_evt.SetLineColor(ROOT.kBlack)

h_evt.Draw("hist")
h_mig_reco.Draw("histsame")
h_mig_true.Draw("histlsame")
h_eff_num.Draw("histlsame")

leg = ROOT.TLegend(.6,.7,.9,.9)
leg.AddEntry(h_evt, "Signal Reco", "l")
leg.AddEntry(h_mig_reco, "Migration Reco", "l")
leg.AddEntry(h_mig_true, "Migration True", "l")
leg.AddEntry(h_eff_num, "Eff Numerator", "l")
leg.Draw()

ROOT.c1.SetLogx()
ROOT.c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/CombinedPlaylists_SelMigEff_evtRate_1.pdf")

#================================
print "Checking N Entries"
print "h_evt: %d"%h_evt.GetEntries()
print "h_mig_reco: %d"%h_mig_reco.GetEntries()
print "h_mig_true: %d"%h_mig_true.GetEntries()
print "h_eff_num: %d"%h_eff_num.GetEntries()

print "Checking Integral"
print "h_evt: %f"%h_evt.Integral()
print "h_mig_reco: %f"%h_mig_reco.Integral()
print "h_mig_true: %f"%h_mig_true.Integral()
print "h_eff_num: %f"%h_eff_num.Integral()



#===============================

h_evt.Divide( h_evt,h_mig_reco )
h_eff_num.Divide(h_eff_num,h_mig_true )

c1 = ROOT.TCanvas("c1","c1")
h_evt.GetYaxis().SetRangeUser(0.9,1.1)
h_evt.Draw("hist")
h_eff_num.Draw("histsame")

leg2 = ROOT.TLegend(.2,.7,.5,.9)
leg2.AddEntry(h_evt, "Signal Reco/Mig Reco", "l")
leg2.AddEntry(h_eff_num, "Eff Num/Mig True", "l")
leg2.Draw()
c1.SetGrid()
c1.SetLogx()
c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/xsec_difference/CombinedPlaylists_SelMigEff_ratio_1.pdf")
