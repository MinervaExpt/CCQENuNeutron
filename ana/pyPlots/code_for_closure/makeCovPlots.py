import ROOT, PlotUtils
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPalette( 53 )
f_xsec = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist_It5_w_lambda_1_550_CombinedPlaylists_hydrogen.root")

h_xs_data = f_xsec.h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section
unfoldingCov = h_xs_data.GetSysErrorMatrix("unfoldingCov", False, False )
unfoldingCov*=10e80

c1 = ROOT.TCanvas("c1","c1")
#unfoldingCov.GetZaxis().SetRangeUser(-8,8)
unfoldingCovHist = ROOT.TH2D( unfoldingCov )
unfoldingCovHist.GetZaxis().SetRangeUser(-8,8)
unfoldingCovHist.Draw("colz")
c1.Print("/minerva/data/users/tejinc/Documents/Neutrons/figures/covariance/unfoldingCov.pdf")
print h_xs_data.GetNbinsX()
print unfoldingCov.GetNrows()
print list(h_xs_data.GetSysErrorMatricesNames() )
print list(h_xs_data.GetCovMatricesNames() )

