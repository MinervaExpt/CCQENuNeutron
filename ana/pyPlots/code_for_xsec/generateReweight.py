import ROOT, PlotUtils
from ROOT import TFile, TChain

ROOT.gROOT.SetBatch()

chain = TChain("Truth")
chain.Add("/pnfs/minerva/persistent/users/drut1186/CCQENu_Anatuples/MuonKludge_ProtonLLR_UpdatedNeutron/MC_Merged/minervame6Dpass1/*root")
print "Chain Entries: ", chain.GetEntries()

f = TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut/CV2ElasticNuwroSF/CrossSectionHist_It4_w_lambda_1_550_minervame6D_hydrogen.root")

h_q2qe_full = f.h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section.GetCVHistoWithStatError()
h_q2qe_full.Reset()
h_q2qe_full.SetName("h_q2qe_full")

h_q2qe_cut = h_q2qe_full.Clone("h_q2qe_cut")




chain.Draw("mc_Q2/1000/1000>>h_q2qe_full", "mc_intType==1 && mc_targetA==1 && mc_current==1 && mc_charm==0 && mc_incoming==-14 && truth_is_fiducial ")
chain.Draw("mc_Q2/1000/1000>>h_q2qe_cut", "mc_intType==1 && mc_targetA==1 && mc_current==1 && mc_charm==0 && truth_is_fiducial && \
    mc_incoming==-14 && mc_primFSLepton[3]/1000 > 1.5 && mc_primFSLepton[3]/1000<20 && \
    TMath::ACos((mc_primFSLepton[0]*0+mc_primFSLepton[1]*-0.058836002+mc_primFSLepton[2]*0.99826766)/(  TMath::Sqrt(mc_primFSLepton[0]*mc_primFSLepton[0]+mc_primFSLepton[1]*mc_primFSLepton[1]+mc_primFSLepton[2]*mc_primFSLepton[2]  ) ) )<20/180*TMath::Pi() " )

c1 = ROOT.TCanvas("c1","c1")
c1.SetLogx()

h_q2qe_full.Draw("hist")
print h_q2qe_full.Integral()
c1.Print("h_q2qe_cut.pdf")

h_q2qe_cut.Divide(h_q2qe_full)
h_q2qe_cut.Draw("hist")
c1.Print("ratio.pdf")

