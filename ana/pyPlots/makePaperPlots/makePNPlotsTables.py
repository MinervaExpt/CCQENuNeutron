import ROOT, PlotUtils
print "Loaded PlotUtils"
ROOT.gROOT.SetBatch()

from ROOT import TFile, TH2D, TH1D

from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

angleLimit=60

ROOT.Cintex.Enable()

plotter = PlotUtils.MnvPlotter

plotter.SetROOT6Palette(55)

fpath ="/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW/CV2ElasticNuwroSFn/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root"

fpath="/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameNu/cuts_np_NoEB_openAngle_newCat5/CV2Elastic/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root"
f = TFile(fpath, "read")

plotUtils = CCQENuPlotUtils()

histoMap2D = dict()


histoMap2D["h_blobDist_pnAngle"] =                 plotUtils.BookHistos(f, "h_blobDist_pnAngle", 1, 1 )
#"h_pnEndDist_pnAngle","h_blobDist_pnAngle", "h_blobDist_pnEndDist" })

categories = [kHistos.kQELike_PN, kHistos.kQELike_NotPN, kHistos.kQELike_QE_NotPN, kHistos.kQELike_2p2h_NotPN, kHistos.kQELike_RES_NotPN, kHistos.kQELike_DIS_NotPN, kHistos.kQELikeNot,kHistos.kMC, kHistos.kData]
categories_names = ["QELike nN", "QELike 0N", "QELike QE 0N", "QELike 2p2h 0N", "QELike RES 0N", "QELike DIS 0N", "QELikeNot", "MC", "Data"]


print "Cat, Total , Angle, Min R, Max R, Combined"
for i in range( len(categories_names) ):
  cat = categories[i]
  name = categories_names[i]
  h = histoMap2D["h_blobDist_pnAngle"][cat]
  xb_min, xb_max = 1, h.GetNbinsX()
  yb_min, yb_max = 1, h.GetNbinsY()

  total = h.Integral()
  angle_bin = h.GetYaxis().FindBin(30)
  angle = h.Integral(0,-1, angle_bin, yb_max)
  blob1_bin = h.GetXaxis().FindBin(400)
  blob2_bin = h.GetXaxis().FindBin(1800)
  blob1 = h.Integral(blob1_bin, xb_max,0,-1)
  blob2 = h.Integral(xb_min, blob2_bin,0,-1)
  allcuts=h.Integral(blob1_bin, blob2_bin,angle_bin,yb_max)

  print "%s,%f,%f,%f,%f,%f"%(name, total, angle, blob1, blob2,allcuts)




































































































































































