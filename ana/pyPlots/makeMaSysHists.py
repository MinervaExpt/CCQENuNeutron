from ROOT import TFile, TCanvas, TLegend, TH1D
import PlotUtils
from CCQENuPlot.CCQENuPlot import CCQENuPlotUtils
from CCQENuPlot.GlobalIncludes import *

savefolder="./plots/MaSystematics/"
folder_input = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_3/CV2ElasticNuwroSF/"

rawFile = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/drut1186/CCQENu_Anatuples/MuonKludge_ProtonLLR_UpdatedNeutron/MC_Merged/minervame1Dpass1/CCQENu_mc_AnaTuple_run00111100.root"


h_q2qe_list = []

def GenerateHist(name, bins):
  h = TH1D(name,name,
