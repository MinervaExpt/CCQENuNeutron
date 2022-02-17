import ROOT
from ROOT import TFile, TH1D, TH2D
import os
import glob
import re
import array

cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_newBlobCut-1"
cuts="_muonPCut-0_maxClusECut-0_vtxCut-0_blobMaxECut-0_n2DCut-0_nBlobCut-1_newBlobCut-1"
cuts="_baseCut_chooseTotalE_cutBlobE_new2"
cuts="_baseCut_fullStat_SigFunc2"
inputfolder="../../rootfiles/minervameAntiNu/cuts%s/CV2ElasticNuwroSFn/"%cuts
fname_signal="%s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal-Angles_CombinedPlaylists.root"%inputfolder
fname_blob="%s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand-Angles_CombinedPlaylists.root"%inputfolder


##------ proton mode
#cuts="_proton"
#inputfolder="../../rootfiles/minervameNu/cuts%s/CV2ElasticNuwroSF/"%cuts
#fname_signal="%s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal-Angles_CombinedPlaylists.root"%inputfolder
#fname_blob="%s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand-Angles_CombinedPlaylists.root"%inputfolder


#files = glob.glob("%s/test/test_*noMP_noNorm_combineSP_blob99.root"%inputfolder);
#files = glob.glob("%s/test/test_*_qe_2p2h_res_snp_scp_noNorm_blob99.root"%inputfolder);
#files = glob.glob("%s/test/test_*_qe_2p2h_snp_scp_noNorm_blob99.root"%inputfolder);
#files = glob.glob("%s/test/test_*_qe_2p2h_scp_snp_mp_noNorm_blob99.root"%inputfolder);
files = sorted(glob.glob("%s/test/test_*_qe_res_2p2h_scp_snp_noNorm_blob99.root"%inputfolder), key=lambda fpath:float(re.findall(r"[-+]?\d*\.\d+|\d+",fpath)[-4]) )
print(files)

f_signal = TFile.Open( fname_signal, "READ")
f_blob = TFile.Open( fname_blob, "READ")

pot = f_signal.Get("pot")
names=["qelike_qe_h", "qelike_qe_oth", "qelike_2p2h", "qelike_res", "qelike_dis","qelike_oth", "qelikenot_singlechargedpion", "qelikenot_singleneutralpion","qelikenot_multipion","qelikenot_nopions"]

CATS=["qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_singlechargedpion", "qelikenot_singleneutralpion"]
catConvert = {
    "qelike_qe_h":"qelike_qe_h",
    "qelike_qe_oth":"qelike_qe_oth",
    "qelike_res":"qelike_res",
    "qelike_dis":"qelike_dis",
    "qelike_oth":"qelike_oth",
    "qelike_2p2h":"qelike_2p2h",
    "qelikenot_singlechargedpion":"qelikenot_scp",
    "qelikenot_singleneutralpion":"qelikenot_snp",
    "qelikenot_multipion":"qelikenot_mp",
    "qelikenot_nopions":"qelikenot_np"}

CATS2=["qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_singlechargedpion", "qelikenot_singleneutralpion"]

region_99_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_99_%s"%name)) for name in names }
h_region_99_data = f_signal.Get("h_q2qe_region_99_data")

region_01_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_01_%s"%name)) for name in names }
h_region_01_data = f_signal.Get("h_q2qe_region_01_data")

region_02_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_02_%s"%name)) for name in names }
h_region_02_data = f_signal.Get("h_q2qe_region_02_data")

region_03_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_03_%s"%name)) for name in names }
h_region_03_data = f_signal.Get("h_q2qe_region_03_data")

region_04_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_04_%s"%name)) for name in names }
h_region_04_data = f_signal.Get("h_q2qe_region_04_data")

region_05_hists_map = { name:ROOT.TH1D(f_signal.Get("h_q2qe_region_05_%s"%name)) for name in names }
h_region_05_data = f_signal.Get("h_q2qe_region_05_data")


print ("Acquired angle hists")

ndata_99,ndata_01, ndata_02, ndata_03, ndata_04, ndata_05 = h_region_99_data.Integral(), h_region_01_data.Integral(),h_region_02_data.Integral(), h_region_03_data.Integral(), h_region_04_data.Integral(), h_region_05_data.Integral()

def CalcChi2( hdata, hmc ):
  nbins = hdata.GetNbinsX()
  scale = pot.X()/pot.Y()
  chi2 = 0
  ndf=0
  for i in range(nbins):
    ibin = i+1
    bc=hdata.GetBinCenter(ibin)
    #if( bc < 0.025 or bc >6) :
    if( bc < 0.05 or bc >1.4) :
    #if( bc < 0.025 or bc >5.9) :
      continue
      #pass
    ndf+=1
    data = hdata.GetBinContent( ibin )
    mc = hmc.GetBinContent( ibin )*scale
    dchi2 = (mc-data)**2/data if data>0 else 0
    chi2+=dchi2
  return chi2/ndf

def GetNDF( hdata ):
  nbins = hdata.GetNbinsX()
  ndf=0
  for i in range(nbins):
    ibin = i+1
    bc=hdata.GetBinCenter(ibin)
    if( bc < 0.025 or bc >6) :
    #if( bc < 0.025 or bc >5.9) :
      continue
      #pass
    ndf+=1
  return ndf



def GetMC( hists_map ):
  ret = hists_map["qelike_qe_h"].Clone("MC")
  ret.Reset()
  for name in  names:
    ret.Add( hists_map[name] )
  return ret



def GetFittedMC( hists_map, fits ):
  ret = hists_map["qelike_qe_h"].Clone("FittedMC")
  ret.Reset()
  #for name in ["qelike_qe_h", "qelike_res", "qelike_dis", "qelike_oth", "qelikenot_multipion", "qelikenot_nopions"]:
  #for name in ["qelike_qe_h", "qelike_res","qelike_dis", "qelike_oth", "qelikenot_nopions"]:
  #for name in ["qelike_qe_h", "qelike_dis", "qelike_oth", "qelikenot_multipion", "qelikenot_nopions"]:
  for name in names:
    if name in CATS: 
      continue
    if name == ["qelike_qe_h"]:
      hists_map[name].Scale(1.08)
    ret.Add( hists_map[name] )

  for i in range( len(CATS) ):
    cat = CATS[i]
    h=hists_map[cat].Clone("%s_fit"%cat)
    h.Multiply( fits[i] )
    ret.Add(h)


  return ret


wArray  =  array.array("f", [0, 0.1,0.3,0.5,0.7,0.9,1,2,3,4,5])
lamArray = array.array("f", [-1,1,10,50,100,500,1000,5000,10000,50000,100000,1000000,10000000,50000000,55000000,60000000,62000000,65000000,70000000,75000000,1000000000])
chi2Hist = ROOT.TH2D("lamvsw","lamvsw", len(wArray)-1, wArray, len(lamArray)-1, lamArray )
chi2LamGraph = ROOT.TGraph2D()
chi2LamGraph.SetName("lamvswGraph")
chi2LamGraph.GetXaxis().SetTitle("w")
chi2LamGraph.GetYaxis().SetTitle("lam")

chi2LamGraph1D = ROOT.TGraph()
chi2LamGraph1D.SetName("lamvswGraph1D")
chi2LamGraph1D.GetXaxis().SetTitle("lam")
chi2LamGraph1D.GetYaxis().SetTitle("w")

out = ROOT.TFile("wlam2.root","recreate")
out.cd()
ROOT.TH1D(h_region_99_data).Write("hdata")
minw,minlam,minchi2 = -1,-1,999999999999


hmc_Nofit = GetMC( region_99_hists_map )
hmc_Nofit.Write("mc_99_orig")

print ("Start Loop")
for i,fpath in enumerate(files):
  f = TFile.Open( fpath ,"READ");

  par=re.findall(r"[-+]?\d*\.\d+|\d+",fpath);
  w, lam = float(par[-5]), float(par[-4])
  #if ( w > 1.5 or w < 0.5 ): continue
  #if (lam < 50000000 or lam > 70000000 ): continue
  print(w,lam)

  fits = []
  #for cat in ["qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_scp", "qelikenot_snp"]:
  for cat in CATS:
    name = catConvert[cat]
    h=f.Get("hs_weights_yvar_bgType_%s"%name)
    fits.append(h)

  #print(fits)

  hdata99 = ROOT.TH1D(h_region_99_data.GetCVHistoWithError() )
  hmc99 =  GetFittedMC( region_99_hists_map, fits )
  chi2_99 = CalcChi2( hdata99, hmc99 )

  hdata01 = ROOT.TH1D(h_region_01_data.GetCVHistoWithError() )
  hmc01 =  GetFittedMC( region_01_hists_map, fits )
  chi2_01 = CalcChi2( hdata01, hmc01 )

  hdata02 = ROOT.TH1D(h_region_02_data.GetCVHistoWithError() )
  hmc02 =  GetFittedMC( region_02_hists_map, fits )
  chi2_02 = CalcChi2( hdata02, hmc02 )

  hdata03 = ROOT.TH1D(h_region_03_data.GetCVHistoWithError() )
  hmc03 =  GetFittedMC( region_03_hists_map, fits )
  chi2_03 = CalcChi2( hdata03, hmc03 )

  hdata04 = ROOT.TH1D(h_region_04_data.GetCVHistoWithError() )
  hmc04 =  GetFittedMC( region_04_hists_map, fits )
  chi2_04 = CalcChi2( hdata04, hmc04 )

  hdata05 = ROOT.TH1D(h_region_05_data.GetCVHistoWithError() )
  hmc05 =  GetFittedMC( region_05_hists_map, fits )
  chi2_05 = CalcChi2( hdata05, hmc05 )

  chi2=( ndata_02*chi2_02 + ndata_03*chi2_03 + ndata_04*chi2_04 )/(ndata_02+ndata_03+ndata_04);

  NDF=GetNDF( hdata01 )

  print(w,lam,chi2, NDF,"ind chi2:", chi2_01, chi2_02, chi2_03, chi2_04, chi2_05, chi2_99)
  if chi2< minchi2:
    minchi2 = chi2
    minw = w
    minlam=lam
  chi2Hist.Fill( w-0.01, lam-0.01, chi2 )
  chi2LamGraph.SetPoint(i,w,lam,chi2)
  chi2LamGraph1D.SetPoint(i,lam,chi2)

  out.cd()
  hmc99.Write("mc_99_%d_%.1f"%(w,lam ) )
  hmc02.Write("mc_02_%d_%.1f"%(w,lam ) )
  hmc03.Write("mc_03_%d_%.1f"%(w,lam ) )
  hmc04.Write("mc_04_%d_%.1f"%(w,lam ) )
  out.cd()
  for i in range( len(fits) ):
    fit = fits[i].GetCVHistoWithStatError()
    cat = CATS[i]
    fit.Write("fit_%s_%d_%.1f"%(cat, w,lam ))

  f.Close()

print("min: ", minw, minlam, minchi2)

out.cd()
chi2Hist.Write()
chi2LamGraph.Write()
chi2LamGraph1D.Write()
out.Close()




