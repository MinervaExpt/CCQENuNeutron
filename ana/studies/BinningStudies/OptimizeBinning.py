import ROOT
from ROOT import TFile, TH1D, TCanvas
from ROOT.TMath import *
from PlotUtils import MnvH2D, MnvH1D
import argparse
import array

ROOT.gROOT.SetBatch()
print "Starting"
def GroupBinsWithError( errThreshold, stat_bins ):
  ret = []
  statErr = []
  
  #0. Reset error , accum, bin count
  #1. loop through each bin, accumulate stat
  #2. If stat satisfy threshold, push bins
  #3. start from beginning
  Reset = True
  binGroup = []
  accum = 0
  for i, n in enumerate(stat_bins):
    if Reset:
      binGroup = []
      accum = 0
      Reset = False

    iBin = i+1
    binGroup.append( iBin )
    accum += n
    fracstatErr = 100/Sqrt(accum)
    if fracstatErr < errThreshold or i == len(stat_bins)-1:
      ret.append( binGroup)
      statErr.append(fracstatErr)
      Reset = True
  return ret,statErr
    


    

if __name__=="__main__":
  print "in main"
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", dest="in_file", default="/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervame1D/bkg_sub_hists_ptmu.root");

  args = parser.parse_args()

  f = TFile(args.in_file,"read")

  mnvh2d_data_nobck = f.Get("h_dthetaP_ptmu_data_2track_weighted_nobck_signal")
  mnvh2d_mc_nobck = f.Get("h_dthetaP_ptmu_mc_2track_weighted_nobck_signal")

  mnvh1d_dthetaP_mc_nobck = mnvh2d_mc_nobck.ProjectionX("h_dthetaP_mc_nobck")
  mnvh1d_dthetaP_data_nobck = mnvh2d_data_nobck.ProjectionX("h_dthetaP_data_nobck")
  
  th1d_dthetaP_mc_nobck_abs_tot_errors = mnvh1d_dthetaP_mc_nobck.GetTotalError()
  th1d_dthetaP_mc_nobck_frac_tot_errors = mnvh1d_dthetaP_mc_nobck.GetTotalError(True, True)

  #Define original thetaP bins
  bins_dthetaP = array.array("d",[-3.14+i*0.157 for i in range(0,41)] )
  nbins_dthetaP = len(bins_dthetaP)
  #choose rebinning algorithm
  # 1. get number of events in each bin
  # 2. find min error

  bins_nEvents = array.array("d", [ mnvh1d_dthetaP_mc_nobck.GetBinContent(i) for i in range(1,mnvh1d_dthetaP_mc_nobck.GetNbinsX()+1 ) ])
  bins_StatErr = array.array("d", [ 100/Sqrt(v) for v in bins_nEvents ] )
  print list(range(1,mnvh1d_dthetaP_mc_nobck.GetNbinsX() ))
  print(bins_StatErr)
  print(bins_nEvents)

  binRange, statErrs = GroupBinsWithError(3, bins_nEvents)
  print binRange

  for bins,err in zip(binRange, statErrs):
    ibin0 = bins[0]
    ibin1 = bins[-1]
    x0 = mnvh1d_dthetaP_mc_nobck.GetBinLowEdge(ibin0)
    x1 = mnvh1d_dthetaP_mc_nobck.GetBinLowEdge(ibin1)
    width = mnvh1d_dthetaP_mc_nobck.GetBinWidth(ibin1)
    print "Bin Range:"
    print round(x0), round(x1+width), err

 
  #rebinning = array.array([4,3,
  #print bins_dthetaP


  
  print "hahaha"
  print mnvh1d_dthetaP_mc_nobck.GetName()

  test_bins = array.array("d",[-180,-125, -80,-55,-35,-25,-15,-10,-5,0,5,10,15,25,35,55,80,125,180])
  th1d_mc = mnvh1d_dthetaP_mc_nobck.Rebin( len(test_bins)-1, "rebinned_mc", test_bins)
  th1d_data = mnvh1d_dthetaP_data_nobck.Rebin( len(test_bins)-1, "rebinned_data", test_bins)

  th1d_mc_err = th1d_mc.Clone("mc_err")
  th1d_mc_err.Reset()
  th1d_data_err = th1d_data.Clone("mc_err")
  th1d_data_err.Reset()
  for i in range(1, th1d_data.GetNbinsX()+1 ):
    data_err = th1d_data.GetBinError(i)
    data_content = th1d_data.GetBinContent(i)
    mc_err = th1d_mc.GetBinError(i)
    mc_content = th1d_mc.GetBinContent(i)

    data_ratio = data_err/data_content if data_content>0 else 0
    mc_ratio = mc_err/mc_content if mc_content>0 else 0
    th1d_data_err.SetBinContent(i, data_ratio )
    th1d_mc_err.SetBinContent(i, mc_ratio )

  
  c = TCanvas("c","c",800,800)
  th1d_mc.Scale(1,"width")
  th1d_data.Scale(1,"width")
  th1d_mc.Draw("histe")
  th1d_data.Draw("pe same")
  c.Print("new_bin.pdf(")
  th1d_data_err.Draw("p")
  th1d_mc_err.Draw("hist same")
  c.Print("new_bin.pdf)")

  f.Close()

