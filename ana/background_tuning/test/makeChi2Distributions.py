import ROOT
import PlotUtils
import glob
import array
import os
import re

basedir="../../rootfiles/minervameAntiNu/CV2ElasticNuwroSF/test/"

tfile_angle=ROOT.TFile.Open( "../../rootfiles/minervameAntiNu/CV2ElasticNuwroSF/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand-Angles_CombinedPlaylists.root")

files = glob.glob(basedir+"*test*w_lambda*noNorm.root")

bin_lambda= array.array("f", [ -1, 0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,70,100,200,300,500,700,1000] )
nbin_lambda = len(bin_lambda)-1
bin_alpha= array.array("f", [0,1,5,10,15,20,25,30,35,40,45,50,60,70,100,200,300,500,700,1000,1500,1800,2000,2500,3000,3500,4000,4500,5000] )
nbin_alpha = len(bin_alpha)-1

hist = ROOT.TH2D("lambda_alpha","lambda_alpha", nbin_alpha, bin_alpha, nbin_lambda, bin_lambda );


#open region 99



print( len(files) )
for fname in files:
  print(fname)
  basename = os.path.basename(fname)
  basename_numbers=re.findall('(\d+)_',basename)
  parameters = [ float(num) for num in basename_numbers ]
  print( parameters )
  tfile_weight = ROOt.TFile.Open( fname ,"read" )
  #weights


  tfile_weight.Close()

tfile_angle.Close()
