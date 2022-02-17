import PlotUtils
import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervame5A"]
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat2_SigFunc2/CV2ElasticNuwroSFn/"

f_xsec = ROOT.TFile(folder+"CrossSectionHist_It4_w_lambda_1_200_CombinedPlaylists_hydrogen.root")

f_xsec.ls()

h_xsec = f_xsec.Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section")
h_xsec.Scale(1,"width")


r=range(1,h_xsec.GetNbinsX()+1)
value_xsec = [h_xsec.GetBinContent(i)*1e40 for i in r]
total_xsec=sum(value_xsec)

for name in h_xsec.GetVertErrorBandNames():

  errband = h_xsec.GetVertErrorBand(name)
  herr = errband.GetErrorBand(True)
  #herr.Fit("pol0")

  errs = [herr.GetBinContent(i)*100 for i in r]   

  err= sum([ e*x for e,x in zip(errs,value_xsec) ])/total_xsec
  
  print name,err
  #l=",".join(map("%.3f",errs))
  #l=",".join(map(lambda x:"%.3f"%x,errs))
  #print l
  

for name in h_xsec.GetLatErrorBandNames():

  errband = h_xsec.GetLatErrorBand(name)
  herr = errband.GetErrorBand(True)
  #herr.Fit("pol0")

  errs = [herr.GetBinContent(i)*100 for i in r]   

  err= sum([ e*x for e,x in zip(errs,value_xsec) ])/total_xsec
  
  print name,err
  #l=",".join(map("%.3f",errs))
  #l=",".join(map(lambda x:"%.3f"%x,errs))
  #print l
 
