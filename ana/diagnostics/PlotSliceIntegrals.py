import PlotUtils
import os,sys
import ROOT

ROOT.TH1.AddDirectory(False)

inputfile = sys.argv[1]
var = sys.argv[2]


myfile = ROOT.TFile(inputfile)
myhist = myfile.Get("h_%s_ptmu_mc_nobck_unfold_effcor_cross_section"%(var))

y_bins = myhist.GetNbinsY()
x_bins = myhist.GetNbinsX()
print y_bins,x_bins
for i in range(0,y_bins+1):
    myproj = myhist.ProjectionX("p_%d"%(i),i,i)
    myintegral_flow = myproj.Integral(0,-1)
    myintegral_noflow = myproj.Integral(1,x_bins)
    print "Pt bin %d\tintegral_flow=%e\tintegral_noflow=%e"%(i,myintegral_flow,myintegral_noflow)

