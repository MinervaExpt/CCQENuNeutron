import ROOT, PlotUtils
from ROOT import TFile, TChain
import array

ROOT.gROOT.SetBatch()

f = TFile("/minerva/app/users/tejinc/Generators/processing/output/genieAntinuHFixed_10000_.gtrac.root");
tree = f.gRooTracker

q2binsEdgeList=[0, 0.00625, 0.0125, 0.025, 0.0375, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 2.0, 4.0, 6.0,10.0]

q2binsEdge=array.array('f',q2binsEdgeList)
histo = ROOT.TH1D("q2","q2", len(q2binsEdge)-1,q2binsEdge)

for evt in tree:
  q2=0
  pm = list(evt.StdHepP4)
  q2+= (pm[0]-pm[2*4+0])**2
  q2+= (pm[1]-pm[2*4+1])**2
  q2+= (pm[2]-pm[2*4+2])**2
  q2+= -(pm[3]-pm[2*4+3])**2
  histo.Fill(q2)

nEntries= tree.GetEntries()
histo.Scale(0.6497159/nEntries,"width");
histo.Draw("histl")
ROOT.c1.SetLogx()
ROOT.c1.SetGrid()
ROOT.c1.Print("ccantinu_evt_rate.pdf")


import ROOT
ROOT.gROOT.SetBatch(True)
import os, array, math
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
def getspline(material):
    tf = ROOT.TFile("/minerva/app/users/kleykamp/cmtuser/Minerva_v22r1p1/MParamFiles/data/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root")
    splinename = None
    if material == "hydrogen":
        splinename = "nu_mu_H1/tot_cc"
    if material == "carbon":
        splinename = "nu_mu_C12/tot_cc"
    if material == "lead":
        splinename = "nu_mu_Pb207/tot_cc"
    if material == "iron":
        splinename = "nu_mu_Fe56/tot_cc"
    if material == "oxygen":
        splinename = "nu_mu_O16/tot_cc"
    assert splinename != None, "Wasn't able to get a spline name from material name %s" % material
    out = tf.Get(splinename)
    assert out != None, "Spline is none %s for spline name %s" % (out, splinename)
    ## Convert to actual function so that integrals are smoother
    funcname = "func_" + material
    func = ROOT.TF1(funcname, "pol6", 0, 100)
    out.Fit(funcname)
    return func



def main():
    filename = "/pnfs/minerva/persistent/users/drut1186/CCQENu_Anatuples/MuonKludge_ProtonLLR_UpdatedNeutron/MC_Merged/minervame6Dpass1/CCQENu_mc_AnaTuple_run00122579.root"
    tf = ROOT.TFile(filename)
    hist = tf.Get("truehydrogen")
    spline = getspline("hydrogen")
    # Now we have the CC rate hist and the spline
    units = 1e-38 # converts to 1e-38cm
    units *= 1 / float(100 * 100)
    bins = range(1, hist.GetNbinsX() + 1)
    tf2 = ROOT.TFile("/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg-14-minervame6A.root"
    fluxhist = tf2.Get("flux_E_cvweighted")
    avg_n_targets = 0
    n = 0
    weighted_avg_n_targets = 0
    total_weight = 0
    xsechist = ROOT.TH1D("xsec", "xsec", 1000, 0, 100)
    ntargetshist = ROOT.TH1D("ntargets", "ntargets", 1000, 0, 100)
    pot_main = 4.94552510035e+21
    potnorm = pot_main
    for b in bins: 
        # N targets = n events / flux / xsec
        binlowedge = hist.GetBinLowEdge(b)
        binwidth = hist.GetBinWidth(b)
        # Splines go from 0 to 100 and there are 1000 bins. So 10x the binlowedge
        #spline_low = int(binlowedge * 10)
        #spline_high = int((binlowedge + binwidth) * 10)
        spline_low = binlowedge
        spline_high = binlowedge + binwidth
        xsecInt = spline.Integral(spline_low, spline_high)
        xsec = xsecInt / binwidth * units
        xsechist.SetBinContent(b, xsec)
        rate = hist.GetBinContent(b)
        rate_err = hist.GetBinError(b)
        flux = fluxhist.GetBinContent(b) # Same binning
        flux_err = fluxhist.GetBinError(b) # Same binning
        foo = xsec * flux * potnorm
        if foo == 0: 
            print "Zero foo:", b, xsec, flux
            continue
        ntargets = rate / foo
        ntargets_err = math.sqrt(rate_err**2 / foo**2 + flux_err**2 * rate / foo)
        ntargetshist.SetBinContent(b, ntargets)
        ntargetshist.SetBinError(b, ntargets_err)
        print "N targets", b, ntargets, ntargets_err
        if ntargets > 0 and b > 10 and b < 200:
            avg_n_targets += ntargets
            n += 1
            weight = 1 / ntargets_err**2
            weighted_avg_n_targets += ntargets * weight
            total_weight += weight
    avg_n_targets /= float(n)
    weighted_avg_n_targets /= float(total_weight)
    print "Avg n targets:", avg_n_targets
    print "Weighted avg n targets:", weighted_avg_n_targets
    c = ROOT.TCanvas("c", "c", 800, 800)
    xsechist.Draw("hist")
    c.Print("test/hydrogen_xsec.png")
    ntargetshist.Draw("hist")
    c.Print("test/hydrogen_ntargets.png")
main()
