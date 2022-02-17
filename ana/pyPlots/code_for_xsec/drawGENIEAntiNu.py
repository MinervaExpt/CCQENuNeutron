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
