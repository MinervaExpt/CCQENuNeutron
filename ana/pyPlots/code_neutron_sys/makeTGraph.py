import ROOT
from ROOT import TGraph
ROOT.gROOT.SetBatch()

fname_new = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/hadronReweight/new_geant4_xsec.root"
fname_hd = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/hadronReweight/cross_section_neutron_carbon_hd.root"

f_new = ROOT.TFile(fname_new)
f_hd = ROOT.TFile(fname_hd)

totalXs=f_hd.totalXS
InelXs =f_hd.inelXS
ElXs =  f_hd.elXS
abf = f_hd.abfalterer

total_xsec=f_new.total_xsec
elastic_xsec=f_new.elastic_xsec
inelastic_xsec=f_new.inelastic_xsec


total_xsec_ratio=ROOT.TGraph()
elastic_xsec_ratio=ROOT.TGraph()
inelastic_xsec_ratio=ROOT.TGraph()

total_xsec_ratio2=ROOT.TGraph()
elastic_xsec_ratio2=ROOT.TGraph()
inelastic_xsec_ratio2=ROOT.TGraph()

for i in range( total_xsec.GetN()+1):
  x,y0=ROOT.Double(), ROOT.Double()
  total_xsec.GetPoint(i,x,y0)

  y1 = elastic_xsec.Eval(x)
  y2 = inelastic_xsec.Eval(x)
  y3 = abf.Eval(x)
  total_xsec_ratio.SetPoint(i,x,(y0-y3)/y3*100)
  elastic_xsec_ratio.SetPoint(i,x,(y1-y3)/y3*100)
  inelastic_xsec_ratio.SetPoint(i,x,(y2-y3)/y3*100)
  if x>550:
    break
  
  print x,y0,y1,y2,y3

for i in range( total_xsec.GetN()+1):
  x,y0=ROOT.Double(), ROOT.Double()
  totalXs.GetPoint(i,x,y0)

  y1 = ElXs.Eval(x)
  y2 = InelXs.Eval(x)
  y3 = abf.Eval(x)
  total_xsec_ratio2.SetPoint(i,x,(y0-y3)/y3*100)
  elastic_xsec_ratio2.SetPoint(i,x,(y1-y3)/y3*100)
  inelastic_xsec_ratio2.SetPoint(i,x,(y2-y3)/y3*100)
  if x>550:
    break
  
  print x,y0,y1,y2,y3



out = ROOT.TFile("out.root","recreate")
total_xsec_ratio.Write("total")
elastic_xsec_ratio.Write("elastic")
inelastic_xsec_ratio.Write("inelastic")

total_xsec_ratio2.Write("total2")
elastic_xsec_ratio2.Write("elastic2")
inelastic_xsec_ratio2.Write("inelastic2")

out.Close()


#for i in range( abf.GetN()+1):
#  x,y0=ROOT.Double(), ROOT.Double()
#  abf.GetPoint(i,x,y0)
#  
#  print x,y0


print "end"
