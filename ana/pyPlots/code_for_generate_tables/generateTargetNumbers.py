import PlotUtils
from PlotUtils import TargetUtils
from decimal import Decimal


target = TargetUtils()

elements = [ "H","C","N", "O","Al","Si","Cl","Ti" ]
Z = [ 1, 6,7, 8, 13, 14,17,22]

print "Element & TargetUtils & DDDB & DDDB/TargetUtils"
nplanes= 2*(80-27+1)
SUM=0
for i in range(len(elements)):
  el, z = elements[i], Z[i]
  data = target.GetTrackerElementNAtoms(z,nplanes,False)
  mc=target.GetTrackerElementNAtoms(z,nplanes,True)
  ratio = "NA"
  SUM+=data
  if data != 0:
    ratio = mc/data
  print "%s & %s & %s & %s"%(el, data, mc, ratio)

NucleonModes=["N atoms", "N protons", "N neutrons", "N nucleons"]

print "\\begin{tabular}{|c|cc|c|}\\hline"
print " & TargetUtils & DDDB & DDDB/TargetUtils \\\\\hline"
print "N atoms & %.4E & %.4E & %.4f \\\\\hline"%( target.GetTrackerNAtoms(nplanes, False), target.GetTrackerNAtoms(nplanes, True) , target.GetTrackerNAtoms(nplanes, True)/target.GetTrackerNAtoms(nplanes, False) )
print "N protons & %.4E & %.4E & %.4f \\\\\hline"%( target.GetTrackerNProtons(nplanes, False), target.GetTrackerNProtons(nplanes, True) , target.GetTrackerNProtons(nplanes, True)/target.GetTrackerNProtons(nplanes, False) )
print "N neutrons & %.4E & %.4E & %.4f\\\\\hline "%( target.GetTrackerNNeutrons(nplanes, False), target.GetTrackerNNeutrons(nplanes, True) , target.GetTrackerNNeutrons(nplanes, True)/target.GetTrackerNNeutrons(nplanes, False) )
print "N nucleons & %.4E & %.4E & %.4f\\\\\hline "%( target.GetTrackerNNucleons(nplanes, False), target.GetTrackerNNucleons(nplanes, True) , target.GetTrackerNNucleons(nplanes, True)/target.GetTrackerNNucleons(nplanes, False) )
print "\\end{tabular}"


for i in range(len(elements)):
  el, z = elements[i], Z[i]
  data = target.GetTrackerElementMassFraction(z,False)
  mc = target.GetTrackerElementMassFraction(z,True)
  a = target.GetTrackerElementA(z)+z
  if data != 0:
    ratio = mc/data
  print "%s & %s & %s & %s"%(el, a,data, mc)

