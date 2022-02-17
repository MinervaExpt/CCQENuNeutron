import os
import sys

if len(sys.argv) != 5:
  print "plot_maker.py Proton/Neutron SignalType input_file.root base_output_folder_name"
  sys.exit()

script_name, mode, SignalType, input_file, base_output_folder_name = sys.argv

if mode not in ["Proton","Neutron", "ProtonBK"]:
  print "Mode must be Proton or Neutron"
  sys.exit()

output_string = "../plots/%sSelection%s/*png"%(mode,SignalType)
function = "./1DME%sSelectionPlots"%mode
if mode == "ProtonBK":
  function = "./1DMEProtonSelectionPlots_bksub"

print "Function is ", function
save_folder = "../plots/%s/"%base_output_folder_name
if not os.path.exists(save_folder):
  os.makedirs(save_folder)

#cmd = ./1DMEProtonSelectionPlot
cmdBase = "%s %s "%(function, input_file)
if mode == "ProtonBK":
  cmdBase = "%s %s %s"%(function, input_file,"../rootfiles/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_minervame1D.root")

Normalization = [0]
PlotStyles = ["QE","QELike","QELike_split","QELike_split_rw","QELike_split_fsi","QELike_split_fsi_rw","QE_PionInFS", "QELike_PionInFS", "QELike_split_PionInFS" ]

for style in PlotStyles:
  path = save_folder+style
  print path
  if not os.path.exists( path ):
    os.makedirs( path )
  for norm in Normalization:
    cmd =" %s %d %s "%(cmdBase, norm, style)
    print cmd
    os.popen( cmd )
  copyCmd = "cp %s %s"%(output_string, path)
  print copyCmd
  os.popen(copyCmd)



