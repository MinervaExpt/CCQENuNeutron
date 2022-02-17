import os,sys
#Assumes SpecialSampleIncluded only

inoutdir = sys.argv[1]

cmd = "./SideBandFitHists %s/SideBandFitHists_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-BlobSideBand_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-MichelSideBand_SpecialSampleIncluded.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_SpecialSampleIncluded.root %s/SideBandFit_SpecialSampleIncluded.root 1 1"%(inoutdir,inoutdir,inoutdir,inoutdir,inoutdir,inoutdir,inoutdir,inoutdir)
os.system(cmd)
