import os,sys
#Assumes SpecialSampleIncluded only

indirs = sys.argv
for ddd in indirs[1:]:
#This is the full combined sample
    d = ddd.rstrip("\n")
    cmd = "./SideBandFit %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-BlobSideBand_CombinedPlaylists.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MichelSideBand_CombinedPlaylists.root %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-MicBlobSideBand_CombinedPlaylists.root %s/SideBandFit_CombinedPlaylists.root"%(d,d,d,d,d)
    os.system(cmd)

