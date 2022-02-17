import os,sys

#assumes SpecialSampleIncluded only

#needoutput, in 2 track, side band, migration , eff,4,1 1

inoutdirs = sys.argv[1:]


samples_to_run = [6,7,8,9,10,11,12,13,14]
samples= {}
samples[6]="dalphat"
samples[7]="dphit"
samples[8]="pn"
samples[9]="dpt"
samples[10]="dptx"
samples[11]="dpty"
samples[12]="signed"
samples[13]="signeddalphat"
samples[14]="signeddphit"


for mydir in inoutdirs:
    d = mydir.rstrip("\n")
    for s in samples_to_run:
        cmd = "./CrossSectionHists %s/CrossSection_CombinedPlaylists.root  %s/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_CombinedPlaylists.root %s/SideBandFit_CombinedPlaylists.root %s/Migration_MakeFlux-1_ApplyFlux-1_Multiplicity-2_%s_CombinedPlaylists.root %s/EffPurity_MakeFlux-1_CombinedPlaylists.root 4 1 1 %d"%(d,d,d,d,samples[s],d,s)
        print cmd
        os.system(cmd)
