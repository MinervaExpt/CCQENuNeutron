import sys,os

basedir = sys.argv[1]
extratag = sys.argv[2]
data = int(sys.argv[3])

#playlists = ["1A","1B","1C","1D","1E","1F","1G","1L","1M","1N","1O","1P","5A","6A","6B"]
playlists = ["5A","6A","6B","6C", "6D", "6E", "6F", "6G"]

for pl in playlists:
    print "I am working on %s"%(pl)
    print "Treating this as data=%d"%(bool(data))

    files = os.popen("ls %s/minervame%s%s/*.root"%(basedir,pl,extratag)).readlines()

    myoutname = "CCQENu_minervame%s_MC_Inextinguishable_merged.txt"%(pl)
    if bool(data): myoutname = myoutname.replace("MC","DATA")
    print myoutname
    myout = open(myoutname,"w")
    for f in files:
        myout.write(f.replace("/pnfs","root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/"))
        



