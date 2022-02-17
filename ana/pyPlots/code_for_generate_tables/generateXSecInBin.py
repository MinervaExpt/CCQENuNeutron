import ROOT

ROOT.gROOT.SetBatch()

playlists=["minervame5A"]
#folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat_SigFunc2/CV2ElasticNuwroSFn/"
folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_fullStat2_SigFunc2/CV2ElasticNuwroSFn/"
#folder = "/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameAntiNu/cuts_baseCut_chooseTotalE_cutBlobE_nuEnergy_FixSel_4_muonPCut_MHRW_FixEff_Lat_SigFunc2_NeutronUncertainty/CV2ElasticNuwroSFn/"
#_baseCut_fullStat2_SigFunc2

fname="CrossSectionHist_It4_w_lambda_1_200_CombinedPlaylists_hydrogen.root"

fpath = "%s/%s"%(folder, fname)
f = ROOT.TFile( fpath, "read")
f.ls()
h = f.h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section
h.Scale(1,"width")
h_totalErr=h.GetCVHistoWithError(True)
h_statErr=h.GetCVHistoWithStatError()


with open("xsec.csv",'w') as fout:

  for i in range(h.GetNbinsX()+2 ):
    ibin = i
    lowEdge = h.GetBinLowEdge(ibin)
    width = h.GetBinWidth(ibin)
    highEdge = lowEdge+width

    center = h.GetBinCenter(ibin)
    value = h.GetBinContent(ibin)
    err = h.GetBinError(ibin)
    total_err = h_totalErr.GetBinError(ibin)
    stat_err= h_statErr.GetBinError(ibin)



    #print "%d & %s & %s & %s & %.2f & %.2f & %.2f &%.2f\\\\\hline" % (ibin, lowEdge, highEdge, center, value, total_err, stat_err, ratio)

    save_values = map(str,(ibin, center, value, total_err, stat_err, err ))
    print save_values
    fout.write(", ".join( save_values ) + '\n')

  fout.close()

f.Close()

cov=ROOT.TH2D(h.GetTotalErrorMatrix())

for y in range( cov.GetNbinsY()+2 ):
  values=[]
  for x in range( cov.GetNbinsX()+2 ):
    values.append( cov.GetBinContent( x, y ) )
  print ",".join(map(str,values))


