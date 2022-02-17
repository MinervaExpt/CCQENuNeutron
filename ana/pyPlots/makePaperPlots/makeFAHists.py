import PlotUtils
import ROOT
import array
import zexp
import BBBA2007
import fitter

ROOT.gROOT.SetBatch()

plotter=PlotUtils.MnvPlotter()

plotter.SetROOT6Palette(55)

#path to the matrix file
path = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/rootfiles/minervameAntiNu/"
plotpath="/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/panel_plotting/plots/minervameAntiNu/"

selver= "cuts_baseCut_fullStat2_SigFunc2/"
tune = "CV2ElasticNuwroSFn/"

nIt=4

fname_Fa_Rebinned=plotpath+selver+tune+"Signal_1_200/"+"Fa_Rebinned_it%d.root"%nIt
fname_Fa_Orig=plotpath+selver+tune+"Signal_1_200/"+"Fa_OriginalBinning_it%d.root"%nIt

print fname_Fa_Orig

fname_fhc_numu = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg14-minervame1D1M1NWeightedAve.root"
fname_fhc_numubar = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg-14-minervame1D1M1NWeightedAve.root"
fname_rhc_numubar = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg-14-minervame6A.root"
fname_rhc_numu = "/minerva/app/users/tejinc/cmtuser/Minerva_v22r1p1_NC_cvmfs/Ana/PlotUtils/data/flux/flux-gen2thin-pdg14-minervame6A.root"



f_fhc_numu = ROOT.TFile( fname_fhc_numu, "read")
f_fhc_numubar = ROOT.TFile( fname_fhc_numubar, "read")

f_rhc_numubar = ROOT.TFile( fname_rhc_numubar, "read")
f_rhc_numu = ROOT.TFile( fname_rhc_numu, "read")

h_fhc_numu = f_fhc_numu.flux_E_cvweighted
h_fhc_numubar = f_fhc_numubar.flux_E_cvweighted
h_rhc_numubar = f_rhc_numubar.flux_E_cvweighted
h_rhc_numu = f_rhc_numu.flux_E_cvweighted

q2map={ 0.003125 : 0.003125 ,
  0.009375 : 0.009375 ,
  0.01875  : 0.01875  ,
  0.03125  : 0.03125  ,
  0.04375  : 0.04375  ,
  0.075    : 0.075    ,
  0.125    : 0.125    ,
  0.175    : 0.175    ,
  0.25     : 0.25     ,
  0.35     : 0.3499   ,
  0.5      : 0.4995   ,
  0.7      : 0.6993   ,
  0.9      : 0.8991   ,
  1.1      : 1.099    ,
  1.6      : 1.582    ,
  3        : 2.815    ,
  5        : 4.768    ,
  8        : 7.596    }




def RebinCrossSection( h ):
  h_xs_rebinned = h.Rebin( n_newbins, "h_hydrogen_xs_rebinned", newbins)
  h_xs_rebinned.GetXaxis().SetTitle("Q^{2} (GeV^{2})")
  h_xs_rebinned.GetYaxis().SetTitle("d#sigma/dQ^{2} (cm^{2}/GeV^{2}/proton)")
  return h_xs_rebinned
h_q2=None
def PrintHisto( h, name , mode=0):
  with open("tables/%s"%name, "w") as f:
    if mode==0: f.write( "Q2,Xs,XsErr\n" )
    elif mode==1 : f.write( "Q2,Fa,FaErr\n" )
    else: f.write("E,Flux,FluxErr\n")
    for i in range(h.GetNbinsX()+2):
      Q2 = h.GetBinCenter(i)
      if mode in [0,1]: 
        if Q2>0 and Q2<10: Q2 = h_q2.GetBinContent( h_q2.FindBin( Q2 ) )
      Xs = h.GetBinContent(i)
      XsErr = h.GetBinError(i)
      f.write("%s,%s,%s\n"%(Q2,Xs, XsErr) )
  return


def FaDipole(Q2, Ma=1.015, gA=-1.2723):
    return gA/(1+Q2/(Ma*Ma))**2

def HistRatioToDipole( hist, h_q2, Ma=1.015, gA=-1.2723 ):
  ret = ROOT.TGraphErrors()
  for xbin in range(1, hist.GetNbinsX()+2):
    Q2 = h_q2.GetBinContent(xbin)
    Fa = hist.GetBinContent(xbin)
    FaErr= hist.GetBinError(xbin)
    dipole = FaDipole(Q2,Ma,gA)
    Fa/=dipole
    FaErr/=dipole
    ret.SetPoint(xbin-1,Q2,abs(Fa))
    ret.SetPointError(xbin-1,0,abs(FaErr))
  return ret

def HistToGraph( hist, h_q2, scale=-1 ):
  ret = ROOT.TGraphErrors()
  for xbin in range(1, hist.GetNbinsX()+2):
    Q2 = h_q2.GetBinContent(xbin)
    Fa = hist.GetBinContent(xbin)
    FaErr= hist.GetBinError(xbin)
    ret.SetPoint(xbin-1,Q2,Fa*scale)
    ret.SetPointError(xbin-1,0,abs(FaErr))
  return ret


transparency=0.2
def FormatGraph( h05,h07, meyer ):
  colors=[ROOT.kBlack, ROOT.kBlue, ROOT.kOrange+1]
  styles=[20,25,0]
  markersizes=[2,2,0]
  linewidths=[2,2,3]
  graphs=[h05,h07,meyer]

  for i in [0,1,2]:
    graphs[i].SetMarkerSize( markersizes[i] )
    graphs[i].SetMarkerStyle(styles[i])
    graphs[i].SetLineWidth( linewidths[i] )
    graphs[i].SetLineColor( colors[i] )
    graphs[i].SetMarkerColor( colors[i] )

  meyer.SetMarkerStyle(0)
  meyer.SetFillColorAlpha(ROOT.kOrange+1, 0.3 )
  meyer.SetFillStyle(1001)

def FormatGraph2( fit, color, style, tran=0.2 ):
  fit.SetMarkerStyle(0)
  fit.SetMarkerSize(0)
  fit.SetFillColorAlpha(color, tran )
  fit.SetLineColor(color)
  fit.SetFillStyle(style)
  fit.SetLineWidth(2);

def ScaleTGraphError( graph, scale ):
  ret = graph.Clone()
  X = graph.GetX() 
  Y = graph.GetY() 
  YErr = graph.GetEY()
  for i in range(graph.GetN() ):
    ret.SetPoint( i, X[i], Y[i]*scale )
    ret.SetPointError(i, 0, YErr[i]*scale)

  return ret
def ScaleTGraph( graph, scale ):
  ret = graph.Clone()
  X = graph.GetX() 
  Y = graph.GetY() 
  for i in range(graph.GetN() ):
    ret.SetPoint( i, X[i], Y[i]*scale )

  return ret


def GetMeyerFit(Q20=0.001, Q21=10, dQ2=0.005):
  ret = ROOT.TGraphErrors()
  Q2=Q20
  i=0
  scale=1
  while Q2<Q21:
    Fa=zexp.MeyerCV(Q2)
    FaErr=zexp.MeyerError(Q2)
    #FaErr=0

    ret.SetPoint(i,Q2,Fa)
    ret.SetPointError(i,0,FaErr)

    i+=1
    Q2+=dQ2
  return ret


def GetCVFit(Q20=0.001, Q21=10, dQ2=0.005,mode="tejin"):
  ret = ROOT.TGraphErrors()
  Q2=Q20
  i=0
  par, cov = zexp.meyerT0TCut
  t0,tcut=zexp.meyerT0TCut
  if mode.lower() == "tejin":
    par,cov = zexp.tejin_par2, zexp.tejin_cov2
    t0,tcut= zexp.tejin2T0TCut
  elif mode.lower() == "joint":
    par,cov = zexp.joint_par, zexp.joint_cov
  elif mode.lower() == "meyer":
    par,cov = zexp.meyer_par, zexp.meyer_cov
    t0,tcut=zexp.meyerT0TCut

  scale=1
  while Q2<Q21:
    if mode.lower() != "meyer":
      Fa=zexp.ZexpN2(Q2,par, t0,tcut,scale)
      FaErr=zexp.FitErrorN2(Q2,cov,t0,tcut, scale)
    else:
      Fa=zexp.MeyerCV(Q2)
      FaErr=zexp.MeyerError(Q2)

    ret.SetPoint(i,Q2,Fa)
    ret.SetPointError(i,0,FaErr)

    i+=1
    Q2+=dQ2
  return ret


def GetDipole(Q20=0.001, Q21=10, dQ2=0.005,Ma=1.015,gA=-1.2723):
  ret = ROOT.TGraphErrors()
  Q2=Q20
  i=0
  while Q2<Q21:
    Fa=FaDipole(Q2,Ma,gA)
    ret.SetPoint(i,Q2,Fa)
    i+=1
    Q2+=dQ2
  return ret



def GraphRatioToDipole( graph, Ma=1.015, gA=-1.2723 ):
  nPoints = graph.GetN()
  ret = ROOT.TGraphErrors()
  for i in range(nPoints):
    Q2=graph.GetX()[i]
    Fa=graph.GetY()[i]
    FaErr=graph.GetEY()[i]

    dipole = FaDipole(Q2,Ma,gA)
    Fa/=dipole
    FaErr/=dipole
    ret.SetPoint(i,Q2,Fa)
    ret.SetPointError(i,0,FaErr)
  return ret




def PrintCov( cov, name ):
  nrows = cov.GetNrows()
  ncols = cov.GetNcols()
  with open("tables/%s"%name, "w") as f:
    for i in range(nrows):
      for j in range(ncols):
        v = cov[i][j]
        f.write(str(v))
        if j!=ncols-1:
          f.write(",")
      f.write("\n")
  return


for i,fname in enumerate([fname_Fa_Orig ]):
#for i,fname in enumerate([fname_Fa_Rebinned]):

  name_app = "_Orig"
  if i == 1: name_app = "_Rebinned"

  f = ROOT.TFile(fname,"read")
  f.ls()
  h_xs = f.data_mnvh1d
  h_q2 = h_xs.Clone("h_q2")
  h_q2.Reset()
  for q20,q21 in q2map.items():
    h_q2.Fill(q20,q21)
  # cross section histogram
  h_xs.Draw() # already scaled

  

  # fa histograms
  fa_bba05 = f.Fa_BBA05
  fa_bbba07 = f.Fa_BBBA07

  fa_bba05G = f.Fa_BBA05G
  fa_bbba07G = f.Fa_BBBA07G

  fa_ratio = fa_bba05.Clone("fa_ratio")
  fa_ratio.Divide(fa_bba05, fa_bbba07)
  fa_ratio.GetYaxis().SetRangeUser(0.8,1.2)
  fa_ratio.Draw("histl")
  ROOT.c1.SetGrid()
  ROOT.c1.SetLogx(1)
  ROOT.c1.Print("plots/FA/fa_ratio%s.pdf"%name_app)
  # covariance
  cov_bba05 = f.Cov_BBA05
  cov_bbba07 = f.Cov_BBBA07


  jac_bba05 = f.Jacobian_BBA05
  jac_bbba07 = f.Jacobian_BBBA07

  jac_bba05.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.Print("plots/FA/jac_bba05%s.pdf"%name_app)

  jac_bbba07.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.Print("plots/FA/jac_bbba07%s.pdf"%name_app)
   

   
  replace=-1.6e38
  #for i in range(h_xs.GetNbinsX()+2):
  #  if abs(jac_bba05[i][i])>5e39: jac_bba05[i][i]=replace
  #  if abs(jac_bbba07[i][i])>5e39: jac_bbba07[i][i]=replace
  



  cov_data = h_xs.GetTotalErrorMatrix()
  corr_data = h_xs.GetTotalCorrelationMatrix()
  corr_sys_data = h_xs.GetTotalCorrelationMatrix(False,False)



  corr_1 = h_xs.GetSysCorrelationMatrix("GENIE_MaCCQEshape")
  corr_1.Draw("colz")
  ROOT.c1.Print("plots/FA/corr_GENIE_MaCCQEshape.pdf")


  h_cov_data = ROOT.TH2D( cov_data )
  h_cov_data.SetMinimum(1e-83)
  h_cov_data.Draw("colz")
  ROOT.c1.SetLogz()
  ROOT.c1.Print("plots/FA/cov_data%s.pdf"%name_app)

  h_corr_data = ROOT.TH2D( corr_data )
  h_corr_data.SetMinimum(-1)
  h_corr_data.SetMaximum(1)
  h_corr_data.Draw("colz")
  ROOT.c1.SetLogz(0)
  ROOT.c1.Print("plots/FA/corr_data%s.pdf"%name_app)

  h_corr_sys_data = ROOT.TH2D( corr_sys_data )
  h_corr_sys_data.SetMinimum(-1)
  h_corr_sys_data.SetMaximum(1)
  h_corr_sys_data.Draw("colz")
  ROOT.c1.SetLogz(0)
  ROOT.c1.Print("plots/FA/corr_sys_data%s.pdf"%name_app)



  h_cov_bba05 = ROOT.TH2D( cov_bba05 )
  h_cov_bba05.SetMinimum(0.0001)
  h_cov_bba05.SetMaximum(5)
  h_cov_bba05.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.SetLogz(1)
  ROOT.c1.Print("plots/FA/cov_bba05%s.pdf"%name_app)


  h_cov_bbba07 = ROOT.TH2D( cov_bbba07 )
  h_cov_bbba07.SetMinimum(0.0001)
  h_cov_bbba07.SetMaximum(5)
  h_cov_bbba07.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.SetLogz(1)
  ROOT.c1.Print("plots/FA/cov_bbba07%s.pdf"%name_app)


  cov_bba05 = jac_bba05*cov_data*jac_bba05
  h_cov_bba05 = ROOT.TH2D( cov_bba05 )
  h_cov_bba05.SetMinimum(0.0001)
  h_cov_bba05.SetMaximum(5)
  h_cov_bba05.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.SetLogz(1)
  ROOT.c1.Print("plots/FA/cov_bba05_redo%s.pdf"%name_app)


  cov_bbba07 = jac_bbba07*cov_data*jac_bbba07
  h_cov_bbba07 = ROOT.TH2D( cov_bbba07 )
  h_cov_bbba07.SetMinimum(0.0001)
  h_cov_bbba07.SetMaximum(5)
  h_cov_bbba07.Draw("colz")
  ROOT.c1.SetLogx(0)
  ROOT.c1.SetLogz(1)
  ROOT.c1.Print("plots/FA/cov_bbba07_redo%s.pdf"%name_app)

  #draw ratio of FA to Dipole

  meyerfit = GetMeyerFit()
  tejinfit = GetCVFit(mode="tejin")
  jointfit = GetCVFit(mode="joint")

  bodek_fit = f.Bodek_Fit


  #Draw raw distribution
  c3=ROOT.TCanvas("c3","c3")
  c3.cd()
  frame=ROOT.TH1D("frame","frame",500,0.01,7)
  frame.GetYaxis().SetRangeUser(0,2)
  frame.GetXaxis().SetTitle("Q^{2} (GeV^{2})")
  frame.GetYaxis().SetRangeUser(0,2)
  frame.GetYaxis().SetTitle("-F_{A}")
  frame.Draw()



  fa_bba05Pos = HistToGraph(fa_bba05,h_q2,-1)
  fa_bbba07Pos = HistToGraph(fa_bbba07,h_q2,-1)

  FormatGraph( fa_bba05Pos, fa_bbba07Pos, meyerfit )
  FormatGraph2( tejinfit, ROOT.kMagenta-1, 1001,.2 )
  FormatGraph2( jointfit, ROOT.kOrange+1,1001,.2 )

  meyerinv = ScaleTGraphError( meyerfit, -1 )
  tejininv = ScaleTGraphError( tejinfit, -1 )
  jointinv = ScaleTGraphError( jointfit, -1 )
  bodekinv = ScaleTGraph( bodek_fit, -1 )
  meyerinv.Draw("l3")
  tejininv.Draw("l3")
  #jointinv.Draw("l3")
  bodekinv.Draw("l")



  ##tejinfit.Draw("l3")
  ##jointfit.Draw("l3")


  ##bodek_fit.Draw("l")

  fa_bba05Pos.Draw("p")
  fa_bbba07Pos.Draw("p")

  leg=ROOT.TLegend(0.45,0.6,0.84,0.9)
  leg.SetBorderSize(0)
  #leg.SetFillStyle(0)

  leg.AddEntry( fa_bba05Pos, "Extracted F_{A} (BBBA05)", "pl")
  leg.AddEntry( fa_bbba07Pos, "Extracted F_{A} (BBBA07)", "pl")
  leg.AddEntry( meyerfit, "Fit: Z-exp, Meyer", "lf")
  leg.AddEntry( tejinfit, "Fit: Z-exp, Hydrogen", "lf")
  #leg.AddEntry( jointfit, "Fit: Z-exp, Joint", "lf")
  leg.AddEntry( bodek_fit,"Fit: BBBA07", "l")
  leg.Draw()
  c3.SetLogx()
  c3.SetLogy(0)
  c3.SetGrid()
  c3.Print("plots/FA/fa_dist_0_2%s.pdf"%name_app)




  #Draw Ratio
  ROOT.c1.cd()
  ratio_meyerfit=GraphRatioToDipole( meyerfit )
  ratio_tejinfit=GraphRatioToDipole( tejinfit )
  ratio_jointfit=GraphRatioToDipole( jointfit )


  ratio_bodekfit=GraphRatioToDipole( bodek_fit)
  ratio_bodekfit.SetLineColor(ROOT.kBlack)
  ratio_bodekfit.SetLineWidth(3)


  ratio_bbba07 = HistRatioToDipole( fa_bbba07, h_q2 )
  ratio_bbba05 = HistRatioToDipole( fa_bba05, h_q2 )

  FormatGraph( ratio_bbba05, ratio_bbba07, ratio_meyerfit )

  FormatGraph2( ratio_tejinfit, ROOT.kMagenta-1, 1001,.2 )
  FormatGraph2( ratio_jointfit, ROOT.kOrange+1,1001,.2 )


  frame=ROOT.TH1D("frame","frame",500,0,7)

  frame.GetYaxis().SetTitle("Ratio to Dipole (M_{A}=1.015)")

  leg=ROOT.TLegend(0.3,0.18,0.75,0.4)
  leg.SetBorderSize(1)
  #leg.SetFillStyle(0)
  frame.GetXaxis().SetTitle("Q^{2} (GeV^{2})")

  frame.GetYaxis().SetRangeUser(0,2)
  frame.Draw()
  ratio_meyerfit.Draw("l3")
  ratio_tejinfit.Draw("l3")
  #ratio_jointfit.Draw("l3")


  ratio_bodekfit.Draw("l")

  ratio_bbba05.Draw("p")
  ratio_bbba07.Draw("p")

  leg.AddEntry( ratio_bbba05, "Extracted F_{A} (BBBA05)", "pl")
  leg.AddEntry( ratio_bbba07, "Extracted F_{A} (BBBA07)", "pl")
  leg.AddEntry( ratio_meyerfit, "Fit: Z-exp, Meyer", "lf")
  leg.AddEntry( ratio_tejinfit, "Fit: Z-exp, Hydrogen", "lf")
  #leg.AddEntry( ratio_jointfit, "Fit: Z-exp, Joint", "lf")
  leg.AddEntry( ratio_bodekfit, "Fit: BBBA07", "l")
  leg.Draw()
  ROOT.c1.SetLogx()
  ROOT.c1.SetLogy(0)
  ROOT.c1.Print("plots/FA/fa_ratio_0_2%s.pdf"%name_app)


  c2=ROOT.TCanvas("c2","c2",400,200)
  c2.cd()

  frame.GetYaxis().SetTitle("")
  frame.GetXaxis().SetTitle("")

  frame.GetYaxis().SetRangeUser(0.1,10)
  frame.GetYaxis().SetNdivisions(404)

  frame.GetYaxis().SetLabelSize(0.1)
  frame.GetXaxis().SetLabelSize(0.1)

  frame.GetYaxis().SetTickSize(0.05)
  frame.GetXaxis().SetTickSize(0.05)


  frame.Draw();
  ratio_meyerfit.Draw("l3")

  ratio_tejinfit.Draw("l3")
  #ratio_jointfit.Draw("l3")

  ratio_bodekfit.Draw("l")

  ratio_bbba05.SetMarkerSize(1)
  ratio_bbba07.SetMarkerSize(1)

  ratio_bbba05.Draw("p")
  ratio_bbba07.Draw("p")

  c2.SetLogx()
  c2.SetLogy(0)
  c2.SetGrid()
  c2.Print("plots/FA/fa_ratio%s.pdf"%name_app)

  #Calculate Chi2:
  xmin,xmax=3,17
  useCov=True
  FA0=1.2723
  FA0=1
  chi2, ndf=fitter.Chisq1DGraph( BBBA2007.GetGraphErr( fa_bbba07G ), fa_bbba07G, cov_bbba07,xmin,xmax, useCov )
  print "BBBA2007 Chi2: ", chi2*FA0
  print "BBBA2007 NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( BBBA2007.GetGraphErr( fa_bbba07G ), fa_bba05G, cov_bba05,xmin,xmax, useCov )
  print "BBBA2005 Chi2: ", chi2*FA0
  print "BBBA2005 NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr2( fa_bbba07G, zexp.meyer_par ), fa_bbba07G, cov_bbba07,xmin,xmax, useCov )
  print "Zexp 07 Chi2: ", chi2*FA0
  print "Zexp NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr2( fa_bbba07G, zexp.meyer_par ), fa_bba05G, cov_bba05,xmin,xmax, useCov )
  print "Zexp 05 Chi2: ", chi2*FA0
  print "Zexp NDf: ", ndf


  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr3( fa_bbba07G, zexp.tejin_par ), fa_bbba07G, cov_bbba07,xmin,xmax, useCov )
  print "Zexp Tejin 07 Chi2: ", chi2*FA0
  print "Zexp Tejin NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr3( fa_bbba07G, zexp.tejin_par ), fa_bba05G, cov_bba05,xmin,xmax, useCov )
  print "Zexp Tejin 05 Chi2: ", chi2*FA0
  print "Zexp Tejin NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr2( fa_bbba07G, zexp.joint_par ), fa_bbba07G, cov_bbba07,xmin,xmax, useCov )
  print "Zexp Joint 07 Chi2: ", chi2*FA0
  print "Zexp Joint NDf: ", ndf

  chi2, ndf=fitter.Chisq1DGraph( zexp.GetGraphErr2( fa_bbba07G, zexp.joint_par ), fa_bba05G, cov_bba05,xmin,xmax, useCov )
  print "Zexp Joint 05 Chi2: ", chi2*FA0
  print "Zexp Joint NDf: ", ndf



  #fitter.PrintGraph( zexp.GetGraphErr( fa_bbba07G ) )





  PrintHisto( h_xs, "xsec"+name_app+".csv" )
  PrintCov( cov_data, "Cov_Xsec"+name_app+".csv")
  PrintCov( corr_data, "Corr_Xsec"+name_app+".csv")


  PrintHisto( fa_bba05, "Fa_BBA05"+name_app+".csv",1 )
  PrintHisto( fa_bbba07, "Fa_BBBA07"+name_app+".csv",1 )

  PrintCov( cov_bba05, "Cov_BBA05"+name_app+".csv" )
  PrintCov( cov_bbba07, "Cov_BBBA07"+name_app+".csv" )

  PrintCov( jac_bba05, "Jac_BBA05"+name_app+".csv" )
  PrintCov( jac_bbba07, "Jac_BBBA07"+name_app+".csv" )




  print "Printed Histo/Cov"
  f.Close()

print "numu"
PrintHisto(h_fhc_numu, "flux_fhc_numu.csv",2)
PrintHisto(h_rhc_numubar, "flux_rhc_numubar.csv",2)
PrintHisto(h_fhc_numubar, "flux_fhc_numubar.csv",2)
PrintHisto(h_rhc_numu, "flux_rhc_numu.csv",2)

