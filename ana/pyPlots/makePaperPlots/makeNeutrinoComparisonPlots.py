import ROOT, PlotUtils
print "Loaded PlotUtils"
ROOT.gROOT.SetBatch()

fpath ="/minerva/data/users/tejinc/CCQENu/CCQENuNeutron/rootfiles/minervameNu/cuts_proton/"

path_fit_sf = fpath+"CV2ElasticNuwroSF/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_1_390_qe_res_2p2h_snp_scp_3.root"
path_fit_norm = fpath+"CV2Elastic/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_1_390_qe_res_2p2h_snp_scp_3.root"

path_hist_sf = fpath+"CV2ElasticNuwroSF/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal-Angles_CombinedPlaylists.root"
path_hist_norm = fpath+"CV2Elastic/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal-Angles_CombinedPlaylists.root"

CATS=["qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_singlechargedpion", "qelikenot_singleneutralpion"]
catConvert = {
    "qelike_qe_h":"qelike_qe_h",
    "qelike_qe_oth":"qelike_qe_oth",
    "qelike_res":"qelike_res",
    "qelike_dis":"qelike_dis",
    "qelike_oth":"qelike_oth",
    "qelike_2p2h":"qelike_2p2h",
    "qelikenot_singlechargedpion":"qelikenot_scp",
    "qelikenot_singleneutralpion":"qelikenot_snp",
    "qelikenot_multipion":"qelikenot_mp",
    "qelikenot_nopions":"qelikenot_np"}

names=["qelike_qe_h", "qelike_qe_oth", "qelike_2p2h", "qelike_res", "qelike_dis","qelike_oth", "qelikenot_singlechargedpion", "qelikenot_singleneutralpion","qelikenot_multipion","qelikenot_nopions"]

def GetFits(fFit):
  fits = dict()
  #for cat in ["qelike_qe_oth", "qelike_res", "qelike_2p2h", "qelikenot_scp", "qelikenot_snp"]:
  fFit.ls("hs_weights_yvar_bgType_*")
  for cat in CATS:
    h=fFit.Get("hs_weights_yvar_bgType_%s"%catConvert[cat])
    fits[cat]=h
  return fits

def GetHists(fHist):
  retDic=dict()
  for region in [0,1,2,3,4,5,99]:
    retDic[region]= { name:fHist.Get("h_q2qe_region_%02d_%s"%(region,name)) for name in names }
    retDic[region]["data"]=fHist.Get("h_q2qe_region_%02d_data"%(region))
    retDic[region]["mc"]=fHist.Get("h_q2qe_region_%02d_mc"%(region))

  print retDic
  return retDic

def ApplyFit(hists, fits): #map of <region, category>
  for region in [0,1,2,3,4,5,99]:
    hists[region]["mc"].Reset()
    for name in names:
      if name in CATS: hists[region][name].Multiply(hists[region][name], fits[name] )
      hists[region]["mc"].Add( hists[region][name] )
  return

def GetHistos( h_path, fit_path ):
  fHist = ROOT.TFile(h_path,"read")
  fFit = ROOT.TFile(fit_path, "read")
  return [ GetHists(fHist), GetFits(fFit) ]

SF=GetHistos( path_hist_sf, path_fit_sf )
NORM=GetHistos( path_hist_norm, path_fit_norm )

ApplyFit( SF[0], SF[1] )
ApplyFit( NORM[0], NORM[1] )

def GenerateCanvasObject(name,title):
  c = ROOT.TCanvas(name,title,800,800)

  xlow1,ylow1,xup1,yup1 = 0.0,0.4,1,1
  xlow2,ylow2,xup2,yup2 = 0, 0.00,1,0.4
  p1 = ROOT.TPad("p1","p1",xlow1,ylow1,xup1,yup1)
  p1.SetBottomMargin(0)
  p1.SetGrid()
  p1.Draw()

  p2 = ROOT.TPad("p2","p2",xlow2,ylow2,xup2,yup2)
  p2.SetTopMargin(0)
  p2.SetBottomMargin(0.3)
  p2.SetGrid()
  p2.Draw()
  return (c,p1,p2)


def DrawHists( SF, NORM, region ):
  sf,nonsf = SF[region], NORM[region]
  colors=[ROOT.kRed, ROOT.kGray+2]

  c1,pad1,pad2 = GenerateCanvasObject("c1","c1");
  pad1.cd()

  data = sf["data"].Clone("data")
  mc = [ sf["mc"].Clone("CV2ElasticSF"), nonsf["mc"].Clone("CV2Elastic") ]
  mc[0].SetLineColor(colors[0])
  mc[1].SetLineColor(colors[1])

  data.GetXaxis().SetTitle("Q^{2}_{QE} (GeV^{2})")
  data.GetYaxis().SetTitle("N Events / Bin")

  data.GetXaxis().SetNdivisions(4)

  labelfont = data.GetYaxis().GetLabelFont()
  titlefont = data.GetYaxis().GetTitleFont()
  data.GetYaxis().SetTitleSize(20)
  data.GetYaxis().SetTitleFont(titlefont+1)
  data.GetYaxis().SetTitleOffset(2.5)
  data.GetYaxis().SetLabelFont(labelfont+1)
  data.GetYaxis().SetLabelSize(15)

  data.Draw()
  
  for h in mc:
    h.SetLineWidth(3)
    h.Draw("histlsame")



  pad2.cd()

  dratio1 = data.Clone("ratio1")
  dratio2 = data.Clone("ratio2")

  dratio1.Divide(data, mc[0])
  dratio2.Divide(data, mc[1])
  ratios = [ dratio1, dratio2 ]
  for ratio in ratios:
    ratio.SetMarkerStyle(0)
    #ratio.SetLineWidth(3)

    labelfont = ratio.GetYaxis().GetLabelFont()
    titlefont = ratio.GetYaxis().GetTitleFont()
    ratio.SetLineWidth(3)

    ratio.GetYaxis().SetNdivisions(204)
    ratio.GetYaxis().SetTitleSize(20)
    ratio.GetYaxis().SetTitleFont(titlefont+1)
    ratio.GetYaxis().SetTitleOffset(2.5)
    ratio.GetYaxis().SetLabelFont(labelfont+1)
    ratio.GetYaxis().SetLabelSize(15)

    ratio.GetXaxis().SetTitleSize(20)
    ratio.GetXaxis().SetTitleFont(titlefont+1)
    ratio.GetXaxis().SetTitleOffset(4.)
    ratio.GetXaxis().SetLabelFont(labelfont+1)
    ratio.GetXaxis().SetLabelSize(15)
    ratio.GetXaxis().SetNdivisions(4)



  dratio1.SetLineColor(colors[0])
  dratio2.SetLineColor(colors[1])
  dratio1.GetYaxis().SetRangeUser(0.2,1.8)
  dratio1.GetYaxis().SetTitle("Data/Post Fit MC")
  dratio1.Draw("histl")
  dratio2.Draw("histlsame")

  pad1.SetLogx()
  pad1.SetLogy()
  pad2.SetLogx()
  pad1.cd()

  leg = ROOT.TLegend(0.18,0.7,0.5,0.9)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.AddEntry( mc[0], mc[0].GetName(),"l" )
  leg.AddEntry( mc[1], mc[1].GetName(),"l" )
  leg.SetTextSize(0.04)
  leg.Draw()
  c1.Print("plots/neutrinoModelComp_%02d.pdf"%region)



DrawHists(SF[0],NORM[0],0)
DrawHists(SF[0],NORM[0],1)
DrawHists(SF[0],NORM[0],2)
#c1 = ROOT.TCanvas("c1","c1",800,400 )
#c1.Divide(2,1)
#c1.cd(1)
#SF[0][0]["data"].Draw()
#SF[0][0]["mc"].SetLineColor( ROOT.kGray )
#NORM[0][0]["mc"].SetLineColor( ROOT.kGray+2 )
#SF[0][0]["mc"].Draw("histlsame")
#NORM[0][0]["mc"].Draw("histlsame")



