import ROOT
import PlotUtils
from GlobalIncludes import *

class CCQENuPlotUtils(object):
  def __init__(self):
    self.QELikeHistos = ( 
        kHistos.kQELike_QE_H, 
        kHistos.kQELike_QE_OTH, 
        kHistos.kQELike_2p2h, 
        kHistos.kQELike_RES, 
        kHistos.kQELike_DIS, 
        kHistos.kQELike_OTH, 
        )
    self.QELikeNotPionHistos = ( 
        kHistos.kQELikeNot_SingleChargedPion,
        kHistos.kQELikeNot_SingleNeutralPion,
        kHistos.kQELikeNot_MultiPion,
        kHistos.kQELikeNot_NoPions
        )

    pass
    
  def BookHistos(self, f, hname, rebinx=1,rebiny=1):
    histos = []
    for name in names:
      h = f.Get( hname+"_"+name )
      #if type(h) in [ type(PlotUtils.MnvH1D), type(ROOT.TH1D) ]:
      #  h.Rebin(rebinx)
      #else:
      #  h.Rebin2D(rebinx,rebiny)
      histos.append( h )
      print h
    return histos

  def Project2DHistos(self, histos, axis="x", amin=-1,amax=0):
    ret = []
    for h in histos:
      hname = h.GetName()+"_"+axis.lower()
      ret.append(
          h.ProjectionX( hname,amin,amax ) if axis.lower() == "x" else
          h.ProjectionY( hname,amin,amax )
          )
    return ret

  def FormatHistosStack( self, histos ):
    for i, color in enumerate(colors):
      histos[i].SetLineColor(ROOT.kBlack)
      histos[i].SetLineWidth( 1 )
      histos[i].SetFillColor(color)

    for ih in self.QELikeNotPionHistos:
      histos[ih].SetFillStyle(3001)

    return

  def FormatHistosLine( self, histos ):
    for i in range(len(histos)):
      if histos[i] is None: continue
      color = colors[i]
      histos[i].SetLineColor( color )
      histos[i].SetLineWidth( 4 )
      histos[i].SetFillStyle( 0 );


  

  def FormatCanvas( self, canvas ):

    return


  def DrawStacked1D(self, canvas, histos, data = None, xtitle ="",ytitle="", opt="hist",ymax=None):
    self.FormatHistosStack( histos )
    self.DUMP = list()
    OPT = opt

    hs = ROOT.THStack("hs","hs")
    legend=ROOT.TLegend(0.7,0.6,0.9,0.9)

    canvas.cd()
    for i in self.QELikeNotPionHistos[::-1]:
      hs.Add(histos[i])

    for i in self.QELikeHistos[::-1]:
      hs.Add(histos[i])
    if ymax is None: hs.SetMaximum( hs.GetMaximum()*1.5)
    else: hs.SetMaximum( ymax )
    hs.Draw(opt)

    if data is not None:
      data.Draw("pesame")
      legend.AddEntry(data,"Data")
    for i in self.QELikeHistos+self.QELikeNotPionHistos:
      legend.AddEntry( histos[i], names[i],"f" )
    legend.Draw()
    hs.GetXaxis().SetTitle(xtitle)
    hs.GetYaxis().SetTitle(ytitle)
    canvas.Update()


    #latex = ROOT.TLatex()
    #latex.DrawLatexNDC(0.5,0.1, xtitle )
    #latex.DrawLatexNDC(0.5,0.1, xtitle )

    self.DUMP.append(legend)
    self.DUMP.append(hs)

    return

  def DrawHistos1D(self, canvas, histos,names, xtitle = "", ytitle="", opt="histl"):
    self.DUMP=list()
    canvas.cd()
    legend=ROOT.TLegend(0.7,0.6,0.9,0.9)
    OPT=opt
    for i, h in enumerate(histos):
      h.GetXaxis().SetTitle(xtitle)
      h.GetYaxis().SetTitle(ytitle)
      h.Draw(OPT)
      OPT=opt+"same"
      legend.AddEntry( h, names[i], "l" )
    legend.Draw()
    self.DUMP.append(legend)
    return

  def DrawHistos1D(self, canvas, histos,names, xtitle = "", ytitle="", opt="histl", doStat=False):
    self.DUMP=list()

    canvas.cd()
    legend=ROOT.TLegend(0.7,0.6,0.9,0.9)
    OPT=opt
    ymax = max([ h.GetMaximum() for h in histos ] )
    ymin = min( min([ h.GetMinimum() for h in histos ] ), 0 )
    for i, h in enumerate(histos):
      h.GetXaxis().SetTitle(xtitle)
      h.GetYaxis().SetTitle(ytitle)
      h.GetYaxis().SetRangeUser(ymin,ymax*1.1)
      h.Draw(OPT)
      OPT=opt+"same"
      legend.AddEntry( h, names[i], "l" )
    legend.Draw()
    
    if doStat:
      ltx = ROOT.TLatex()
      ltx.SetTextSize(0.02)
      ymax = 0.58
      xmin = 0.7
      ltx.DrawLatexNDC(xmin,ymax,"histo \t mean \t rms")
      for i,h in enumerate(histos):
        ypos = ymax-(i+1)*0.05
        ltx.DrawLatexNDC(xmin,ypos,"%s \t %.3f \t %.2f"%(names[i],h.GetMean(), h.GetRMS() ) )
      self.DUMP.append(ltx)




    self.DUMP.append(legend)
    return


  DUMP = list()

  




