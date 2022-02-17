import ROOT
ROOT.TH1.SetDefaultSumw2()
def Chisq1D(hmc0, hdata0,xmin=0,xmax=-9):
  hmc = hmc0.Clone("chisqMC")
  hdata = hdata0.Clone("chisqData")
  for h in [hmc,hdata]:
    h.GetXaxis().SetRangeUser(xmin,xmax)

  Bmin,Bmax = 1, hdata.GetNbinsX()+1
  if (xmax < xmin):
    bmin, bmax = 1, hdata.GetNbinsX()+1
  else:
    bmin, bmax = hdata.FindBin(xmin), hdata.FindBin(xmax)
  ndof = bmax-bmin-1

  hdiff = hdata-hmc
  hdiverr = hdata.Clone("herr")
  hdiverr.Reset()

  for ibin in range(hdata.GetNbinsX() ):
    num_bin = ibin+1
    err = hdata.GetBinError( num_bin )
    diverr = 0 if err == 0 else 1/err
    hdiverr.SetBinContent(num_bin, diverr )
  hchisq = (hdiff*hdiverr)*(hdiff*hdiverr)
  return hchisq.Integral(), ndof


def Chisq1DVec(mc, data,covariance,xmin=0,xmax=-9):
  #covariance is TMatrixT<double>
  #print hmc0
  d = mc-data
  Cov=covariance
  svd = ROOT.TDecompSVD(Cov)
  svd.Invert(Cov)
  invCov=Cov
  res = d*invCov
  res = d*res

  return res

def Chisq1DGraph(mc, data,covariance,xmin=0,xmax=-9, useCov=True):
  #covariance is TMatrixT<double>
  #print hmc0
  vmc=ROOT.TVectorD( mc.GetN(), mc.GetY() )
  x = ROOT.TVectorD( mc.GetN(), mc.GetX() )
  vdata=ROOT.TVectorD( data.GetN(), data.GetY() )
  d = vmc-vdata

  ndf=mc.GetN()
  Cov=covariance.Clone("cov")
  if not useCov:
    for i in range(Cov.GetNrows() ):
      for j in range(Cov.GetNcols() ):
        if i!=j: Cov[i][j]=0


  svd = ROOT.TDecompSVD(Cov)
  svd.Invert(Cov)
  invCov=Cov

  if xmin <= xmax and xmin>=0:
    if xmax > mc.GetN()-1: xmax=mc.GetN()-1
    for i in range(0,xmin): d[i]=0
    for i in range(xmax+1,mc.GetN() ): d[i]=0
    for i in range(Cov.GetNrows()):
      for j in range(Cov.GetNcols() ):
        if i<xmin or i>xmax:
          invCov[i][j]=0
        if j<xmin or j>xmax:
          invCov[i][j]=0
    ndf=xmax-xmin+1
  #vdata.Print()
  #d.Print()
  #invCov.Print()

  res = d.Clone("res")
  res*=invCov
  #res.Print()

  #res.Print()
  #x.Print()
  return res*d, ndf





def Chisq1DHisto(hmc0, hdata0,covariance,xmin=0,xmax=-9):
  #covariance is TMatrixT<double>
  #print hmc0
  hmc = hmc0.Clone("chisqMC")
  hdata = hdata0.Clone("chisqData")
  #hu = hmc - hdata
  hu=hdata.Clone("hu")
  hu.Reset()
  for i in range( hdata.GetNbinsX() ):
    ibin = i+1
    x = hdata.GetBinCenter(ibin)
    dy = hmc.GetBinContent( hmc.GetXaxis().FindBin(x))-hdata.GetBinContent(ibin)
    hu.SetBinContent(ibin, dy )
  n_elements = hdata.GetNbinsX()

  highBin, lowBin = hdata0.GetXaxis().GetLast(),hdata0.GetXaxis().GetFirst()
  Nbins = highBin-lowBin+1

  #cov = ROOT.TMatrixD(Nbins,Nbins)
  #for i in range(Nbins):
  #  for j in range(Nbins):
  #    cov[i][j]=covariance[i+lowBin][j+lowBin]


  Imin, Imax = 0, n_elements
  IMIN,IMAX = Imin, Imax
  if (xmin < xmax):
    Imin=hdata.FindBin(xmin)-1
    Imax=hdata.FindBin(xmax)-1

  Cov = covariance*1e80
  #Cov *=1e80
  #print Cov
  svd = ROOT.TDecompSVD(Cov)
  svd.Invert(Cov)

  invCov = Cov*1e80
  #invCov *=0e80

  chisq = 0
  ndf = Imax-Imin
  for i in range(Imin, Imax):
    for j in range(Imin,Imax):
      dchisq= invCov[i+1][j+1]*hu.GetBinContent(i+1)*hu.GetBinContent(j+1)
      chisq+=dchisq
      print i+1, j+1
      print hu.GetBinCenter(i+1),hu.GetBinCenter(j+1)
      print invCov[i+1][j+1]
      print hu.GetBinContent(i+1),hu.GetBinContent(j+1)
      print dchisq, chisq
      print hu.GetBinError(i+1), hu.GetBinError(j+1)

  return chisq, ndf

  
  


def Chisq1DHisto2(hmc0, hdata0,covariance,xmin=0,xmax=-9):
  #covariance is TMatrixT<double>
  hmc = hmc0.Clone("chisqMC")
  hdata = hdata0.Clone("chisqData")
  hu = hmc - hdata
  n_elements = hdata.GetNbinsX()
  cov = []
  cov_arr = covariance.GetMatrixArray()
  for i in range(n_elements):
    row = []
    for j in range(n_elements):
      row.append(cov_arr.At(i*n_elements+j))
    cov.append(row)

  uovererr = [ ]
  for i in range(n_elements):
    diff = hu.GetBinContent(i+1)
    err = hu.GetBinError(i+1)
    v = diff/err if err>0 else 0
    uovererr.append(v)


  Imin, Imax= 0, n_elements
  IMIN, IMAX= 0, n_elements
  if (xmax > xmin ):
    Imin = hu.FindBin(xmin)-1
    Imax = hu.FindBin(xmax)-1
  chisq = 0
  ndf = n_elements - 1
  for i in range( Imin, Imax ):
    for j in range( IMIN, IMAX ):
      chisq+=uovererr[i]*uovererr[j]*cov[i][j]

  return chisq, ndf

#=============== USE THIS ==============
#calculate Chi2 given the covariance matraix
def ChisqInRangesCov( hmc0, hdata0,cov, ranges ):
  #ranges = [(xmin0, xmax0),(xmin1,xmax1)...(xminN,xmaxN)]
  ret = []
  for xmin, xmax in ranges:
    chisq,ndof = Chisq1DHisto(hmc0,hdata0,cov, xmin, xmax )
    ret.append((chisq,ndof))
  return ret


def ChisqInRanges( hmc0, hdata0, ranges ):
  #ranges = [(xmin0, xmax0),(xmin1,xmax1)...(xminN,xmaxN)]
  ret = []
  for xmin, xmax in ranges:
    chisq,ndof = Chisq1D(hmc0,hdata0, xmin, xmax )
    ret.append((chisq,ndof))
  return ret


#Additional Matrix Operations
# Symmetric Folding Operator
# Difference
# | 1      -1|  a1     a1-a5                                       
# |   1  -1  |  a2     a2-a4                                      
# |     0    |  a3 ===>a3-a3                                            
# |  -1   1  |  a4     a4-a2                                      
# |-1       1|  a5     a5-a1                                       
#
#    | 1       1|  a1     .5(a1+a5)                                  
#  1 |   1   1  |  a2     .5(a2+a4)                                  
#  - |     2    |  a3 ===>.5(a3+a3)                                      
#  2 |   1   1  |  a4     .5(a4+a2)                                  
#    | 1       1|  a5     .5(a5+a1)                                  
#
#
#
#

def GenerateFoldingAvgMatrix( nElements, hasOverflow = False ):# with overflow
  if hasOverflow: nElements+=2
  matrix = ROOT.TMatrixD(nElements, nElements)        #diagnals with 0.5
  for i in range(nElements):
    matrix[i][i] += 1
    matrix[i][nElements-1-i] += 1
  return matrix

def GenerateFoldingDiffMatrix( nElements, hasOverflow = False ):# with overflow
  if hasOverflow: nElements+=2
  matrix = ROOT.TMatrixD(nElements, nElements)        #diagnals with 0.5
  for i in range(nElements):
    matrix[i][i] += 1
    matrix[i][nElements-1-i] += -1
  return matrix

#def CovAsymA( h1 ):
#  #TH1D h1
#  #TMatrixD cov
#  nCols = h1.GetNbinsX()+2
#  A = GenerateFoldingDiffMatrix(nCols)*-1
#  for i in range(nCols):
#    vi = h1.GetBinContent(i)
#    for j in range(nCols):
#      vj = h1.GetBinContent(j)
#      s = (vi+vj)
#      A[i][j] =(2*vj/s/s)*A[i][j] if s != 0 else 0
#  return A

def CovAsymA( h1 ):  # Populate a folding covariance matrix
  #TH1D h1
  #TMatrixD cov
  nCols = h1.GetNbinsX()+2
  A = ROOT.TMatrixD(nCols, nCols)
  for i in range(nCols):
    vi = h1.GetBinContent(i)
    j = nCols-1-i
    vj = h1.GetBinContent(j)
    s = vi+vj
    A[i][i] = -2*vj/s/s if s!= 0 else 0 
    A[i][j] = 2*vi/s/s if s!= 0 else 0 
  return A


def VecAsymA( h1 ):
  #TH1D h1
  #TMatrixD cov
  nCols = h1.GetNbinsX()+2
  A = GenerateFoldingDiffMatrix(nCols)
  for i in range(nCols):
    vi = h1.GetBinContent(i)
    j = nCols - 1 - i
    vj = h1.GetBinContent(j)
    s = (vi+vj)
    A[i][i]+= (vi/s/s) if s != 0 else 0
    A[i][j]+= -(vj/s/s) if s != 0 else 0
  return A

def SumHist( h1, cov, xmin=0,xmax=9999 ):
  total, V = 0., 0.
  Range = range( cov.GetNcols() )
  for i in Range:
    xi = h1.GetBinCenter(i)
    if xi<xmin or xi>xmax: continue
    print h1[i], xi,"================"
    total+=h1[i]
    for j in Range:
      xj = h1.GetBinCenter(j)
      if xj<xmin or xj>xmax: continue
      print cov[i][j],xi,xj
      V+=cov[i][j]
  return (total, ROOT.TMath.Sqrt(V))

#=================== Matrix Transformation ================
def U(V,A):
  # A.V.A^T
  ret = ROOT.TMatrixD(A)
  ret*=V
  AT = ROOT.TMatrixD(A)
  AT.Transpose(A)
  ret*=AT
  return ret

#==================== USE THIS =============================
# Subtract left from the right, getting asymmetry histogram
# i.e. folding the distribution at origin
def GetAsymHist(h, covariance, foldingMatrix=None):
  
  ret = h.Clone( h.GetName()+"_fold" )
  ret.Reset()

  nElements = h.GetNbinsX()+2
  vec = ROOT.TVectorD( nElements )
  for i in range( nElements ): vec[i]=h[i]  # populate a vector of histogram content

  CovA = CovAsymA( h )            #covariance matrix 
  covFold = U( covariance, CovA ) #folded covariance, 
 
  foldingMatrix = VecAsymA(h)   # calculating folding matrix
  vec*=foldingMatrix            # fold data
  errVec = ROOT.TVectorD(vec)
  for i in range(errVec.GetNrows()): errVec[i] = covFold[i][i]**.5
  for i in range(nElements):
    ret[i]=vec[i]
    ret.SetBinError(i, errVec[i])
  return ret

def GetAsymCov(h, covariance):
  CovA = CovAsymA(h)
  return U( covariance, CovA )

def SumMatrix( matrix ):
  identity = ROOT.TMatrixD( matrix.GetNcols(),1 )
  for i in range( matrix.GetNcols() ): identity[i][0] = 1
  identity = identity.Transpose(identity)
  s = ROOT.TMatrixD( matrix )
  s*=identity
  s*=identity
  return s[0][0]

def GetFoldedHist(h, covariance, foldingMatrix):
  nBins = h.GetNbinsX()
  nElements = nBins+2
  
  ret = h.Clone( h.GetName()+"_fold" )
  ret.Reset()

  vec = ROOT.TVectorD( nElements )
  for i in range( nElements ): vec[i]=h[i]

  covFolded = U( covariance, foldingMatrix )

  vec*= foldingMatrix
  errVec = ROOT.TVectorD(vec)
  for i in range(errVec.GetNrows()): errVec[i] = covFolded[i][i]**.5
  for i in range(nElements):
    ret[i]=vec[i]
    ret.SetBinError(i, errVec[i])
  return ret


def GetFoldedHistChisq(hdata, cov, foldingMatrix, ranges, hmc=None, vectorFoldMatrix = None):
  if vectorFoldMatrix is None: vectorFoldMatrix = foldingMatrix
  nElements = cov.GetNcols()
  
  hMeasure = hdata.Clone( hdata.GetName()+"fold" )
  hMeasure.Reset()
  if hmc is None:
    hmc = hMeasure.Clone("hbase")

  vec = ROOT.TVectorD( nElements  )
  for i in range( nElements ): vec[i]=hdata[i]

  covFolded = U( cov, foldingMatrix )

  vec*= vectorFoldMatrix
  errVec = ROOT.TVectorD(vec)
  for i in range(errVec.GetNrows()): errVec[i] = covFolded[i][i]**.5
  for i in range(nElements):
    xi = hMeasure.GetBinCenter(i)
    if xi > 0:
      hMeasure[i]=vec[i]
      hMeasure.SetBinError(i, errVec[i])
    for j in range(nElements):
      xj = hMeasure.GetBinCenter(j)
      if xi>0 and xj>0: continue
      covFolded[i][j]=0
  

  return ChisqInRangesCov( hmc, hMeasure, covFolded, ranges )



def CalculateAvg(h, xmin=0, xmax=-1 ):
  IbinMin, IbinMax = 1, h.GetNbinsX()
  if (xmin< xmax):
    IbinMin, IbinMax = h.FindBin(xmin),h.FindBin(xmax)

  num, denum = 0, 0
  err = 0
  for ibin in range(IbinMin, IbinMax+1 ):
    x, v, e = h.GetBinCenter( ibin ), h.GetBinContent( ibin ), h.GetBinError(ibin)
    num+= x*v
    denum+= v
  mean= num/denum if denum != 0 else -99999
  return mean, sigma


def CalculateAvgCov(h, cov=None, xmin=0, xmax=-1 ):
  nrow = h.GetNbinsX()+2
  imin, imax = 0, h.GetNbinsX()+1
  if (xmin < xmax ):
    imin, imax = h.FindBin(xmin),h.FindBin(xmax)
  # calculate <x>
  if cov is None:
    Cov = ROOT.TMatrixD( nrow, nrow)
    for i in range(nrow):
      Cov[i][i]=(h.GetBinError(i)**2)
      if "gibuu" in h.GetName():
        Cov[i][i]=(h.GetBinError(i)*1e-20)**2
  else:
    Cov = cov

  num = 0
  variance = 0

  Positions = [h.GetBinCenter(i ) for i in range(nrow) ]
  Xsecs = [h.GetBinContent(i ) for i in range(nrow) ]
  den = sum( Xsecs[imin:imax+1] )

  weights = [ x/den for x in Positions ]

  for ibin in range( imin, imax+1):
    xs = Xsecs[ibin]
    wi = weights[ibin]
    num += xs*wi

    for jbin in range(imin,imax+1):
      wj = weights[jbin]
      covij = Cov[ibin][jbin]
      variance+= wi*wj*covij

  if den == 0:
    return -9999,-9999
  x_avg = num
  x_err = ROOT.TMath.Sqrt( variance )
  return x_avg, x_err

def PrintGraph( graph ):
  y =ROOT.TVectorD( graph.GetN(), graph.GetY() )
  x = ROOT.TVectorD( graph.GetN(), graph.GetX() )

  x.Print()
  y.Print()
