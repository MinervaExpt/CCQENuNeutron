import ROOT
gA=-1.267
Mn = 0.93956542052

XiK=[0,1/6., 1/3., 1/2., 2/3., 5/6., 1 ]
parameters=dict()
parameters["AEp"]=[1.,0.9927,0.9898,0.9975,0.9812,0.9340,1.]
parameters["AMp"]=[1.,1.0011,0.9992,0.9974,1.0010,1.0003,1.]
parameters["AEpd"]=[1.,0.9839,0.9632,0.9748,0.9136,0.5447,-0.2682]
parameters["AMpd"]=[1.,0.9916,0.9771,0.9801,1.0321,1.0429,0.5084]
parameters["AMn25"]=[1.,0.9958,0.9877,1.0193,1.0350,0.9164,0.7300]
parameters["AMn43"]=[1.,0.9958,0.9851,1.0187,1.0307,0.9080,0.9557]
parameters["AEn25"]=[1.,1.1011,1.1392,1.0203,1.1093,1.5429,0.9706]
parameters["AEn43"]=[1.,1.1019,1.1387,1.0234,1.1046,1.5395,1.2708]
parameters["AFA25D"]=[1.,0.9207,0.9795,1.0480,1.0516,1.2874,0.7707]



def AN( name, xi ):
  A=0
  par= parameters[name]
  nPar = len(par)
  for i in range(nPar):
    j=i+1
    Pj = par[i]
    for k in range(1,nPar+1):
      Pj*= 1 if (j==k) else (xi-XiK[k-1])/(XiK[j-1] - XiK[k-1])
    A+=Pj
  return A

def GD(Q2,Ma=1.015):
  return gA/(1+Q2/Ma/Ma)**2

def tau(Q2):
  return Q2/4/Mn/Mn

def Xi(Q2):
  if Q2==0: return 0
  return 2/ (1+(1+1/tau(Q2))**0.5)

def FA(Q2):
  xi = Xi(Q2)
  gd = GD(Q2)
  return AN("AFA25D",xi)*gd

def GetGraphErr( graph ):
  n = graph.GetN()
  ret= ROOT.TGraphErrors()
  for i in range(n):
    Q2=graph.GetX()[i]
    ret.SetPoint(i, Q2, FA(Q2) )
  
  return ret

