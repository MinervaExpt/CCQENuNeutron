
gA = -1.2723
mpi=.1395702
tcut=9*mpi*mpi

A0=[1,2.30,-.6,-3.8,2.3]
A1=[-0.759, 2.30, -0.6, -3.8, 2.3, 2.16, -0.896, -1.58, 0.823]

Q2max=1
t0 = tcut*(1- (1+Q2max/tcut)**.5)

print tcut, t0

def z(q2):
  t1 = (tcut-q2)**.5
  t2 = (tcut-t0)**.5
  return (t1-t2)/(t1+t2)

def FA(Q2,A):
  q2 = -Q2
  ret = 0
  for i in range(len(A)):
    ret+=A[i]*z(q2)**i
  return ret


def GetA0(A):
  ret=0
  q2=0
  for i in range(1,len(A)):
    ret+=A[i]*z(q2)**i
  a0= gA-ret
  return a0

A0[0]=GetA0(A0)
print A0[0]


def FADipole(Q2,Ma):
  return gA/(1+Q2/Ma/Ma)**2

Q2List=[ i*0.1 for i in range(10)]

#A[0]=FADipole(t0,1.014)

for Q2 in Q2List:
  print Q2, -FA(Q2,A0), -FA(Q2,A1), -FADipole(Q2,1.014), (FA(Q2,A0) - FA(Q2,A1))/FA(Q2,A1)*100
 
