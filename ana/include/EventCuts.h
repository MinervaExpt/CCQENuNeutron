#ifndef _EVENT_CUT_H
#define _EVENT_CUT_H

#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
using namespace CCQENU_ANA;


//===========Muon P Cut===============
bool useMuonPCut = true;

double MuonPMin = 1.5;
double MuonPMax = 20;

//===========Modify Eb===============
bool modify_Eb = false;
int mode_eb = 1;

//===========Total Blob Energy Cut===============
bool useTotalECut = false;

//===========Explicit Enu Cut===============// Never applied
bool useEnuCut = false;

bool isDebug = false;

double MuonCylinderCut = 100;//mm
double NeutronZCut = 100;//mm
double NeutronDistCut = 150;//mm
double MuonConeAngleDeg = 15;//degree
bool cutForwardNeutron = true;

bool useRecoilCut = false;
string orderPar = "ClusTotalE";


bool dropClusters = false;
double discardFraction = 0.2;
double discardThreshold = 10;

int pass_neutronClusMaxE = -1;
int pass_neutronQsqCut = 1;// model this a part of the recoil cut

//===================================================
std::vector<bool> CommonNeutronCuts( CCQENuEvent* evt, NeutronBlobCuts *ncutter, string par="ClusTotalE" )
{
  std::vector<bool> ret;
  bool hasMainCandidate, muonConeCut, muonCylinderCut, forward, backward,outsideVtx, passClusMaxE;
  bool passQsq;
  hasMainCandidate = ncutter->GetMainCandidate( evt, par );//ClusMaxE parameter currently not implemented
  if (hasMainCandidate)
  {
    muonCylinderCut = ncutter->DistanceToMuonCut( evt, MuonCylinderCut );
    muonConeCut = ncutter->MuonConeCut( evt, ncutter->MainNeutronCandidate(), MuonConeAngleDeg, true, 0 );
    forward = ncutter->ForwardNeutronCut( evt, ncutter->MainNeutronCandidate(), cutForwardNeutron,NeutronZCut );
    backward = ncutter->ForwardNeutronCut( evt, ncutter->MainNeutronCandidate(), !cutForwardNeutron,NeutronZCut );
    outsideVtx = ncutter->DistanceToVtxCut(evt, NeutronDistCut, false );
    passClusMaxE = ncutter->MainNeutronCandidate()->TotalE > 10 ;
    passQsq = ncutter->CandidateQsqCut( evt, 0, 0 );
  }

  ret.push_back( hasMainCandidate );
  //ret.push_back( (muonCylinderCut && muonConeCut) );
  ret.push_back( (muonCylinderCut) );
  ret.push_back( forward );
  ret.push_back( backward );
  ret.push_back( outsideVtx );
  ret.push_back( passClusMaxE );
  ret.push_back( passQsq );
  return ret;
}


//===================================================


bool CutSelector( bool AcceptTrueOrFalse, bool TrueOrFalse )
{
  if (AcceptTrueOrFalse == 1 ) return TrueOrFalse;
  else return !TrueOrFalse;
}

bool NucleonKinematicsCut( XYZVector beam, XYZVector vec )
{
  double angle = ACos(beam.Unit().Dot(vec.Unit()) )*180/TMath::Pi();
  bool passAngle = angle < 70;
  double p = vec.R();
  bool passMomentum = ( p> 0.8 && p<1.2 );

  return (passAngle && passMomentum);
}


bool EventCutBase( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, map<int,int> &cuts_counter, bool isData = true, int multiplicity=1, int pass_signalFunc = -1, int pass_michel = -1, int pass_nblobs =-1, int pass_blobHitEnergy=-1, int pass_singleProton = -1, int pass_extraProtons = -1, int pass_neutronCut = -1, int pass_neutronMuonCut = -1, int pass_neutronVtxCut = -1, int pass_neutronClusMaxE = -1, string par = "ClusTotalE" )//-1 not in use, 0 must be false, 1 must be true
{
  if (isDebug) cout<<"enter eventcutbase"<<endl;
  if (multiplicity != 0 ) 
  { 
    if (multiplicity != evt->multiplicity)  return false;
  }
 
  //Common event cuts:
  if( !cutter->passInteractionVertex( evt ) ) return false;   ++cuts_counter[0+multiplicity*10];
  if( !cutter->passDeadTimeCuts( evt ) ) return false;        ++cuts_counter[1+multiplicity*10];
  if( !cutter->passNuHelicityCut( evt ) ) return false;       ++cuts_counter[2+multiplicity*10];       
  if( !cutter->passTrackAngleCut( evt ) ) return false;       ++cuts_counter[8+multiplicity*10];


  //differentiated cuts
  //if(isDebug)cout<<"2"<<endl;
  if ( pass_michel!=-1 ) 
  {
    if  ( !CutSelector(pass_michel, cutter->passImprovedMichelCut( evt )) )     return false; ++cuts_counter[5+multiplicity*10];
      if(isDebug)cout<<"pass_michel, "<<endl;
  }

  //if(isDebug)cout<<"3"<<endl;
  if ( pass_nblobs!=-1 )
  {
      if( !CutSelector(pass_nblobs , cutter->passNBlobs( evt ) ) )               return false; ++cuts_counter[6+multiplicity*10];
      //if( !CutSelector(pass_nblobs , cutter->passBlobTotalEnergyNHits( evt ) ) )               return false; ++cuts_counter[6+multiplicity*10];
      //if( !CutSelector(pass_nblobs , cutter->passNBlobs( evt ) ) )               return false; ++cuts_counter[6+multiplicity*10];

      if(isDebug)cout<<"pass_nblobs, "<<endl;
  }
  if ( pass_blobHitEnergy!=-1 )
  {
      if( !CutSelector(pass_nblobs , cutter->passBlobTotalEnergyNHits( evt ) ) ) return false;
  }

  //if(isDebug)cout<<"4"<<endl;
  //Some functions about neutrons here
  std::vector<bool> commonNeutronCuts = CommonNeutronCuts( evt, ncutter, par );
  if(isDebug)
  {
    //cout<<"5"<<endl;
    cout<<"got vector of commonNeutronCuts"<<endl;
    cout<<commonNeutronCuts[0]<<", "<<commonNeutronCuts[1]<<", "<<commonNeutronCuts[2]<<", "<<commonNeutronCuts[3]<<endl;
  }

  if ( pass_neutronCut !=-1 )
  {
    //cout<<"doing neutron?"<<endl;
      if (!CutSelector(pass_neutronCut, commonNeutronCuts[0]) )  return false; ++cuts_counter[9+multiplicity*10];
      //cout<<"pass_neutron, "<<endl;
  }
  if ( pass_neutronMuonCut!=-1)
  {
      if (!CutSelector(pass_neutronMuonCut, commonNeutronCuts[1]) )  return false; 
  }
  if ( pass_neutronVtxCut !=-1)
  {
    if(!CutSelector(pass_neutronVtxCut, commonNeutronCuts[4] )) return false;
  }
  if (pass_neutronClusMaxE != -1 )
  {
    if(!CutSelector( pass_neutronClusMaxE, commonNeutronCuts[5] ) ) return false;
  }

  if ( pass_signalFunc!=-1 ) 
  {
    if ( !CutSelector(pass_signalFunc, cutter->passSignalFunction( evt, 0, 0 )) )  return false; ++cuts_counter[7+multiplicity*10];
    if(isDebug)cout<<"pass_signalFunc, "<<endl;

    if (pass_neutronQsqCut != -1 )
    {
      //cout<<"do neutronQsqCut"<<endl;
      if(!CutSelector( pass_neutronQsqCut, commonNeutronCuts[6] ) ) return false;
      //cout<<"passed"<<endl;
    }

  }



  //Some functions about protons here
  if (multiplicity == 1) return true;

  if ( pass_singleProton != -1 )
  {
    if (pass_singleProton != cutter->passSingleProtonCut( evt, 0, 0 ) ) return false; 
    if (pass_singleProton != cutter->passHybridProtonNodeCut(evt, 10 ) ) return false;
  }

  if ( pass_extraProtons != -1 )
  {
    if ( !CutSelector(pass_extraProtons, cutter->passExtraProtonsCut( evt, 0, 0 )) ) return false; ++cuts_counter[4+multiplicity*10];
      //cout<<"pass_extraproton, "<<endl;
    if ( !CutSelector(pass_extraProtons, cutter->passExtraTracksProtonsCut( evt, NULL, 0 )) ) return false;;
      //cout<<"pass_extratrackproton, "<<endl;
  }


  return true;
}

bool EventCut( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, map<int,int> &cuts_counter, int multiplicity=1, string sample = "Signal", bool isData = true, int npMode = 1, bool useRecoilCut = true, string par="ClusTotalE")
{
  //np mode: 
  //0: neutron
  //1: proton
  //10: neutron proton
  int pass_signalFunc = -1,  
      pass_michel = -1,  
      pass_nblobs =-1,  
      pass_blobHitEnergy = -1,
      pass_singleProton = -1,  
      pass_extraProtons = -1,  
      pass_hasNeutronCandidate = -1,  
      pass_neutronCut = -1,
      pass_neutronMuonCut = -1,
      pass_neutronVtxCut = -1,
      pass_neutronClusMaxE = -1;


  if(useRecoilCut) pass_signalFunc = 1;
  else pass_signalFunc = -1;

 
  if(npMode == 1 || npMode == 10)
  {
    pass_singleProton = 1;
    pass_extraProtons = 1;
  }
  if(npMode == 0 || npMode ==10)
  {
    pass_neutronCut = 1;
    pass_neutronMuonCut = 1;
  }
  if( npMode == 10 )
  {
    pass_neutronVtxCut = 1;
  }

  int neutronOnlyMode = (npMode == 0);
  if(neutronOnlyMode)
  {
    //pass_neutronClusMaxE = 1;

  }
  if (sample == "Signal")
  {
    if( npMode != neutronOnlyMode )
    {
      pass_michel = 1;
      pass_nblobs = 1;
      pass_blobHitEnergy = 1;
    }
    else
    {
      pass_michel = 1;
      pass_nblobs = 1;
      pass_blobHitEnergy = 1;
    }
  } 
  else if ( sample == "MichelSideBand" )
  {
    if( npMode != neutronOnlyMode )
    {
      pass_michel = 0;
      pass_nblobs = 1;
      pass_blobHitEnergy = 1;
    }
    else
    {
      pass_michel = 0;
      pass_nblobs = 1;
      pass_blobHitEnergy = 1;
    }
  } 
  else if ( sample == "BlobSideBand" )
  {
    if( npMode != neutronOnlyMode )
    {
      pass_michel = 1;
      pass_nblobs = 0;
      pass_blobHitEnergy = 0;
    }
    else
    {
      pass_michel = 1;
      pass_nblobs = 0;
      pass_blobHitEnergy = 0;
    }
  } 
  else if ( sample == "MicBlobSideBand" )
  {
    if( npMode != neutronOnlyMode )
    {
      pass_michel = 0;
      pass_nblobs = 0;
      pass_blobHitEnergy = 0;
    }
    else
    {
      //pass_neutronVtxCut = 1;
      pass_michel = 0;
      pass_nblobs = 0;
      pass_blobHitEnergy = 0;
    }
  } else if (sample == "SignalTest")
  {
    if( npMode != neutronOnlyMode )
    {
      pass_michel = 1;
      pass_nblobs = 1;
    }
    else
    {
      //pass_neutronVtxCut = 1;
      pass_michel = 1;
      pass_nblobs = -1;
    }
    pass_signalFunc = -1;
  } 

  else if (sample == "SideBand")
  {
    pass_signalFunc = 0;
  } 
  else
  {
    cout<<"Wrong Sample Selection."<<endl;
    cout<<"Select: Signal, MichelSideBand, BlobSideBand, MicBlobSideBand"<<endl;
    return false;
  }
    return EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE, par );
}

bool EventCutWithoutNeutron( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, map<int,int> &cuts_counter, int multiplicity=1, string sample = "Signal", bool isData = true, int npMode = 1, bool useRecoilCut = true)
{
  //np mode: 
  //0: neutron
  //1: proton
  //10: neutron proton
  int pass_signalFunc = -1,  
      pass_michel = -1,  
      pass_nblobs =-1,  
      pass_blobHitEnergy = -1,
      pass_singleProton = -1,  
      pass_extraProtons = -1,  
      pass_hasNeutronCandidate = -1,  
      pass_neutronCut = -1,
      pass_neutronMuonCut = -1,
      pass_neutronVtxCut = -1;


  if(useRecoilCut) pass_signalFunc = 1;
  else pass_signalFunc = -1;

  pass_blobHitEnergy = 1;
 
  if (sample == "Signal")
  {
    pass_michel = 1;
    pass_nblobs = 1;
  } 
  else if ( sample == "MichelSideBand" )
  {
    pass_michel = 0;
    pass_nblobs = 1;
  } 
  else if ( sample == "BlobSideBand" )
  {
    pass_michel = 1;
    pass_nblobs = 0;
  } 
  else if ( sample == "MicBlobSideBand" )
  {
    pass_michel = 0;
    pass_nblobs = 0;
  } else if (sample == "SignalTest")
  {
    pass_michel = 1;
    pass_nblobs = -1;
    pass_signalFunc = -1;
  } 

  else if (sample == "SideBand")
  {
    pass_signalFunc = 0;
  } 
  else
  {
    cout<<"Wrong Sample Selection."<<endl;
    cout<<"Select: Signal, MichelSideBand, BlobSideBand, MicBlobSideBand"<<endl;
    return false;
  }
    return EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_blobHitEnergy, pass_singleProton, pass_extraProtons , pass_neutronCut, pass_neutronMuonCut, pass_neutronVtxCut, pass_neutronClusMaxE );
}

template< class T >
bool IgnoreEvent(T* mc, bool isData, bool print = false )
{
  if (isData) return false;
  if(mc->truth_genie_wgt_MFP_N[1]>50. || mc->truth_genie_wgt_MFP_N[2]>50. || mc->truth_genie_wgt_MFP_N[3]>50. || mc->truth_genie_wgt_MFP_N[4]>50. || mc->truth_genie_wgt_MFP_N[5]>50. || mc->truth_genie_wgt_MFP_N[6]>50. || mc->truth_genie_wgt_MFP_N[7]>50.) {
    if(print)
    {
      cout << "Event ignored - truth_genie_wgt_MFP_N[1] == " << mc->truth_genie_wgt_MFP_N[1]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[2] == " << mc->truth_genie_wgt_MFP_N[2]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[3] == " << mc->truth_genie_wgt_MFP_N[3]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[4] == " << mc->truth_genie_wgt_MFP_N[4]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[5] == " << mc->truth_genie_wgt_MFP_N[5]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[6] == " << mc->truth_genie_wgt_MFP_N[6]<< endl;
      cout << "Event ignored - truth_genie_wgt_MFP_N[7] == " << mc->truth_genie_wgt_MFP_N[7]<< endl;
    }
        return true;
  }
    if(mc->truth_genie_wgt_MFP_pi[1]>50. || mc->truth_genie_wgt_MFP_pi[2]>50. || mc->truth_genie_wgt_MFP_pi[3]>50. || mc->truth_genie_wgt_MFP_pi[4]>50. || mc->truth_genie_wgt_MFP_pi[5]>50. || mc->truth_genie_wgt_MFP_pi[6]>50. || mc->truth_genie_wgt_MFP_pi[7]>50.) {
      if(print)
      {
        cout << "Event ignored - truth_genie_wgt_MFP_pi[1] == " << mc->truth_genie_wgt_MFP_pi[1]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[2] == " << mc->truth_genie_wgt_MFP_pi[2]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[3] == " << mc->truth_genie_wgt_MFP_pi[3]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[4] == " << mc->truth_genie_wgt_MFP_pi[4]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[5] == " << mc->truth_genie_wgt_MFP_pi[5]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[6] == " << mc->truth_genie_wgt_MFP_pi[6]<< endl;
        cout << "Event ignored - truth_genie_wgt_MFP_pi[7] == " << mc->truth_genie_wgt_MFP_pi[7]<< endl;
      }
      return true;
    }
    return false;
}

template bool IgnoreEvent(CCQENuEvent* mc, bool isData, bool print );
template bool IgnoreEvent(CCQENuTruth* mc, bool isData, bool print );

//============================================

bool WriteHistos( CCQENuUtils* utils , TFile* f)
{
  f->cd();
  for( map<std::string, MnvH1D**>::iterator it  = utils->histos1D.begin(); it != utils->histos1D.end(); ++it )
  {
    MnvH1D**h = (MnvH1D**)(it->second);
    for ( int i = kData; i< nHistos; ++i ) {
      h[i]->Write();
    }
  }
  for( map<std::string, MnvH2D**>::iterator it  = utils->histos2D.begin(); it != utils->histos2D.end(); ++it )
  {
    MnvH2D**h = (MnvH2D**)(it->second);
    for ( int i = kData; i< nHistos; ++i ) h[i]->Write();
  }
  for( map<std::string, MnvH3D**>::iterator it  = utils->histos3D.begin(); it != utils->histos3D.end(); ++it )
  {
    MnvH3D**h = (MnvH3D**)(it->second);
    for ( int i = kData; i< nHistos; ++i ) h[i]->Write();
  }
  return true;
}

bool PassRecoilCuts( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts *ncutter )
{
  if ( !CutSelector(1, cutter->passSignalFunction( evt, 0, 0 )) )  return false; 
  return true;
}
bool PassRecoilBlobQsqCutsCV( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts *ncutter )
{
  if ( !CutSelector(1, cutter->passSignalFunction( evt, 0, 0 )) )  return false; 
  if(!CutSelector(1, ncutter->CandidateQsqCut(evt,0,0) ) ) return false;
  return true;
}
bool PassRecoilBlobQsqCutsLoose( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts *ncutter )
{
  double recoilShift = -100;//MeV
  double q2shift = 0.2;//GeV^2 ~ E = Q^2/2M ~ 100 MeV
  if ( !CutSelector(1, cutter->passSignalFunction( evt, 0, recoilShift )) )  return false; 
  if(!CutSelector(1, ncutter->CandidateQsqCut(evt,0,q2shift) ) ) return false;
  return true;
}


bool muonPCut( bool use, double muon_p )
{
  if(!use) return true;
  if( muon_p > 10 ) return false;
  return true;

}

#endif
