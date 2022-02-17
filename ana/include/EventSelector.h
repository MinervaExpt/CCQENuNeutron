#ifndef _EVENT_SEL_H
#define _EVENT_SEL_H
#include "include/CCQENuUtils.h"
#include "include/EventCuts.h"


bool PassReco( CCQENuEvent* evt, 
    NeutronBlobCuts* ncutter,
    KinematicsModifier* kineMod,
    PhysicsUtils* physicsUtils,
    GeoUtils* geoUtils,
    bool isMC, int npMode , vector<XYZTVector> &particles )
{
  bool ret = true;
  particles.clear();
  XYZTVector muon(evt->CCQENu_leptonE[0]/1000.,
                     evt->CCQENu_leptonE[1]/1000.,
                     evt->CCQENu_leptonE[2]/1000.,
                     evt->CCQENu_leptonE[3]/1000.);

  if( modify_Eb && isMC ) muon = kineMod->GetCorrectedMuon( evt, muon, mode_eb );
  
  double muon_p = muon.P();

  //CUT1
  if( useMuonPCut ) ret*=(muon_p>=MuonPMin && muon_p<=MuonPMax); // a muon momentum cut


  //================ Single Proton ======================
  //CUT2 for proton analysis
  bool hasProton = (evt->CCQENu_proton_Px_fromdEdx !=  -9999);
  if( npMode != 0 )
  {
    ret*=hasProton;
  }
  double Mp=0.93827231;
  XYZVector proton3D(evt->CCQENu_proton_Px_fromdEdx,evt->CCQENu_proton_Py_fromdEdx,evt->CCQENu_proton_Pz_fromdEdx);//MeV
  proton3D/=1000.;
  XYZTVector proton(proton3D.X(), proton3D.Y(), proton3D.Z(), TMath::Sqrt(Mp*Mp + proton3D.R()*proton3D.R() ) );

  if(isMC && modify_Eb && npMode != 0) proton = kineMod->GetCorrectedNucleon( evt,proton , mode_eb );

  //================ Neutron ======================
  std::vector<bool> commonNeutronCuts = CommonNeutronCuts(evt,ncutter, orderPar);
  bool hasNC = commonNeutronCuts[0] && commonNeutronCuts[1];
  bool isForward = commonNeutronCuts[2] && hasNC;
  bool isBackward = commonNeutronCuts[3] && hasNC;

  //CUT3 Neutron Cut
  if( npMode != 1 )  ret*=hasNC;

  XYZVector neutronDir(-1,-1,-1);
  if(hasNC)
  {
    neutronDir = ncutter->MainNeutronCandidate()->BlobVtxDir;
    neutronDir/=neutronDir.R();
  }
  double Mn = 0.9395654133;
  XYZTVector neutron(neutronDir.X(), neutronDir.Y(), neutronDir.Z(), TMath::Sqrt(Mn*Mn + neutronDir.R()*neutronDir.R() ) );

  //==================Neutrino====================
  //CUT4 Enu Cut

  XYZVector beam = geoUtils->BeamAxis( 0 ); 
  //if (useEnuCut)
  //{
  //}
  double binding_e = (npMode==0)? 0 : bindingE;
  int charge = (neutrinoMode)? -1:1;
  beam = beam.Unit()*physicsUtils->nuEnergyCCQE( beam, muon, charge, binding_e );
  XYZTVector neutrino(beam.X(), beam.Y(), beam.Z(), beam.R());
  

  particles.push_back( neutrino );
  particles.push_back( muon );
  particles.push_back( proton );
  particles.push_back( neutron );

  return ret;
}


template<class T> 
bool GetTrueParticles( T* evt, 
    CCQENuUtils* utils,
    KinematicsModifier* kineMod,
    PhysicsUtils* physicsUtils,
    GeoUtils* geoUtils,
    bool isMC, int npMode , vector<XYZTVector> &particles )
{
  bool ret = true;
  particles.clear();
  XYZTVector muon(   evt->mc_primFSLepton[0]/1000.,
                     evt->mc_primFSLepton[1]/1000.,
                     evt->mc_primFSLepton[2]/1000.,
                     evt->mc_primFSLepton[3]/1000.);

  if( modify_Eb && isMC ) muon = kineMod->GetCorrectedMuon( evt, muon, mode_eb );
  
  TLorentzVector proton_lorentz = utils->GetHighestTruePart4P( evt, 2212 );
  TLorentzVector neutron_lorentz = utils->GetHighestTruePart4P( evt, 2112 );
  XYZTVector proton = Convert( proton_lorentz );
  if( modify_Eb && isMC ) proton = kineMod->GetCorrectedNucleon( evt, proton, mode_eb );


  XYZTVector neutron = Convert( neutron_lorentz );
  

  XYZVector beam = geoUtils->BeamAxis( 0 ); 
  double binding_e = (npMode==0)? 0 : bindingE;
  int charge = (neutrinoMode)? -1:1;
  beam = beam.Unit()*physicsUtils->nuEnergyCCQE( beam, muon, charge, binding_e );
  XYZTVector neutrino(beam.X(), beam.Y(), beam.Z(), beam.R());

  bool hasProton = proton.E()>0;
  bool hasNeutron = proton.E()>0;
  
  if( modify_Eb && hasProton ) proton = kineMod->GetCorrectedNucleon( evt, proton, mode_eb );
  if( modify_Eb && hasNeutron ) neutron = kineMod->GetCorrectedNucleon( evt, neutron, mode_eb );

  particles.push_back( neutrino );
  particles.push_back( muon );
  particles.push_back( proton );
  particles.push_back( neutron );

  return ret;
}

#endif
