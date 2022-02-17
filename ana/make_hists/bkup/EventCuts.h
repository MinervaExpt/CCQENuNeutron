#ifndef _EVENT_SEL_H
#define _EVENT_SEL_H

#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
using namespace CCQENU_ANA;
using namespace Neutron_ANA;

bool EventCutBase( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, map<int,int> &cuts_counter, bool isData = true, int multiplicity=1, int pass_signalFunc = -1, int pass_michel = -1, int pass_nblobs =-1, int pass_singleProton = -1, int pass_extraProtons = -1, int pass_hasNeutronCandidate = -1, int pass_neutronCut = -1 )//-1 not in use, 0 must be false, 1 must be true
{
  //Common event cuts:
  if( !cutter->passInteractionVertex( evt ) ) return false;   ++cuts_counter[0+multiplicity*10];
  if( !cutter->passDeadTimeCuts( evt ) ) return false;        ++cuts_counter[1+multiplicity*10];
  if( !cutter->passNuHelicityCut( evt ) ) return false;       ++cuts_counter[2+multiplicity*10];       
  if (multiplicity != 0 && (multiplicity != evt->multiplicity) ) return false;
                                                              ++cuts_counter[3+multiplicity*10];
  if( !cutter->passTrackAngleCut( evt ) ) return false;       ++cuts_counter[8+multiplicity*10];


  if ( pass_signalFunc!=-1 && 
      pass_signalFunc != cutter->passSignalFunction( evt, 0, 0 ) )  return false; ++cuts_counter[7+multiplicity*10];
  if ( pass_michel!=-1 && 
      pass_michel     != cutter->passImprovedMichelCut( evt ) )     return false; ++cuts_counter[5+multiplicity*10];
  if ( pass_nblobs!=-1 && 
      pass_nblobs     != cutter->passNBlobs( evt ) )                return false; ++cuts_counter[6+multiplicity*10];

  //Some functions about neutrons here
  if (multiplicity == 1) return true;

  if ( pass_singleProton != -1 &&
      pass_singleProton != cutter->passSingleProtonCut( evt, 0, 0 ) ) return false; 
  if ( pass_extraProtons != -1 &&
      pass_extraProtons != cutter->passExtraProtonsCut( evt, 0, 0 ) ) return false; ++cuts_counter[4+multiplicity*10];
  return true;
}

bool EventCut( CCQENuEvent* evt, CCQENuCuts *cutter, NeutronBlobCuts* ncutter, map<int,int> &cuts_counter, int multiplicity=1, string sample = "Signal", bool isData = true)
{
     int pass_signalFunc = -1,  
         pass_michel = -1,  
         pass_nblobs =-1,  
         pass_singleProton = -1,  
         pass_extraProtons = -1,  
         pass_hasNeutronCandidate = -1,  
         pass_neutronCut = -1;
    if (sample == "Signal")
    {
      pass_michel = 1;
      pass_nblobs = 1;
    } else if ( sample == "MichelSideBand" )
    {
      pass_michel = 0;
      pass_nblobs = 1;
    } else if ( sample == "BlobSideBand" )
    {
      pass_michel = 1;
      pass_nblobs = 0;
    } else if ( sample == "BlobSideBand" )
    {
      pass_michel = 0;
      pass_nblobs = 0;
    } else if (sample == "SideBand")
    {
      pass_signalFunc = 0;
    } else
    {
      cout<<"Wrong Sample Selection."<<endl;
      cout<<"Select: Signal, MichelSideBand, BlobSideBand, MicBlobSideBand"<<endl;
      return false;
    }

    if (isData)
    {
      pass_singleProton = 1;
      pass_extraProtons = 1;
      pass_signalFunc = 1;
    }

   return EventCutBase( evt, cutter, ncutter, cuts_counter, isData, multiplicity, pass_signalFunc , pass_michel , pass_nblobs , pass_singleProton, pass_extraProtons , pass_hasNeutronCandidate , pass_neutronCut );
}




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
#endif
