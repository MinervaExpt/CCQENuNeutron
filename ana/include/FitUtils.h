#ifndef __Fit_UTILS_H
#define __Fit_UTILS_H
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"
#include "include/GeneralFunc.h"


vector<int> categories_to_fit({ kQELike_QE_OTH,  kQELike_RES, kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion} );
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_res", "qelike_2p2h","qelikenot_scp", "qelikenot_snp"});
vector<string> categories_to_fit_title_names({"QELike && QE OTH",  "QELike && RES", "QELike && 2p2h", "QELikeNot && Single Charged Pion", "QELikeNot && Single Neutral Pion"});

vector<int> histosUsed({ kData,  kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_2p2h, kQELike_RES, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } ); 

//Get and apply fits
vector<MnvH2D*> GetFitsHisto( TFile* f, string xvar )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH2D*) f->Get( Form( "hs_weights_%s_yvarbins_bgType_%s", xvar.c_str(), categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}
vector<MnvH1D*> GetFitsHisto( TFile* f  )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i<categories_to_fit_names.size();i++)
  {
    ret.push_back( (MnvH1D*) f->Get( Form( "hs_weights_yvar_bgType_%s", categories_to_fit_names[i].c_str() ) ) );
  }
  return ret;
}

template<class MnvHND>
void ApplyFit( MnvHND**h, vector<MnvHND*>&fits )
{
  for( int i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    h[cat]->Multiply( h[cat],fits[i] );
  }

  //h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
  //h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
 

  h[kMC]->Reset();
  h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  h[kQELike_QE]->Reset();
  h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}
template void ApplyFit( MnvH1D** h, vector<MnvH1D*>&fits );
template void ApplyFit( MnvH2D** h, vector<MnvH2D*>&fits );

void ApplyFit( MnvH2D**h, vector<MnvH1D*>&fits )
{
  for( int i = 0; i< categories_to_fit.size(); i++ )
  {
    int cat = categories_to_fit[i];
    cout<<"apply fit: "<<names[cat]<<endl;
    cout<<"weight: "<<fits[i]->GetName()<<endl;
    MnvH2D* fit = (MnvH2D*) h[kData]->Clone("fit");
    cout<<"fit nx: "<<fit->GetNbinsX()<<", ny: "<<fit->GetNbinsY()<<endl;
    fit->Reset();
    fit->ClearAllErrorBands();

    ExpandHisto( fits[i], fit, 1 );
    //for( int j = 0; j< fits[i]->GetNbinsX(); j++ )
    //{
    //  cout<<fits[i]->GetBinContent(j+1);
    //  cout<<", "<<fit->GetBinContent(1, j+1)<<endl;
    //}
    h[cat]->Multiply( h[cat],fit );
    delete fit;
  }

  //h[kQELikeNot_SingleNeutralPion]->Multiply(h[kQELikeNot_SingleNeutralPion], fits[2]);
  //h[kQELikeNot_SingleChargedPion]->Multiply(h[kQELikeNot_SingleChargedPion], fits[2]);
 

  h[kMC]->Reset();
  h[kQELike]->Reset();
  h[kQELikeNot]->Reset();


  h[kQELike_QE]->Reset();
  h[kQELike_QE]->Add(h[kQELike_QE_OTH]);
  h[kQELike_QE]->Add(h[kQELike_QE_H]);
 
  vector<int> qelike({ kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });
  for( auto cat: qelike ) h[kQELike]->Add( h[cat] );
  
  vector<int> qelikenot({ kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions } );
  for( auto cat: qelikenot ) h[kQELikeNot]->Add( h[cat] );

  h[kMC]->Add( h[kQELike] );
  h[kMC]->Add( h[kQELikeNot] );
}

template<class T>
void SubtractBackground( T** h, bool subqelike=false )
{
  h[kData]->AddMissingErrorBandsAndFillWithCV( *h[kMC] );
  h[kData]->Add( h[kQELikeNot], -1 );
  h[kMC]->Add( h[kQELikeNot], -1 );

  h[kQELikeNot]->Reset();
  h[kQELikeNot_SingleNeutralPion]->Reset();
  h[kQELikeNot_SingleChargedPion]->Reset();
  h[kQELikeNot_MultiPion]->Reset();
  h[kQELikeNot_NoPions]->Reset();

  if( subqelike )
  {
    h[kData]->Add( h[kQELike_QE_H] );
    h[kMC]->Add( h[kQELike_QE_H] );
    h[kData]->Add( h[kQELike], -1 );
    h[kMC]->Add( h[kQELike], -1 );

    h[kQELike_QE_OTH]->Reset();
    h[kQELike_2p2h]->Reset();
    h[kQELike_DIS]->Reset();
    h[kQELike_RES]->Reset();

    h[kQELike] = h[kQELike_QE_H];

  }
}

template void SubtractBackground(MnvH1D** h, bool subqelike );
template void SubtractBackground(MnvH2D** h, bool subqelike );

template<class T>
void RatioToMC( T** h )
{
  T* hmc = (T*) h[kMC]->Clone("hmc");
  vector<int> categories{ kData, kMC, kQELike, kQELike_QE_H, kQELike_QE_OTH, kQELike_QE, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH, kQELikeNot, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions };

  for( auto cat : categories ) h[cat]->Divide(h[cat], hmc);
}





#endif
