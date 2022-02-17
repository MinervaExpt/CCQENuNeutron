#include "include/GlobalIncludes.h"
#include "include/CCQENuUtils.h"
#include "include/CCQENuPlotUtils.h"

#include "include/GeneralFunc.h"
#include "include/CommonBins.h"


#include "PlotUtils/HyperDimLinearizer.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFitter.h"
using namespace CCQENU_ANA;
using namespace std;

bool doErrorBands = false;
bool normalizeData = false;

bool forceEqual = false;

double highweight = 5;
double lambda = 1;
//Global histograms because TMinuit needs them global...



//string fit_xvar = "h_ptmu";
//const int n_fit_bins = 9;
//int fit_bins_start[n_fit_bins] = {1,2,3,4,5,6,7,8,14};
//int fit_bins_end[n_fit_bins] = {1,2,3,4,5,6,7,13,14};

string fit_xvar = "h_q2qe";
vector<string> fit_xvars({ 
    "h_q2qe_region_00", 
    "h_q2qe_region_01", 
    "h_q2qe_region_02", 
    "h_q2qe_region_03", 
    "h_q2qe_region_04", 
    "h_q2qe_region_05",
    "h_q2qe_region_99",
    "h_q2qe_vtx_region_99",
    });
const int nSamplesPerSideBand = fit_xvars.size();
const int NumberOfSideBands = 4;
//1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18
//
//
//const int n_fit_bins = 1;
//int fit_bins_start[n_fit_bins] = {1};
//int fit_bins_end[n_fit_bins] =   {18};

//const int n_fit_bins = 2;
//int fit_bins_start[n_fit_bins] = {1,9};
//int fit_bins_end[n_fit_bins] =   {8,18};
//
//const int n_fit_bins = 3;
//int fit_bins_start[n_fit_bins] = {1,7,14};
//int fit_bins_end[n_fit_bins] =   {6,13,18};

//const int n_fit_bins = 5;
//int fit_bins_start[n_fit_bins] = {1,5,7,9,11};
//int fit_bins_end[n_fit_bins] =   {4,6,8,10,18};

//const int n_fit_bins = 15; //include overflow bins
//int fit_bins_start[n_fit_bins] = {0,1,6,7,8,9,10,11,12,13,14,15,17,19};
//int fit_bins_end[n_fit_bins] =   {0,5,6,7,8,9,10,11,12,13,14,16,18,19};

const int n_fit_bins = 16; //include overflow bins
int fit_bins_start[n_fit_bins] = {0,1,3,5,6,7,8,9,10,11,12,13,14,15,17,19};
int fit_bins_end[n_fit_bins] =   {0,2,4,5,6,7,8,9,10,11,12,13,14,16,18,19};


//const int n_fit_bins = 11;
//int fit_bins_start[n_fit_bins] = {1,6,7,8,9,10,11,12,13,14,15};
//int fit_bins_end[n_fit_bins] =   {5,6,7,8,9,10,11,12,13,14,18};

vector<int> Samples({ //0,1,2,3,4,5,6,7,8,
                      //9,10,11,12,13,14,15,16,17,
                      //18,19,20,21,22,23,24,25,26,
                      //27,28,29,30,31,32,33,34,35
                      //});
                      0, 1, 2, 3, 4, 5, 6, 7,
                      8, 9, 10, 11, 12, 13, 14, 15,
                      16, 17, 18, 19, 20, 21, 22, 23,
                      24, 25, 26, 27, 28, 29, 30, 31


                      //0,1,2,3,4,5,6,
                      //7, 8, 9, 10, 11, 12, 13,
                      //14, 15, 16, 17, 18, 19, 20,
                      //21, 22, 23, 24, 25, 26, 27});

                      //6,7,8,9,10,11, 
                      //12,13,14,15,16,17, 
                      //18,19,20,21,22,23});
                      //24,25,26,27});
                      });
//vector<int> Sidebands({ 0,1,2,3,4,5, 
//                        6,7,8,9,10,11, 
//                        12,13,14,15,16,17, 
//                        18,19,20,21,22,23});
//                        //24,25,26,27});
vector<int> Sidebands({1,5, 6,
                       14});
                       //13,14,15,16,17, 
                       //19,20,21,22,23});
const int nInput = Samples.size();
vector<double> RegionWeights(nInput,1);

//categories_to_fit: histos with scaling, types_to_stay_const, no scaling, needed to make kMC = sum(fit+const)
//vector<int> categories_to_fit({ kQELike_QE_OTH,kQELike_RES,kQELike_2p2h, kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion    });
//vector<string> categories_to_fit_names({"qe_oth", "res","2p2h","qelikenot_scp","qelikenot_snp","qelikenot_mp"});
//vector<int> categories_stay_const({ kQELike_QE_H, kQELike_DIS,kQELike_OTH, kQELikeNot_NoPions  });

vector<int> categories_to_fit({ kQELike_QE_OTH,kQELike_2p2h, kQELikeNot_SingleNeutralPion , kQELikeNot_SingleChargedPion  });
vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_snp", "qelikenot_scp"});
//vector<int> categories_to_fit({ kQELike_QE_OTH,kQELike_2p2h, kQELikeNot_SinglePion  });
//vector<string> categories_to_fit_names({"qelike_qe_oth", "qelike_2p2h","qelikenot_sp"});

vector<int> categories_stay_const({ kQELike_QE_H, kQELike_RES, kQELike_DIS,kQELike_OTH, kQELikeNot_MultiPion,kQELikeNot_NoPions  });


//only the background components included here. These histos makes up kQELikeNot
vector<int> bg_categories_to_scale = categories_to_fit;
vector<string> bg_categories_to_scale_names=categories_to_fit_names;
vector<int> bg_categories_stay_const= categories_stay_const;

//vector<int> all_universe_histos({kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH,kQELikeNot_SingleChargedPion, kQELikeNot_SingleNeutralPion, kQELikeNot_MultiPion, kQELikeNot_NoPions    });
vector<int> all_universe_histos({kQELike_QE_H, kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH,kQELikeNot_SinglePion, kQELikeNot_MultiPion, kQELikeNot_NoPions    });

const int FitCateogry=kMC;//kQELikeNot

//vector<TH1D*> input_histos(24, new TH1D[nHistos]);
vector<vector<TH1D>> input_histos(nInput, vector<TH1D>(nHistos) );

//--------------------------------------------
// all functions
//-------------------------------------------
double chi2Function(double *scale );
void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg);
void setScaleAndError(TH1D*h,vector<double> scales, vector<double> err);


void setScaleAndError(TH1D*h,vector<double> scales, vector<double> err);
void setScaleAndError(TH1D&h,vector<double> scales, vector<double> err);
void setScaleAndError(TH2D*h,vector<double> scales, vector<double> err);
void setScaleAndError(TH2D&h,vector<double> scales, vector<double> err);

void setScaleAndError( TH1D *h, TH1D* scale );
void setScaleAndError( TH2D *h, TH1D* scale );



void getProjection(MnvH2D**h, MnvH1D**proj,bool projX);
vector<TH1D*> GetSingleTypeFitScalesAndErrors();
vector<TH1D*> GetBackgroundScales(vector<TH1D*> &scale_histos, int background_type_id=kMC);


void initializeVector( vector<MnvH1D**> &h);
void initializeVector( vector<MnvH2D**> &h);


int SideBandFit( vector<string> filenames );
//scale = int[nFitHisto*nBins]
//==============================================
// Define chi2 function
//==============================================


void PrintMnvH1D( MnvH1D* h )
{
  cout<<"Print ======== "<<h->GetName()<<"=========="<<endl;
  cout<<"cv: "<<h->Integral();

  vector<string> vertNames = h->GetVertErrorBandNames();
  for( auto name : vertNames )
  {
    cout<<"VertErr: "<<name<<", "<<h->GetVertErrorBand( name )->Integral()<<endl;
    unsigned int nUniverses = h->GetVertErrorBand( name )->GetNHists();
    for( unsigned int iUniv = 0 ; iUniv<nUniverses; iUniv++ )
      cout<<"Univ "<<iUniv<<": "<<h->GetVertErrorBand( name )->GetHist(iUniv)->Integral()<<endl;
  }
  vector<string> latNames = h->GetLatErrorBandNames();
  for( auto name : latNames )
  {
    cout<<"LatErr: "<<name<<", "<<h->GetLatErrorBand( name )->Integral()<<endl;
    unsigned int nUniverses = h->GetLatErrorBand( name )->GetNHists();
    for( unsigned int iUniv = 0 ; iUniv<nUniverses; iUniv++ )
      cout<<"Univ "<<iUniv<<": "<<h->GetLatErrorBand( name )->GetHist(iUniv)->Integral()<<endl;
  }
}

void DrawNHists( vector<TH1D*> &hists, string basename )
{
  string folder = "/minerva/app/users/tejinc/cmtuser/Minerva_v21r1p1_NC_cvmfs/Ana/CCQENuNeutron/ana/background_tuning/test/";
  TCanvas *c = new TCanvas("c","c");
  hists[kData]->Draw("pe");
  hists[kMC]->Draw("histsame");
  hists[kQELike]->Draw("histlsame");
  hists[kQELikeNot]->Draw("histcsame");
  c->Print( Form( "%s/%s.png", folder.c_str(), basename.c_str() ) );
  delete c;
}

double chi2Function(double *scale )
{
  double chi2 = 0;
  double chi2_shape = 0;
  int nSidebands = Sidebands.size();
  int nFitCategories = categories_to_fit.size();
  double mc = 0;
  double data = 0;
  double chi2weight = 1;
  //cout<<"data chi2 = "<<endl;
  for( int i = 0; i< n_fit_bins; i++ )
  {
    chi2weight = 1;
    for( vector<int>::iterator itSideband = Sidebands.begin(); itSideband!=Sidebands.end();++itSideband )
    {
      //cout<<*itSideband<<endl;
      //index of sideband sample
      if( *itSideband == 16 ) chi2weight = highweight;
      //if( *itSideband == 1 ) chi2weight = highweight;
      else chi2weight = 1;
      int iSideband = itSideband - Sidebands.begin();
      chi2weight*=RegionWeights[iSideband];
      double sum = 0;
      //loop through the fitting histos, with scaling
      for( vector<int>::iterator itCat = categories_to_fit.begin(); itCat!=categories_to_fit.end();++itCat )
      {
        //cout<<"Histo Type: "<<names[*itCat]<<endl;
        int iCat = itCat - categories_to_fit.begin();
        double v_scale = scale[i*nFitCategories+iCat];
        if( iCat == 3 && forceEqual) scale[i*nFitCategories+iCat] = scale[i*nFitCategories+iCat-1];
        double v_sideband = input_histos[*itSideband][*itCat].Integral(fit_bins_start[i], fit_bins_end[i]);
        
        sum+= v_scale*v_sideband;
        mc+= v_scale*v_sideband;
        //cout<<"v: "<<v_sideband<<", "
        //  <<"Scale: "<<v_scale<<endl;
      }

      //loop through the const histos, without scaling
      for( vector<int>::iterator itCat = categories_stay_const.begin(); itCat!=categories_stay_const.end();++itCat )
      {
        //cout<<"Histo Type: "<<names[*itCat]<<endl;
        double v_sideband = input_histos[*itSideband][*itCat].Integral(fit_bins_start[i], fit_bins_end[i]);
        sum+=v_sideband; 
        mc+=v_sideband;
        //cout<<"v: "<<v_sideband<<", "
        //  <<"Scale: "<<"1"<<endl;
      }

      //return 0;
      double v_data = input_histos[*itSideband][kData].Integral(fit_bins_start[i], fit_bins_end[i]);
      data+=v_data;


      double dchi2 = (v_data>0)? pow(sum-v_data, 2)/v_data*chi2weight:0;
      chi2+= dchi2;
      //cout<< i<<", "<<iSideband<<" "<<sum<<", "<<v_data<<"--chi2: "<<dchi2<<"...."<<endl;

      ////chi2_shape
      //for( int ibin = fit_bins_start[i]; ibin <= fit_bins_end[i]; ibin++ )
      //{
      //  double sum_reg = 0;
      //  for( vector<int>::iterator itCat = categories_to_fit.begin(); itCat!=categories_to_fit.end();++itCat )
      //    sum_reg+= scale[i*nFitCategories+(itCat-categories_to_fit.begin())]*input_histos[*itSideband][*itCat].Integral(ibin,ibin);
      //    //sum_reg+= pow(scale[i*nFitCategories+(itCat-categories_to_fit.begin())], 2);
      //  for( vector<int>::iterator itCat = categories_stay_const.begin(); itCat!=categories_stay_const.end();++itCat )
      //    sum_reg+= input_histos[*itSideband][*itCat].Integral(ibin,ibin);
      //  double data_reg = input_histos[*itSideband][kData].Integral(ibin,ibin);
      //  chi2_shape+=pow(sum_reg-data_reg,2)*chi2weight;
      //}
    }
  }
  
  return chi2;//+chi2_shape;
}

double regularizer(double *scale )
{
  int nFitCategories = categories_to_fit.size();
  double reg = 0;
  //for( vector<int>::iterator itCat = categories_to_fit.begin(); itCat!=categories_to_fit.end();++itCat )
  for( int iCat = 0; iCat < nFitCategories; ++iCat )
  {

    for( int i = 0; i<n_fit_bins-2; i++ )
    {

      double s0 = scale[i*nFitCategories+iCat];
      double s1 = scale[(i+1)*nFitCategories+iCat];
      double s2 = scale[(i+2)*nFitCategories+iCat];
      reg+=pow( s0+s2-2*s1, 2);
    }
  }
  return reg;
}


double result = 0;
void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg)
{
  result = chi2Function(par)+lambda*regularizer(par);
  //cout<<result<<endl;
}


//==============================================
// Define setScaleAndErr
//==============================================
void setScaleAndError(TH1D*h,vector<double> scales, vector<double> err)
{
  for(int i=0;i<n_fit_bins;i++)
  {
    int bin_begin = fit_bins_start[i], bin_end = fit_bins_end[i];
    for( int ibin = bin_begin; ibin <=bin_end; ++ibin )
    {
      h->SetBinContent(ibin,scales[i]);
      h->SetBinError( ibin,err[i]);
    }
  }
}


void setScaleAndError(TH2D*h,vector<double> scales, vector<double> err)
{
  //fitbins is here in case you want a different application (but the code currently sets this as false
  //cout << "Fit bin = " << fitBins << endl;
  for( int i = 0; i< n_fit_bins; i++ )
  {
    int bin_begin = fit_bins_start[i], bin_end = fit_bins_end[i];
    for( int xbin = bin_begin; xbin <=bin_end; ++xbin )
    {
      for( int ybin = 0; ybin< h->GetNbinsY()+1; ++ybin )
      {
        h->SetBinContent(xbin,ybin,scales[i]);
        h->SetBinError( xbin,ybin,err[i]);
      }
    }
  }
}

void setScaleAndError(TH1D&h,vector<double> scales, vector<double> err)
{
  for(int i=0;i<n_fit_bins;i++)
  {
    int bin_begin = fit_bins_start[i], bin_end = fit_bins_end[i];
    for( int ibin = bin_begin; ibin <=bin_end; ++ibin )
    {
      h.SetBinContent(ibin,scales[i]);
      h.SetBinError( ibin,err[i]);
    }
  }
}


void setScaleAndError(TH2D&h,vector<double> scales, vector<double> err)
{
  //fitbins is here in case you want a different application (but the code currently sets this as false
  //cout << "Fit bin = " << fitBins << endl;
  for( int i = 0; i< n_fit_bins; i++ )
  {
    int bin_begin = fit_bins_start[i], bin_end = fit_bins_end[i];
    for( int xbin = bin_begin; xbin <=bin_end; ++xbin )
    {
      for( int ybin = 0; ybin< h.GetNbinsY()+1; ++ybin )
      {
        h.SetBinContent(xbin,ybin,scales[i]);
        h.SetBinError( xbin,ybin,err[i]);
      }
    }
  }
}

void setScaleAndError( TH1D *h, TH1D* scale )
{
  for( int i = 0; i<h->GetNbinsX()+1 ; i++)
  {
    //y is ptbin
    double s = scale->GetBinContent( i+1 );
    double err= scale->GetBinError( i+1 );
    h->SetBinContent(i+1, s );
    h->SetBinError(i+1, err );
  }
}

void setScaleAndError( TH2D *h, TH1D* scale )
{
  for( int i = 0; i<h->GetNbinsY()+1 ; i++)
  {
    //y is ptbin
    double s = scale->GetBinContent( i+1 );
    double err= scale->GetBinError( i+1 );
    for ( int j = 0; j< h->GetNbinsX()+1; j++ )
    {
      h->SetBinContent(j+1,i+1, s );
      h->SetBinError(j+1,i+1, err );
    }

  }
}

//==============================================
// Define getProjection
//==============================================
void getProjection(MnvH2D**h, MnvH1D**proj,bool projX)
{
  for(uint i=0;i<nHistos;i++)
  {
    string name = h[i]->GetName();
    if(projX)
    {
      name+=Form("_%d_px",i);
      proj[i] = h[i]->ProjectionX(name.c_str(),0,-1);
    }
    else
    {
      name+=Form("_%d_py",i);
      proj[i] = h[i]->ProjectionY(name.c_str(),0,-1);
    }
  }
}

vector<TH1D*> GetSingleTypeFitScalesAndErrors()
{
  vector<vector<double>> fit_scales, fit_errors;

  vector<TH1D*> scale_histos;

  for(int i = 0; i< categories_to_fit.size(); ++i)
  {
    fit_scales.push_back( vector<double>() );
    fit_errors.push_back( vector<double>() );
  }
  //start to fit
  double print_val = -1;  
  int nFitCategories = categories_to_fit.size();
  int nSidebands = Sidebands.size();
  //TFitter* minimizer = new TFitter(n_fit_bins*nFitCategories);
  TFitter* minimizer = new TFitter(100);
  minimizer->ExecuteCommand("SET PRINTOUT",&print_val,1);
  minimizer->SetFCN(minuitFunction);
  //minimizer->GetMinuit()->SetPrintLevel(1);

  double v = .8, verr=0.0000001, vmin = 0, vmax = 10;
  for(int i=0;i<n_fit_bins;i++)
  {
    //for( int j = 0; j<nFitCategories;j++) 
    //{
    //  int hist_type = categories_to_fit[j];
    //  string fit_name = names[hist_type];
    for( int j = 0; j<nFitCategories; j++ )
      minimizer->SetParameter(i*nFitCategories+j,Form("scale_%s_%d_%d",categories_to_fit_names[j].c_str(), fit_bins_start[i],fit_bins_end[i]), v,verr,vmin,vmax);
      //minimizer->SetParameter(i*nFitCategories+1,Form("scale_%s_%d_%d","pizero", fit_bins_start[i],fit_bins_end[i]), v,verr,vmin,vmax);
      //minimizer->SetParameter(i*nFitCategories+2,Form("scale_%s_%d_%d","npi", fit_bins_start[i],fit_bins_end[i]), v,verr,vmin,vmax);
    //}
  }  

  minimizer->ExecuteCommand("SIMPLEX",0,0);
  minimizer->ExecuteCommand("MIGRAD",0,0);

  for( int i = 0;i<nFitCategories; i++ ) 
  {
    int hist_type = categories_to_fit[i];
    string fit_histo_name = names[hist_type];

    for( int j = 0; j< n_fit_bins; j++ )
    {
      fit_scales[i].push_back(minimizer->GetParameter(j*nFitCategories+i));
      fit_errors[i].push_back(minimizer->GetParError(j*nFitCategories+i));
    }
    //string clone_name = Form("h_scale_%s",fit_histo_name.c_str() ); 
    TH1D* h_scale = (TH1D*) input_histos[0][kData].Clone( Form("h_scale_%s",fit_histo_name.c_str() ) );
    scale_histos.push_back( h_scale ) ;
    setScaleAndError( scale_histos[i], fit_scales[i], fit_errors[i] );
  }
  delete minimizer;


  //cout<<"Result is "<<result<<endl;
  //for( int i = 0; i< fit_scales.size(); ++ i ) 
  //{
  //  cout<<"Looking at "<<i<<endl;
  //  for( int j = 0; j< fit_scales[i].size(); j++ )
  //  cout<<fit_scales[i][j]<<endl;
  //}


  return scale_histos;
}

vector<TH1D*> GetBackgroundScales(vector<TH1D*> &scale_histos, int background_type_id)
{
  //return type vector (sample, side1 side2... )
  // samples -- background scale
  vector<TH1D*> background_scales;
  for( int i= 0; i<Samples.size(); ++i)
  {
    int sample = Samples[i];
    TH1D* background = (TH1D*) input_histos[sample][background_type_id].Clone( Form("background_angle_%02d", sample) );

    TH1D* tmp_background = (TH1D*) background->Clone( Form("tmp_%s", background->GetName()) );
    TH1D* background_scale = (TH1D*) background->Clone(Form("%s_scale",background->GetName()) );
    background_scale->Reset();
    //set errors to 0 
    for(int j = 0; j< tmp_background->GetNbinsX()+1; ++j ) 
    {
      tmp_background->SetBinError( j+1, 0 );
    }
    //now add scaled histos
    //cout<<"look at background_scale"<<endl;
    //for(int i = 0; i<background_scale->GetNbinsX()+1;i++) cout<<background_scale->GetBinContent(i+1)<<endl;
    for(int iCat = 0; iCat < categories_to_fit.size(); ++iCat)
    {
      int type_id = categories_to_fit[iCat];
      TH1D* h_type = (TH1D*) input_histos[sample][type_id].Clone( Form("tmp_%02d_%s",sample, names[type_id].c_str() ) );
      for(int j = 0; j< h_type->GetNbinsX()+1; ++j ) h_type->SetBinError( j+1, 0 );

      h_type->Multiply( scale_histos[iCat]) ;
      background_scale->Add( h_type );
    }
    //add bg components
    for(int iCat = 0; iCat < bg_categories_stay_const.size(); ++iCat)
    {
      int type_id = bg_categories_stay_const[iCat];
      TH1D* h_type = (TH1D*) input_histos[sample][type_id].Clone( Form("tmp_%s",input_histos[sample][type_id].GetName() ) );
      for(int j = 0; j< h_type->GetNbinsX()+1; ++j ) h_type->SetBinError( j+1, 0 );
      
      background_scale->Add( h_type );
    }
    //for(int i = 0; i<input_histos[sample][background_type_id].GetNbinsX()+1;i++) cout<<input_histos[sample][background_type_id].GetBinContent(i+1)<<endl;
    //cout<<"before division"<<endl;
    //for(int i = 0; i<background_scale->GetNbinsX()+1;i++) cout<<background_scale->GetBinContent(i+1)<<endl;
    //cout<<"tmp"<<endl;
    //for(int i = 0; i<tmp_background->GetNbinsX()+1;i++) cout<<tmp_background->GetBinContent(i+1)<<endl;
    
    background_scale->Divide(tmp_background);

    //cout<<"after division"<<endl;
    //for(int i = 0; i<background_scale->GetNbinsX()+1;i++) cout<<background_scale->GetBinContent(i+1)<<endl;
    
    background_scales.push_back(background_scale);

    //cout<<"background_scales.push_back(background_scale);"<<endl;

    //ratio new/old so we can then apply in cross section histos (old)*weight
  }

  return background_scales;
}




void getAllUniverseScales( vector<MnvH1D**>& tracks,vector<MnvH2D*> &track_scale )
{
  //vector<MnvH2D*> ret(nHistos);

  //CV
  for( int i = 0; i<nHistos; i++ )
  {
    for(int j=0;j<Samples.size();j++) input_histos[j][i] = *(new TH1D( tracks[j][i]->GetCVHistoWithStatError() ) );
  }

  vector<TH1D*> cv_scales_histos = GetSingleTypeFitScalesAndErrors();
  vector<TH1D*> cv_background_scale_histos = GetBackgroundScales( cv_scales_histos, kMC);

  for( int s = 0; s < Samples.size(); s++ ) setScaleAndError( track_scale[s], cv_background_scale_histos[s] );



  cout << "Fitting the Vertical errors" << endl;
  vector<string> vertNames = tracks[0][kMC]->GetVertErrorBandNames();
  cout<< vertNames.size() << endl;

  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    cout<<"Working on "<<vertNames[k] <<endl;
    MnvVertErrorBand *errBand = tracks[0][kMC]->GetVertErrorBand( vertNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH2D*>> VertErrHists;
    //Signal, Blob, Michel and MicBlob
    for( int iSample = 0; iSample < Samples.size(); iSample++)
    {
      VertErrHists.push_back( vector<TH2D*>(nUniverses) );
    }
    
    for ( int i_uni = 0; i_uni < nUniverses; i_uni++ ) //loop over all universes in a latteral error band
    {
      for( int iSample = 0; iSample < Samples.size(); iSample++)//loop over each sample
      {
        int sample = Samples[iSample];
        VertErrHists[iSample][i_uni] = (TH2D*) track_scale[0]->Clone(Form("vertband_%s_%d_%02d",vertNames[k].c_str(),i_uni, sample));
      }

      for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
      { // set input histo to the current errorband universe
        string histo_name = names[*ittype];
        for( int iH = 0; iH < Samples.size(); iH++ )
        {
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, histo_name.c_str(), vertNames[k].c_str(),i_uni)));
        }
      }
      //now, fit this universe
      //0:signal, 1:blob, 2:michel,3:micblob
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kMC);
      for( int s = 0; s < Samples.size(); s++ ) setScaleAndError( VertErrHists[s][i_uni], uni_background_scales_histos[s] );

    }



    for( int s = 0; s < Samples.size(); s++ ) track_scale[s]->AddVertErrorBand( vertNames[k], VertErrHists[s] );
    for( int s = 0; s < Samples.size(); s++ ) 
    {
      for( int i = 0; i<nUniverses;i++) delete VertErrHists[s][i];
    }

  }


  
  cout << "Fitting the Lateral errors" << endl;
  vector<string> latNames = tracks[0][kMC]->GetLatErrorBandNames();
  cout<< latNames.size() << endl;
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    cout<<"Working on "<<latNames[k] <<" ...";
    MnvLatErrorBand *errBand = tracks[0][kMC]->GetLatErrorBand( latNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH2D*>> LatErrHists;
    //Signal, Blob, Michel and MicBlob
    for( int iSample = 0; iSample < Samples.size(); iSample++)
    {
      LatErrHists.push_back( vector<TH2D*>(nUniverses) );
    }
    cout<<"LatErrHists created ...";
    
    cout<<"N_uni: "<<nUniverses<<endl;
    for ( int i_uni = 0; i_uni < nUniverses; i_uni++ )
    {
      cout<<endl<<i_uni;
      //cout<<"__set sample ";
      for( int iSample = 0; iSample < Samples.size(); iSample++)
      {
        int sample = Samples[iSample];
        LatErrHists[iSample][i_uni] = (TH2D*) track_scale[0]->Clone(Form("latband_2track_%s_%d_%02d",latNames[k].c_str(),i_uni, sample));
        //cout<<sample<<", ";
      }

      //cout<<"ittype: ";
      for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
      {
        string histo_name = names[*ittype];
        //cout<<histo_name<<", ";
        for( int iH = 0; iH < Samples.size(); iH++ )
        {
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetLatErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, histo_name.c_str(), latNames[k].c_str(),i_uni)));
        }
      }
      //now, fit this universe
      //0:signal, 1:blob, 2:michel,3:micblob
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kMC);
      //cout<<"fit: ";
      for( int s = 0; s < Samples.size(); s++ ) 
      {
        //cout<<s<<", ";
        setScaleAndError( LatErrHists[s][i_uni], uni_background_scales_histos[s] );
      }
      uni_scales_histos.clear();uni_scales_histos.shrink_to_fit();
      uni_background_scales_histos.clear();uni_background_scales_histos.shrink_to_fit();
    }
    //cout<<endl;

    for( int s = 0; s < Samples.size(); s++ ) track_scale[s]->AddLatErrorBand( latNames[k], LatErrHists[s] );
    cout<<"lat added"<<endl;
    for( int s = 0; s < Samples.size(); s++ ) 
    {
      for( int i = 0; i<nUniverses;i++) delete LatErrHists[s][i];
      LatErrHists[s].clear(); LatErrHists[s].shrink_to_fit();
    }
    LatErrHists.clear(); LatErrHists.shrink_to_fit();
    cout<<endl;
  }//end k 
}


vector<MnvH1D*> expandCategoryScaleToSample( vector<vector<MnvH1D*>>&sample_tracks, vector<MnvH1D*> & track_scale, int Category )
{
  vector<MnvH1D*> ret;
  for( int i = 0; i< Samples.size(); i++ )
  {
    string hname = Form("h_weights_q2qe_%s_%02d",names[Category].c_str(), Samples[i] );
    MnvH1D* ratio = (MnvH1D*) sample_tracks[i][Category]->Clone( hname.c_str() );
    ratio->SetTitle( hname.c_str() );
    ratio->Reset();
    //ratio->AddMissingErrorBandsAndFillWithCV( *sample_tracks[i][Category] );
    MnvH1D* tmp_hist = (MnvH1D*) sample_tracks[i][Category]->Clone("tmp") ;
    for(int iCat = 0; iCat < categories_to_fit.size(); ++iCat)
    {
      int type_id = categories_to_fit[iCat];
      MnvH1D* tmp_cat = (MnvH1D*) sample_tracks[i][type_id]->Clone("tmp_cat") ;
      tmp_cat->Multiply( tmp_cat, track_scale[iCat] );
      ratio->Add( tmp_cat );
      delete tmp_cat;
    }
    //add bg components
    for(int iCat = 0; iCat < bg_categories_stay_const.size(); ++iCat)
    {
      int type_id = categories_stay_const[iCat];
      MnvH1D* tmp_cat = (MnvH1D*) sample_tracks[i][type_id]->Clone("tmp_cat") ;
      ratio->Add( tmp_cat );
      delete tmp_cat;
    }
    ratio->Divide( ratio, tmp_hist );
    ret.push_back(ratio);
    delete tmp_hist;
  }
  return ret;
}


void getAllUniverseScales1D( vector<vector<MnvH1D*>> &tracks,vector<MnvH1D*> &track_scale )
{
  //vector<MnvH1D*> ret(nHistos);

  //CV
  for( int i = 0; i<Samples.size(); i++ )
  {
    vector<TH1D*> hists;
    for( int j = 0; j< nHistos;j++ )
    {
      input_histos[i][j] = *(new TH1D( tracks[i][j]->GetCVHistoWithStatError() ));
      input_histos[i][j].SetName( Form("%s_%02d_%02d", input_histos[i][j].GetName(), i, j ) );
      cout<<i<<" "<<j<<endl;
      hists.push_back( new TH1D( input_histos[i][j] ) );
    }
    string basename= Form("%02d",i );
    DrawNHists( hists,basename);
    //delete[] &hists[0];
  }

  //for( int i = 0; i<nHistos; i++ )
  //{
  //  for(int j = 0; j<Samples.size();++j) 
  //  {
  //    input_histos[j][i] = *(new TH1D( tracks[j][i]->GetCVHistoWithStatError() ));
  //    input_histos[j][i].SetName( Form("%s_%02d_%02d", input_histos[j][i].GetName(), j, i ) );
  //  }
  //}

  vector<TH1D*> cv_scales_histos = GetSingleTypeFitScalesAndErrors();// scales for pizero, piplus and multipi
  vector<TH1D*> cv_background_scale_histos = GetBackgroundScales( cv_scales_histos, kMC);

  //setting cv bg scale to 1D CV
  for( int s = 0; s < Samples.size(); s++ ) setScaleAndError( track_scale[s], cv_background_scale_histos[s] );


  if(!doErrorBands) return;

  cout << "Fitting the Vertical errors" << endl;
  vector<string> vertNames = tracks[0][kMC]->GetVertErrorBandNames();
  cout<< vertNames.size() << endl;

  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    cout<<"Working on "<<vertNames[k] <<endl;
    MnvVertErrorBand *errBand = tracks[0][kMC]->GetVertErrorBand( vertNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH1D*>> VertErrHists;
    //Signal, Blob, Michel and MicBlob
    for( int iSample = 0; iSample < Samples.size(); iSample++)
    {
      VertErrHists.push_back( vector<TH1D*>(nUniverses) );
    }
    
    for ( int i_uni = 0; i_uni < nUniverses; i_uni++ ) //loop over all universes in a vertical error band
    {
      for( int iSample = 0; iSample < Samples.size(); iSample++)//loop over each sample
      {
        int sample = Samples[iSample];
        VertErrHists[iSample][i_uni] = (TH1D*) track_scale[0]->Clone(Form("vertband_%s_%d_%02d",vertNames[k].c_str(),i_uni, sample));
      }

      for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
      { // set input histo to the current errorband universe
        string histo_name = names[*ittype];
        for( int iH = 0; iH < Samples.size(); iH++ )
        {
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, names[*ittype].c_str(), vertNames[k].c_str(),i_uni)));
        }
      }

      for( int iH = 0; iH < Samples.size(); iH++ )
      {
        cout<<Form("%s has %d, %d events", input_histos[iH][kQELike_QE_OTH].GetName(), (int) input_histos[iH][kQELike_QE_OTH].Integral(), (int) tracks[iH][kQELike_QE_OTH]->GetVertErrorBand( vertNames[k] )->Integral() )<<endl;
        cout<<Form("%s has %d, %d events", input_histos[iH][kQELike_RES].GetName(),    (int) input_histos[iH][kQELike_RES].Integral(),    (int) tracks[iH][kQELike_RES]->GetVertErrorBand( vertNames[k] )->Integral() )<<endl;
        cout<<Form("%s has %d, %d events", input_histos[iH][kQELike_2p2h].GetName(),   (int) input_histos[iH][kQELike_2p2h].Integral(),   (int) tracks[iH][kQELike_2p2h]->GetVertErrorBand( vertNames[k] )->Integral() )<<endl;
      }

      //now, fit this universe
      //0:signal, 1:blob, 2:michel,3:micblob
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kMC);
      for( int s = 0; s < Samples.size(); s++ ) setScaleAndError( VertErrHists[s][i_uni], uni_background_scales_histos[s] );

    }
    for( int s = 0; s < Samples.size(); s++ ) track_scale[s]->AddVertErrorBand( vertNames[k], VertErrHists[s] );
    for( int s = 0; s < Samples.size(); s++ ) 
    {
      for( int i = 0; i<nUniverses;i++) delete VertErrHists[s][i];
    }

  }




  cout << "Fitting the Lateral errors" << endl;
  vector<string> latNames = tracks[0][kMC]->GetLatErrorBandNames();
  cout<< latNames.size() << endl;
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    cout<<"Working on "<<latNames[k] <<endl;
    MnvLatErrorBand *errBand = tracks[0][kMC]->GetLatErrorBand( latNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH1D*>> LatErrHists;
    //Signal, Blob, Michel and MicBlob
    for( int iSample = 0; iSample < Samples.size(); iSample++)
    {
      LatErrHists.push_back( vector<TH1D*>(nUniverses) );
    }
    
    for ( int i_uni = 0; i_uni < nUniverses; i_uni++ )
    {
      for( int iSample = 0; iSample < Samples.size(); iSample++)
      {
        int sample = Samples[iSample];
        LatErrHists[iSample][i_uni] = (TH1D*) track_scale[0]->Clone(Form("latband_2track_%s_%d_%02d",latNames[k].c_str(),i_uni, sample));
      }

      for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
      {
        string histo_name = names[*ittype];
        for( int iH = 0; iH < Samples.size(); iH++ )
        {
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, histo_name.c_str(), latNames[k].c_str(),i_uni)));
        }
      }
      //now, fit this universe
      //0:signal, 1:blob, 2:michel,3:micblob
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kMC);
      for( int s = 0; s < Samples.size(); s++ ) setScaleAndError( LatErrHists[s][i_uni], uni_background_scales_histos[s] );

    }
    for( int s = 0; s < Samples.size(); s++ ) track_scale[s]->AddLatErrorBand( latNames[k], LatErrHists[s] );
    for( int s = 0; s < Samples.size(); s++ ) 
    {
      for( int i = 0; i<nUniverses;i++) delete LatErrHists[s][i];
      LatErrHists[s].clear(); LatErrHists[s].shrink_to_fit();
    }
    LatErrHists.clear(); LatErrHists.shrink_to_fit();
  }
}

void getAllUniverseScalesSingleType1D( vector<vector<MnvH1D*>> &tracks, vector<MnvH1D*> &track_scale)
{
  int nTypes = bg_categories_to_scale.size();
  cout<<"getAllUniverseScalesSingleType1D"<<endl;
  // updating the global histo holder

  //TCanvas c("c","c",300,300);
  //c.cd();
  //
  //
  cout<<"track size: "<<tracks.size()<<endl;
  for( int i = 0; i<tracks.size(); i++ )
  {
    cout<<i<<endl;
    for(int j = 0; j < nHistos; ++j ) 
    {
      cout<<i<<", "<<j<<endl;
      cout<<tracks[i][j]->GetName()<<"---"<<names[j]<<endl;
      input_histos[i][j] = *(new TH1D( tracks[i][j]->GetCVHistoWithStatError() ));
      cout<<Form("%02d %02d %d", i,j,(int)tracks[i][j]->Integral() )<<endl;
    }
  }

  cout<<"Transferred cv histos"<<endl;
  vector<TH1D*> cv_scales_histos = GetSingleTypeFitScalesAndErrors();// scales for pizero, piplus and multipi
  for( int t = 0; t< nTypes ; t++) 
  {
    setScaleAndError( track_scale[t], cv_scales_histos[t] );
    //cv_scales_histos[t]->Draw("histe");
    //c.Print( Form("./test_%02d.png", t ));
  }
  cout<<"Got CV"<<endl;


  if(!doErrorBands) return;
  //Vert Error
  cout<<"Working on Vert Error Bands"<<endl;
  vector<string> vertNames = tracks[0][kMC]->GetVertErrorBandNames();
  for( unsigned int k= 0; k < vertNames.size(); ++k )
  {
    cout<<vertNames[k]<<endl;
    MnvVertErrorBand *errBand = tracks[0][kMC]->GetVertErrorBand( vertNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH1D*>> VertErrHists;// Type -- Universe
    for( int iCat = 0; iCat < nTypes; iCat++ ) VertErrHists.push_back( vector<TH1D*>(nUniverses) );
    //loop over each universe in vert err band
    for( int i_uni = 0; i_uni < nUniverses; i_uni++ )
    {
      //loop over each  type
      for ( int it = 0; it < nTypes; it++ )
      {
        VertErrHists[it][i_uni] = (TH1D*) track_scale[0]->Clone(Form("vertband_%s_%d_%d", vertNames[k].c_str(), bg_categories_to_scale[it], i_uni) );
      }

      //for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
      for( int iH = 0; iH <Samples.size(); iH++)
      {
        //TCanvas c("c","c",500,500);
        input_histos[iH][kData] = *((TH1D*) tracks[iH][kData]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, "data", vertNames[k].c_str(),i_uni)));
        input_histos[iH][kQELike] = *((TH1D*) tracks[iH][kQELike]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, "qelike", vertNames[k].c_str(),i_uni)));
        //input_histos[iH][kData].Draw("hist");
        for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
        { 
          string histo_name = names[*ittype];
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, names[*ittype].c_str(), vertNames[k].c_str(),i_uni)));
          //input_histos[iH][*ittype].Draw("histlsame");
          //cout<<Form("vert_2track_angle_%02d_%s_%s_%d integral: %d", iH, names[*ittype].c_str(), vertNames[k].c_str(), i_uni, (int) input_histos[iH][*ittype].Integral() )<<endl;
        }

        //cout<<Form("%s integral: %d", input_histos[iH][kData].GetName(), (int) input_histos[iH][kData].Integral() )<<endl;
        //cout<<Form("%s integral: %d", input_histos[iH][kQELike].GetName(), (int) input_histos[iH][kQELike].Integral() )<<endl;
        //cout<<Form("%s integral: %d", input_histos[iH][kQELike_RES].GetName(), (int) input_histos[iH][kQELike_RES].Integral() )<<endl;
        //cout<<Form("vert_1track_angle_%02d_%s_%s_%d integral: %d", iH, "data", vertNames[k].c_str(), i_uni, (int) input_histos[iH][kData].Integral() )<<endl;
        //cout<<Form("vert_2track_angle_%02d_%s_%s_%d integral: %d", iH, "qelike", vertNames[k].c_str(), i_uni, (int) input_histos[iH][kQELike].Integral() )<<endl;
        //c.Print(Form("plots/vert_2track_angle_%02d_%s_%s_%d.pdf", iH, "data", vertNames[k].c_str(), i_uni));
      }

      //now fit this universe
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      for( int t = 0; t<nTypes; t++) setScaleAndError( VertErrHists[t][i_uni], uni_scales_histos[t] );
    }//done looping through all uni in a sideband
    //push sideband 
    for( int t = 0; t<nTypes;t++ ) 
    {
      track_scale[t]->AddVertErrorBand( vertNames[k], VertErrHists[t] );
      for( int i = 0; i< nUniverses; i++ ) delete VertErrHists[t][i];
    }
  }// done Vert Err


  //Lat Error
  vector<string> latNames = tracks[0][kMC]->GetLatErrorBandNames();
  cout<<"Working on Lat Error Bands"<<endl;
  for( unsigned int k= 0; k < latNames.size(); ++k )
  {
    cout<<latNames[k]<<endl;
    MnvLatErrorBand *errBand = tracks[0][kMC]->GetLatErrorBand( latNames[k] );
    int nUniverses = errBand->GetNHists();
    vector<vector<TH1D*>> LatErrHists;// Type -- Universe
    for( int iCat = 0; iCat < nTypes; iCat++ ) LatErrHists.push_back( vector<TH1D*>(nUniverses) );
    //loop over each universe in lat err band
    for( int i_uni = 0; i_uni < nUniverses; i_uni++ )
    {
      //loop over each  type
      for ( int it = 0; it < nTypes; it++ )
      {
        LatErrHists[it][i_uni] = (TH1D*) track_scale[0]->Clone(Form("latband_%s_%d_%d", latNames[k].c_str(), bg_categories_to_scale[it], i_uni) );
      }

      for( int iH = 0; iH <Samples.size(); iH++)
      {
        input_histos[iH][kData] = *((TH1D*) tracks[iH][kData]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, "data", latNames[k].c_str(),i_uni)));
        input_histos[iH][kQELike] = *((TH1D*) tracks[iH][kQELike]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, "qelike", latNames[k].c_str(),i_uni)));
        //cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, "data", latNames[k].c_str(), i_uni, (int) input_histos[iH][kData].Integral() )<<endl;
        //cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, "qelike", latNames[k].c_str(), i_uni, (int) input_histos[iH][kQELike].Integral() )<<endl;
        for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
        { 
          string histo_name = names[*ittype];
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, names[*ittype].c_str(), latNames[k].c_str(),i_uni)));
          //cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, names[*ittype].c_str(), latNames[k].c_str(), i_uni, (int) input_histos[iH][*ittype].Integral() )<<endl;
        }
      }

      //now fit this universe
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      for( int t = 0; t<nTypes; t++) setScaleAndError( LatErrHists[t][i_uni], uni_scales_histos[t] );
    }//done looping through all uni in a sideband
    //push sideband 
    for( int t = 0; t<nTypes;t++ ) 
    {
      track_scale[t]->AddLatErrorBand( latNames[k], LatErrHists[t] );
      for( int i = 0; i< nUniverses; i++ ) delete LatErrHists[t][i];
    }
  }// done Lat Err
}
//////////////////////////////
// This is the main function//
//////////////////////////////

void initializeVector( vector<MnvH1D**> &h){
  const int nelements = h.size();
  for (int i=0; i<nelements; ++i)
  {
    h[i] = new MnvH1D*[nHistos];
  }
}

void initializeVector( vector<MnvH2D**> &h){
  const int nelements = h.size();
  for (int i=0; i<nelements; ++i)
  {
    h[i] = new MnvH2D*[nHistos];
  }
}

int SideBandFit( vector<string> filenames )
{
  //-----------------------------------------
  // Read input files to get histograms
  //-----------------------------------------
  //File order = signal 2 , blob side 2, michel side 2 , micblob side 2 , outputfile
  const int n_samples = Samples.size();
  

  cout << "Loading file: " << filenames[1] << endl;
  TFile* f_signal = new TFile( filenames[1].c_str() , "READ" );
  TFile* f_blob = new TFile( filenames[2].c_str() , "READ" );
  TFile* f_michel = new TFile( filenames[3].c_str() , "READ" );
  TFile* f_micblob = new TFile( filenames[4].c_str() , "READ" );
  if (f_signal->IsZombie() || f_signal->GetListOfKeys()->IsEmpty())
  {
    Error("SideBandFit","Could not get Signal histogram ROOT file or it was empty.");
    return 1;
  }
  if (f_blob->IsZombie() || f_blob->GetListOfKeys()->IsEmpty())
  {
    Error("SideBandFit","Could not get BlobSideBand histogram ROOT file or it was empty.");
    return 1;
  }
  if (f_michel->IsZombie() || f_michel->GetListOfKeys()->IsEmpty())
  {
    Error("SideBandFit","Could not get MichelSideBand histogram ROOT file or it was empty.");
    return 1;
  }
  if (f_micblob->IsZombie() || f_micblob->GetListOfKeys()->IsEmpty())
  {
    Error("SideBandFit","Could not get MicBlobSideBand histogram ROOT file or it was empty.");
    return 1;
  }
 

  vector<string> multiplicity;
  multiplicity.push_back("2_track");
  //check that each is what it is supposed to be (sample, track)
  //TMap *sampleMap = (TMap*)fin->Get("sample");
  //TObjString *multiplicity_obj = (TObjString*)sampleMap->GetValue("multiplicity");
  //TString multiplicity_label = multiplicity_obj->GetString();    
  //TObjString *sample_obj = (TObjString*)sampleMap->GetValue("sample");
  //TString sample = sample_obj->GetString();
  //cout<<"Samples: "<<sample<<endl;
  //cout<<"Multiplicty: "<<multiplicity_label<<endl;
  //
  //delete sampleMap;
  //delete multiplicity_obj;
  //delete sample_obj;


  //Get QELikenot weights
  
  bool fluxHistoExists = true; //Hack, Previous loop determined all 6 input files have flux histo
  cout << "SideBandFit: fluxHistoExists " << fluxHistoExists << endl; 

  //Get the tools
  CCQENuPlotUtils *plotutils    = new CCQENuPlotUtils( fluxHistoExists );
  //CCQENuUtils     *histoutils   = new CCQENuUtils();
  //AnaBinning      *binner       = new AnaBinning();
  //CCQENuBinning   *minmodbinner = new CCQENuBinning();
  //NeutronBlobBinning *neutbinner = new NeutronBlobBinning();
  CCQENuCuts      *cutter       = new CCQENuCuts();
  MnvPlotter      *plotter      = plotutils->getPlotter();
  //CCQENuUtils     *utils        = new CCQENuUtils( false, fluxHistoExists );

  //--------------------------------------------------
  // Get PT, PZ binning for the weights histograms 
  //--------------------------------------------------
  cout << "Done with setting up general functions" << endl;
  
  //---------------------------------------------------------------
  // set which bins to use
  //--------------------------------------------------------------
  //axis_binning ybins = muonPtbins;
  axis_binning ybins = Q2bins;
  if (fit_xvar == "h_ptmu" ) ybins = muonPtbins;

  //---------------------------------------------------------------
  // Create vector of pointers to pointers
  // Each element of vector is a pointer to an array of MnvH1D 
  // sample - categories
  //---------------------------------------------------------------
  // vector < vector < MnvH1D*[nHistos]> > 
  vector<vector<MnvH1D*>>  h_fit_varbins_samples;
  for( int i = 0; i<n_samples; i++ ) 
  {
    h_fit_varbins_samples.push_back( vector<MnvH1D*>(nHistos) );
  }
  cout<<"h_fit_varbins_samples size: "<<h_fit_varbins_samples.size()<<endl;

  //signal, blob, michel, micblob
  //for( int i = 0; i<n_samples; i++ )
  //{
  //  vector<MnvH1D**> vec( n_fitvarbins_2track );
  //  initializeVector( vec ); 
  //  h_fit_varbins_samples.push_back( vec );
  //}

  double pot_data = plotutils->getPOTData( f_signal );
  double pot_mc = plotutils->getPOTMC( f_signal );
  double pot_scale = plotutils->getPOTNormFactor( f_signal );
  
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "---------------------------------------"  << endl;
  cout<< "POT INFORMATION:"  << endl;
  cout<< "POT Data = " << pot_data << endl;
  cout<< "POT MC   = " << pot_mc << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;

  for( int j = 0; j < nSamplesPerSideBand; j++ )
  {
    cout<<"Getting Region "<<j<<endl;
    cout<<"signal"<<endl;
    plotutils->bookHistos( f_signal, &(h_fit_varbins_samples[j][0]), fit_xvars[j] );      
    cout<<"blob"<<endl;
    plotutils->bookHistos( f_blob, &(h_fit_varbins_samples[j+nSamplesPerSideBand][0]), fit_xvars[j] );      
    cout<<"michel"<<endl;
    plotutils->bookHistos( f_michel, &(h_fit_varbins_samples[j+nSamplesPerSideBand*2][0]), fit_xvars[j] );      
    cout<<"micblob"<<endl;
    plotutils->bookHistos( f_micblob, &(h_fit_varbins_samples[j+nSamplesPerSideBand*3][0]), fit_xvars[j] );      
    cout<<h_fit_varbins_samples[j][kData]<<endl;
    cout<<"sample name: "<<h_fit_varbins_samples[j][kData]->GetName()<<endl;
    cout<<"sample bins: "<<h_fit_varbins_samples[j][kData]->GetNbinsX()<<endl;
    cout<<"sample name: "<<h_fit_varbins_samples[j+nSamplesPerSideBand][kData]->GetName()<<endl;
    cout<<"sample bins: "<<h_fit_varbins_samples[j+nSamplesPerSideBand][kData]->GetNbinsX()<<endl;
    cout<<"sample name: "<<h_fit_varbins_samples[j+nSamplesPerSideBand*2][kData]->GetName()<<endl;
    cout<<"sample bins: "<<h_fit_varbins_samples[j+nSamplesPerSideBand*2][kData]->GetNbinsX()<<endl;
    cout<<"sample name: "<<h_fit_varbins_samples[j+nSamplesPerSideBand*3][kData]->GetName()<<endl;
    cout<<"sample bins: "<<h_fit_varbins_samples[j+nSamplesPerSideBand*3][kData]->GetNbinsX()<<endl;
    plotutils->scaleMCHistos( &h_fit_varbins_samples[j][0], pot_scale);
    plotutils->scaleMCHistos( &h_fit_varbins_samples[j+nSamplesPerSideBand][0], pot_scale);
    plotutils->scaleMCHistos( &h_fit_varbins_samples[j+nSamplesPerSideBand*2][0], pot_scale);
    plotutils->scaleMCHistos( &h_fit_varbins_samples[j+nSamplesPerSideBand*3][0], pot_scale);

    MnvH1D *hMCForSideBands = h_fit_varbins_samples[j][kMC];
    h_fit_varbins_samples[j][kData]->AddMissingErrorBandsAndFillWithCV( *hMCForSideBands );
    h_fit_varbins_samples[j+nSamplesPerSideBand][kData]->AddMissingErrorBandsAndFillWithCV( *hMCForSideBands );
    h_fit_varbins_samples[j+nSamplesPerSideBand*2][kData]->AddMissingErrorBandsAndFillWithCV( *hMCForSideBands );
    h_fit_varbins_samples[j+nSamplesPerSideBand*3][kData]->AddMissingErrorBandsAndFillWithCV( *hMCForSideBands );



    if(normalizeData)
    {
      for( int k = 0; k < 4; k++ ) 
      {
        int iRegion = j+nSamplesPerSideBand*k;
        double ndata = h_fit_varbins_samples[iRegion][kData]->Integral(0,-1);
        double ndata0 = h_fit_varbins_samples[0][kData]->Integral(0,-1);
        RegionWeights[iRegion] = ndata0/ndata;
        //for( int l = 0; l < nHistos; l++) h_fit_varbins_samples[iRegion][l]->Scale(1./ndata,"", true);
      }
    }

    //cout<< "Subtracting Data"<<endl;
    //MnvH1D* h_bck = (MnvH1D*) h_fit_varbins_samples[j][kQELikeNot]->Clone("h_bck");
    //h_bck->Multiply( h_bck, h_bck_weights );
    //cout<<j<<" has "<<h_fit_varbins_samples[j][kData]->Integral()<<" original data"<<endl;
    //h_fit_varbins_samples[j][kData]->Add( h_bck, -1 );
  
  }
  //for( int i = 0; i<n_samples;i++)
  //{
  //  //cout<<"Hists "<<i<<"==========="<<endl;
  //  //for( int j = 0; j< nHistos; j++ ) cout<<Form("fill %02d %s %d", i, names[j].c_str(), (int) h_fit_varbins_samples[i][j]->Integral() )<<endl;
  //  //PrintMnvH1D(h_fit_varbins_samples[i][kQELike_QE_OTH] );
  //}


  cout << "Done with Tejin's track" << endl;
  f_signal->Close();
  f_blob->Close();
  f_michel->Close();
  f_micblob->Close();
  cout << "Done booking histograms" << endl;


  //vector<MnvH1D*> h_vect_yvar_scales;
  //vector<MnvH2D*> h_vect_dthetaPdthetaR_scales;

  //  for( int i = 0; i<n_samples; i++ ) 
  //{
  //  string hname = Form( "h_weights_dthetaPdthetaR_yvarbins_qelikenot_%s", Samples[i].c_str() ) ;
  //  h_vect_dthetaPdthetaR_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), global_x.nbins, &(global_x.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );

  //}

  //vector<MnvH2D*> h_dalphat_scales, h_dthetaR_scales, h_dthetaP_scales;
  //getAllUniverseScales1D(h_fit_varbins_samples, h_vect_yvar_scales);
  
  //ExpandHistos( h_vect_yvar_scales, h_vect_dthetaP_scales, 1);
  cout << "*********************************" <<endl;
  cout << "GetSingleType" << endl;
  cout << "*********************************" << endl;
  vector<MnvH1D*> hs_vect_yvar_scales;
  vector<MnvH2D*> hs_vect_dthetaPpos_scales;
  vector<MnvH2D*> hs_vect_dthetaPdthetaR_scales;

  vector<MnvH2D*> hs_vect_dpt_scales;
  vector<MnvH2D*> hs_vect_dptx_scales;
  vector<MnvH2D*> hs_vect_dpty_scales;

  vector<MnvH2D*> hs_vect_dalphat_scales;
  vector<MnvH2D*> hs_vect_dphit_scales;
  vector<MnvH2D*> hs_vect_pn_scales;

  vector<MnvH2D*> hs_vect_ptheta_scales;
  vector<MnvH2D*> hs_vect_dthetaP_scales;
  vector<MnvH2D*> hs_vect_dthetaR_scales;

  vector<MnvH2D*> hs_vect_recoil_scales;
  vector<MnvH2D*> hs_vect_dtheta2D_scales;


  vector<MnvH2D*> hs_vect_blobDist_scales;
  vector<MnvH2D*> hs_vect_blobEnergy_scales;
  vector<MnvH2D*> hs_vect_nBlobs_scales;



  axis_binning hdlx = dthetaPerpbins; 
  axis_binning hdlz = dthetaReactbins;
  axis_binning global_x;
  GetGlobal_x( hdlx, ybins, hdlz, global_x );



  int nTypes = bg_categories_to_scale.size();
  for( int i = 0; i<nTypes; i++ ) 
  {
    string typeName = bg_categories_to_scale_names[i];
    string hname;
    axis_binning xbins;

    hname = Form("hs_weights_yvar_bgType_%s", typeName.c_str() );
    hs_vect_yvar_scales.push_back( new MnvH1D( hname.c_str(), hname.c_str(), ybins.nbins, &(ybins.bin_edges[0]) ) ) ;

    xbins = dptbins;
    hname = Form( "hs_weights_dpt_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dpt_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dptxbins;
    hname = Form( "hs_weights_dptx_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dptx_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dptybins;
    hname = Form( "hs_weights_dpty_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dpty_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dalphatbins;
    hname = Form( "hs_weights_dalphat_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dalphat_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dphitbins;
    hname = Form( "hs_weights_dphit_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dphit_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = pnbins;
    hname = Form( "hs_weights_pn_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_pn_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

   
    xbins = dthetaReactbins;
    hname = Form( "hs_weights_dthetaR_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dthetaR_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dthetaPerpbins;
    hname = Form( "hs_weights_dthetaP_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dthetaP_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = dthetaPerpbinsPositive;
    hname = Form( "hs_weights_dthetaPpos_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_dthetaPpos_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    xbins = protonThetabins;
    hname = Form( "hs_weights_ptheta_yvarbins_bgType_%s", typeName.c_str() ) ;
    hs_vect_ptheta_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    
    hname = Form( "hs_weights_dthetaPdthetaR_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = global_x;
    hs_vect_dthetaPdthetaR_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );


    hname = Form( "hs_weights_recoil_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = UniformVisibleEBins;
    hs_vect_recoil_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );

    hname = Form( "hs_weights_dtheta2D_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = UniformDThetaBins;
    hs_vect_dtheta2D_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );

    hname = Form( "hs_weights_blobDist_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = BlobDistBins;
    hs_vect_blobDist_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );

    hname = Form( "hs_weights_blobEnergy_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = BlobEnergyBins;
    hs_vect_blobEnergy_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );

    hname = Form( "hs_weights_nBlobs_yvarbins_bgType_%s", typeName.c_str() ) ;
    xbins = BlobNumBins;
    hs_vect_nBlobs_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );



  }


  //vector<MnvH2D*> hs_dalphat_scales, hs_dthetaR_scales, hs_dthetaP_scales;

  getAllUniverseScalesSingleType1D(h_fit_varbins_samples, hs_vect_yvar_scales);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaPdthetaR_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaPpos_scales, 1);

  ExpandHistos( hs_vect_yvar_scales, hs_vect_ptheta_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaP_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaR_scales, 1);

  ExpandHistos( hs_vect_yvar_scales, hs_vect_dalphat_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_pn_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dphit_scales, 1);

  ExpandHistos( hs_vect_yvar_scales, hs_vect_dpt_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dptx_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dpty_scales, 1);

  ExpandHistos( hs_vect_yvar_scales, hs_vect_recoil_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_dtheta2D_scales, 1);

  ExpandHistos( hs_vect_yvar_scales, hs_vect_blobDist_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_blobEnergy_scales, 1);
  ExpandHistos( hs_vect_yvar_scales, hs_vect_nBlobs_scales, 1);







  cout<< " Done with single category "<<endl;
  cout<< " Expanding To Full Category: "<<names[FitCateogry]<<endl;
  vector<MnvH1D*> h_vect_yvar_scales = expandCategoryScaleToSample( h_fit_varbins_samples, hs_vect_yvar_scales, FitCateogry );
  //vector<MnvH2D*> hs_vect_dthetaPdthetaR_scales;
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaPdthetaR_scales, 1);

  cout << "*********************************" << endl;
  cout << "Done Saving Ouput" << endl;
  cout << "*********************************" << endl;
  TFile *tmp = new TFile(filenames[5].c_str(),"RECREATE"); //micblob
  tmp->cd();
  for(int i=0;i<n_samples;i++){
    int iR = i;
    //if(i==n_samples-1) iR=99;
    cout<<"Saving... "<<Samples[i]<<endl;
    //h_vect_yvar_scales[i]->Write();
    h_vect_yvar_scales[i]->Write();
    MnvH1D *h = h_fit_varbins_samples[i][FitCateogry]->Clone( Form("h_q2qe_%02d_%s_orig",iR,names[FitCateogry].c_str() ) );
    MnvH1D *data = h_fit_varbins_samples[i][kData]->Clone( Form("h_q2qe_%02d_data",iR ) );
    h->Write();
    h->Multiply(h, h_vect_yvar_scales[i] );
    h->Write( Form("h_q2qe_%02d_%s_fit",iR,names[FitCateogry].c_str() ) );
    data->Write();
    delete data;
    delete h;
    
    //h_vect_dthetaP_scales[i]->Write();
    //h_vect_dthetaR_scales[i]->Write();
    //h_vect_dalphat_scales[i]->Write();
    //h_vect_pn_scales[i]->Write();
    //h_vect_dpt_scales[i]->Write();
    //h_vect_dthetaPq2qe_scales[i]->Write();
    //h_vect_dthetaRq2qe_scales[i]->Write();
  }

  for(int i=0;i<nTypes;i++){
    cout<<"Saving... "<<bg_categories_to_scale_names[i]<<endl;
    hs_vect_yvar_scales[i]->Write();
    hs_vect_dthetaPdthetaR_scales[i]->Write();
    hs_vect_dthetaPpos_scales[i]->Write();
                                   
    hs_vect_dpt_scales[i]->Write();
    hs_vect_dptx_scales[i]->Write();
    hs_vect_dpty_scales[i]->Write();
                                   
    hs_vect_dalphat_scales[i]->Write();
    hs_vect_dphit_scales[i]->Write();
    hs_vect_pn_scales[i]->Write();
                                   
    hs_vect_ptheta_scales[i]->Write();
    hs_vect_dthetaP_scales[i]->Write();
    hs_vect_dthetaR_scales[i]->Write();

    hs_vect_recoil_scales[i]->Write();
    hs_vect_dtheta2D_scales[i]->Write();

    hs_vect_blobDist_scales[i]->Write();
    hs_vect_blobEnergy_scales[i]->Write();
    hs_vect_nBlobs_scales[i]->Write();



  }

  tmp->Close();
  delete plotutils;
  delete minmodbinner; 
  delete cutter; 
  delete plotter; 
  delete neutbinner;

  return 0;
}

int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./SideBandFit  Name_and_path_to_SignalSelectedEvents_file  Name_and_path_to_BlobSideBandSelectedEvents_file Name_and_path_to_MichelSideBandeSelectedEvents_file \n\n"<<
      "\t-Name_and_path_to_background_templates_file\t =\t Name of and path to MuonSelection selecting the signal events 2 track\n"<<
      "\t-Name_and_path_to_Output_file\t =\t Name of and path to the Output ROOT file that will be created \n"<<
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("SideBandFit");
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_Signal.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_BlobSideBand.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_MichelSideBand.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/MuonEvent_MicBlobSideBand.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/SignalFitWeightsAngle.root",getenv("CCQENUROOT") ) ); //output file

  const int nArgs = par.size();

  //! Set user parameters
  for( int i=0; i<nArgs; ++i){
    par.at(i) = argv[i];
  }
  if( argc > nArgs ) 
  {
    highweight = (int) std::stoi( argv[nArgs] );
    if(argc > nArgs+1) lambda = (double) std::stod( argv[nArgs+1] );
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  std::cout<<"highweight =  "<< highweight << std::endl;
  std::cout<<"lambda =  "<< lambda << std::endl;

  
  return SideBandFit(par);
}

