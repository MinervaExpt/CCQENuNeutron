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

bool doErrorBands = true;
//Global histograms because TMinuit needs them global...



//string fit_xvar = "h_ptmu";
//const int n_fit_bins = 9;
//int fit_bins_start[n_fit_bins] = {1,2,3,4,5,6,7,8,14};
//int fit_bins_end[n_fit_bins] = {1,2,3,4,5,6,7,13,14};

string fit_xvar = "h_q2qe";
vector<string> fit_xvars({ "h_q2qe_region_00", 
    "h_q2qe_region_01", 
    "h_q2qe_region_02", 
    "h_q2qe_region_03", 
    "h_q2qe_region_04", 
    "h_q2qe_region_05"});

//1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18
//
//const int n_fit_bins = 5;
//int fit_bins_start[n_fit_bins] = {1,5,7,9,11};
//int fit_bins_end[n_fit_bins] =   {4,6,8,10,18};
//const int n_fit_bins = 1;
//int fit_bins_start[n_fit_bins] = {1};
//int fit_bins_end[n_fit_bins] =   {18};
//const int n_fit_bins = 2;
//int fit_bins_start[n_fit_bins] = {1,9};
//int fit_bins_end[n_fit_bins] =   {8,18};
//const int n_fit_bins = 3;
//int fit_bins_start[n_fit_bins] = {1,6,13};
//int fit_bins_end[n_fit_bins] =   {5,12,18};


const int n_fit_bins = 15;
int fit_bins_start[n_fit_bins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
int fit_bins_end[n_fit_bins] =   {1,2,3,4,5,6,7,8,9,10,11,12,13,14,18};


vector<int> Samples({0,1,2,3,4,5});
vector<int> Sidebands({0,1,2, 3,4,5});//different angular bins

//categories_to_fit: histos with scaling, types_to_stay_const, no scaling, needed to make kMC = sum(fit+const)
vector<int> categories_to_fit({ kQELike_QE_OTH, kQELike_2p2h, kQELike_RES});
vector<string> categories_to_fit_names({"qe_oth", "2p2h","res"});
vector<int> categories_stay_const({ kQELike_DIS,kQELike_OTH  });

//only the background components included here. These histos makes up kQELikeNot
vector<int> bg_categories_to_scale = categories_to_fit;
vector<string> bg_categories_to_scale_names=categories_to_fit_names;
vector<int> bg_categories_stay_const= categories_stay_const;

vector<int> all_universe_histos({ kQELike_QE_OTH, kQELike_RES, kQELike_2p2h, kQELike_DIS, kQELike_OTH });

vector<TH1D*> input_histos(1+5, new TH1D[nHistos]);

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
vector<TH1D*> GetBackgroundScales(vector<TH1D*> &scale_histos, int background_type_id=kQELikeNot);


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
    cout<<"VertErr: "<<name<<endl;
    unsigned int nUniverses = h->GetVertErrorBand( name )->GetNHists();
    for( unsigned int iUniv = 0 ; iUniv<nUniverses; iUniv++ )
      cout<<"Univ "<<iUniv<<": "<<h->GetVertErrorBand( name )->GetHist(iUniv)->Integral()<<endl;
  }
  vector<string> latNames = h->GetLatErrorBandNames();
  for( auto name : latNames )
  {
    cout<<"LatErr: "<<name<<endl;
    unsigned int nUniverses = h->GetLatErrorBand( name )->GetNHists();
    for( unsigned int iUniv = 0 ; iUniv<nUniverses; iUniv++ )
      cout<<"Univ "<<iUniv<<": "<<h->GetLatErrorBand( name )->GetHist(iUniv)->Integral()<<endl;
  }
}

double chi2Function(double *scale )
{
  double chi2 = 0;
  double chi2_shape = 0;
  int nSidebands = Sidebands.size();
  int nFitCategories = categories_to_fit.size();
  double mc = 0;
  double data = 0;
  //cout<<"data chi2 = "<<endl;
  for( int i = 0; i< n_fit_bins; i++ )
  {
    for( vector<int>::iterator itSideband = Sidebands.begin(); itSideband!=Sidebands.end();++itSideband )
    {
      //cout<<*itSideband<<endl;
      //index of sideband sample
      int iSideband = itSideband - Sidebands.begin();
      double sum = 0;
      //loop through the fitting histos, with scaling
      for( vector<int>::iterator itCat = categories_to_fit.begin(); itCat!=categories_to_fit.end();++itCat )
      {
        //cout<<"Histo Type: "<<names[*itCat]<<endl;
        int iCat = itCat - categories_to_fit.begin();
        double v_sideband = input_histos[*itSideband][*itCat].Integral(fit_bins_start[i], fit_bins_end[i]);
        double v_scale = scale[i*nFitCategories+iCat];
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


      double dchi2 = pow(sum-v_data, 2);
      chi2+= dchi2;
      //cout<< i<<", "<<iSideband<<" "<<sum<<", "<<v_data<<"--chi2: "<<dchi2<<"...."<<endl;

      //chi2_shape
      for( int ibin = fit_bins_start[i]; ibin <= fit_bins_end[i]; ibin++ )
      {
        double sum_reg = 0;
        for( vector<int>::iterator itCat = categories_to_fit.begin(); itCat!=categories_to_fit.end();++itCat )
          sum_reg+= scale[i*nFitCategories+(itCat-categories_to_fit.begin())]*input_histos[*itSideband][*itCat].Integral(ibin,ibin);
          //sum_reg+= pow(scale[i*nFitCategories+(itCat-categories_to_fit.begin())], 2);
        for( vector<int>::iterator itCat = categories_stay_const.begin(); itCat!=categories_stay_const.end();++itCat )
          sum_reg+= input_histos[*itSideband][*itCat].Integral(ibin,ibin);
        double data_reg = input_histos[*itSideband][kData].Integral(ibin,ibin);
        chi2_shape+=pow(sum_reg-data_reg,2);
      }

    }
  }
  
  return chi2+chi2_shape;
}


double result = 0;
void minuitFunction( int &nDim, double* gout, double &result, double par[], int flg)
{
  result = chi2Function(par);
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
  vector<TH1D*> cv_background_scale_histos = GetBackgroundScales( cv_scales_histos, kQELikeNot);

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
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kQELikeNot);
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
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kQELikeNot);
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


void getAllUniverseScales1D( vector<MnvH1D**> tracks,vector<MnvH1D*> &track_scale )
{
  //vector<MnvH1D*> ret(nHistos);

  //CV
  for( int i = 0; i<nHistos; i++ )
  {
    for(int j = 0; j<Samples.size();++j) input_histos[j][i] = *(new TH1D( tracks[j][i]->GetCVHistoWithStatError() ));
  }

  vector<TH1D*> cv_scales_histos = GetSingleTypeFitScalesAndErrors();// scales for pizero, piplus and multipi
  vector<TH1D*> cv_background_scale_histos = GetBackgroundScales( cv_scales_histos, kQELikeNot);

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
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, histo_name.c_str(), vertNames[k].c_str(),i_uni)));
          cout<<Form("Vert_%d_%s_%d has %d events", iH, names[*ittype].c_str(), i_uni, input_histos[iH][*ittype].Integral() )<<endl;
        }

      }
      //now, fit this universe
      //0:signal, 1:blob, 2:michel,3:micblob
      vector<TH1D*> uni_scales_histos = GetSingleTypeFitScalesAndErrors();
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kQELikeNot);
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
      vector<TH1D*> uni_background_scales_histos = GetBackgroundScales( uni_scales_histos, kQELikeNot);
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
  for( int i = 0; i<tracks.size(); i++ )
  {
    for(int j = 0; j < nHistos; ++j ) input_histos[i][j] = *(new TH1D( tracks[i][j]->GetCVHistoWithStatError() ));
    cout<<"CV kData "<<input_histos[i][kData].Integral()<<endl;
    cout<<"CV kQELike_QE_OTH "<<input_histos[i][kQELike_QE_OTH].Integral()<<endl;
    cout<<"CV kQELike_RES "<<input_histos[i][kQELike_RES].Integral()<<endl;
    cout<<"CV kQELike_2p2h "<<input_histos[i][kQELike_2p2h].Integral()<<endl;
    cout<<"CV kQELike_OTH "<<input_histos[i][kQELike_OTH].Integral()<<endl;
    //cout<<"Signal Histo: "<<i<<" has "<<input_histos["Signal"][i].Integral()<<" events"<<endl;
    //input_histos["Signal"][i].Draw("histe");
    //c.Print( Form("./test/Signal_%02d_%s.png", i, names[i].c_str() ) );
    //input_histos["BlobSideBand"][i].Draw("histe");
    //c.Print( Form("./test/BlobSideBand_%02d_%s.png", i, names[i].c_str() ) );
    //input_histos["MichelSideBand"][i].Draw("histe");
    //c.Print( Form("./test/MichelSideBand_%02d_%s.png", i, names[i].c_str() ) );
    //input_histos["MicBlobSideBand"][i].Draw("histe");
    //c.Print( Form("./test/MicBlobSideBand_%02d_%s.png", i, names[i].c_str() ) );
  }

  vector<TH1D*> cv_scales_histos = GetSingleTypeFitScalesAndErrors();// scales for pizero, piplus and multipi
  for( int t = 0; t< nTypes ; t++) 
  {
    setScaleAndError( track_scale[t], cv_scales_histos[t] );
    //cv_scales_histos[t]->Draw("histe");
    //c.Print( Form("./test_%02d.png", t ));
  }


  if(!doErrorBands) return;
  //Vert Error
  vector<string> vertNames = tracks[0][kMC]->GetVertErrorBandNames();
  for( unsigned int k= 0; k < vertNames.size(); ++k )
  {
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
        input_histos[iH][kQELike] = *((TH1D*) tracks[iH][kQELike]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, "data", vertNames[k].c_str(),i_uni)));
        //input_histos[iH][kData].Draw("hist");
        cout<<Form("vert_2track_angle_%02d_%s_%s_%d integral: %d", iH, "data", vertNames[k].c_str(), i_uni, (int) input_histos[iH][kData].Integral() )<<endl;
        cout<<Form("vert_2track_angle_%02d_%s_%s_%d integral: %d", iH, "qelike", vertNames[k].c_str(), i_uni, (int) input_histos[iH][kQELike].Integral() )<<endl;
        for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
        { 
          string histo_name = names[*ittype];
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetVertErrorBand( vertNames[k] )->GetHist(i_uni)->Clone(Form("vert_2track_angle_%02d_%s_%s_%d", iH, names[*ittype].c_str(), vertNames[k].c_str(),i_uni)));
          //input_histos[iH][*ittype].Draw("histlsame");
          cout<<Form("vert_2track_angle_%02d_%s_%s_%d integral: %d", iH, names[*ittype].c_str(), vertNames[k].c_str(), i_uni, (int) input_histos[iH][*ittype].Integral() )<<endl;
        }
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
  for( unsigned int k= 0; k < latNames.size(); ++k )
  {
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
        input_histos[iH][kQELike] = *((TH1D*) tracks[iH][kQELike]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, "data", latNames[k].c_str(),i_uni)));
        cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, "data", latNames[k].c_str(), i_uni, (int) input_histos[iH][kData].Integral() )<<endl;
        cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, "qelike", latNames[k].c_str(), i_uni, (int) input_histos[iH][kQELike].Integral() )<<endl;
        for( vector<int>::iterator ittype = all_universe_histos.begin(); ittype!=all_universe_histos.end();++ittype )
        { 
          string histo_name = names[*ittype];
          input_histos[iH][*ittype] = *((TH1D*) tracks[iH][*ittype]->GetLatErrorBand( latNames[k] )->GetHist(i_uni)->Clone(Form("lat_2track_angle_%02d_%s_%s_%d", iH, names[*ittype].c_str(), latNames[k].c_str(),i_uni)));
          cout<<Form("lat_2track_angle_%02d_%s_%s_%d integral: %d", iH, names[*ittype].c_str(), latNames[k].c_str(), i_uni, (int) input_histos[iH][*ittype].Integral() )<<endl;
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
  TFile *fin;
  TFile *fweight;
  

  cout << "Loading file: " << filenames[1] << endl;
  fin = new TFile( filenames[1].c_str() , "READ" );
  fweight = new TFile( filenames[2].c_str(), "READ");
  if (fin->IsZombie() || fin->GetListOfKeys()->IsEmpty())
  {
    Error("SideBandFit","Could not get histogram ROOT file or it was empty.");
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
  MnvH1D* h_bck_weights = (MnvH1D*) fweight->Get("h_weights_q2qe_qelikenot_Signal");
  
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
  CCQENuUtils     *utils        = new CCQENuUtils( false, fluxHistoExists );

  //--------------------------------------------------
  // Get PT, PZ binning for the weights histograms 
  //--------------------------------------------------
  cout << "Done with setting up general functions" << endl;
  const int n_fitvarbins_2track = 1; // only 1 version of the binning exists

  
  //---------------------------------------------------------------
  // set which bins to use
  //--------------------------------------------------------------
  //axis_binning ybins = muonPtbins;
  axis_binning ybins = Q2bins;
  if (fit_xvar == "h_ptmu" ) ybins = muonPtbins;

  //---------------------------------------------------------------
  // Create vector of pointers to pointers
  // Each element of vector is a pointer to an array of MnvH1D 
  //---------------------------------------------------------------
  // vector < vector < MnvH1D*[nHistos]> > 
  vector<vector<MnvH1D*>>  h_fit_varbins_samples;
  for( int i = 0; i<n_samples; i++ ) 
  {
    h_fit_varbins_samples.push_back( vector<MnvH1D*>(nHistos) );
  }

  //signal, blob, michel, micblob
  //for( int i = 0; i<n_samples; i++ )
  //{
  //  vector<MnvH1D**> vec( n_fitvarbins_2track );
  //  initializeVector( vec ); 
  //  h_fit_varbins_samples.push_back( vec );
  //}

  double pot_data = plotutils->getPOTData( fin );
  double pot_mc = plotutils->getPOTMC( fin );
  double pot_scale = plotutils->getPOTNormFactor( fin );
  
  if(pot_scale==0) pot_scale = 1.0;
  cout<< "---------------------------------------"  << endl;
  cout<< "POT INFORMATION:"  << endl;
  cout<< "POT Data = " << pot_data << endl;
  cout<< "POT MC   = " << pot_mc << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;

  for( int j = 0; j < n_samples; j++ )
  {
    plotutils->bookHistos( fin, &h_fit_varbins_samples[j][0], fit_xvars[j] );      
    cout<<"sample name: "<<h_fit_varbins_samples[j][kData]->GetName()<<endl;
    cout<<"sample bins: "<<h_fit_varbins_samples[j][kData]->GetNbinsX()<<endl;
    plotutils->scaleMCHistos( &h_fit_varbins_samples[j][0], pot_scale);

    cout<< "Subtracting Data"<<endl;
    MnvH1D* h_bck = (MnvH1D*) h_fit_varbins_samples[j][kQELikeNot]->Clone("h_bck");
    h_bck->Multiply( h_bck, h_bck_weights );
    cout<<j<<" has "<<h_fit_varbins_samples[j][kData]->Integral()<<" original data"<<endl;
    h_fit_varbins_samples[j][kData]->AddMissingErrorBandsAndFillWithCV( *h_bck_weights );
    h_fit_varbins_samples[j][kData]->Add( h_bck, -1 );
  
  }
  for( int i = 0; i<n_samples;i++)
  {
    cout<<i<<" has "<<h_fit_varbins_samples[i][kData]->Integral()<<" subtracted data"<<endl;
    cout<<i<<" has "<<h_fit_varbins_samples[i][kQELike]->Integral()<<" qelike event"<<endl;
    cout<<i<<" has "<<h_fit_varbins_samples[i][kQELike_QE_OTH]->Integral()<<" qe_oth event"<<endl;
    cout<<i<<" has "<<h_fit_varbins_samples[i][kQELike_RES]->Integral()<<" res event"<<endl;
    cout<<i<<" has "<<h_fit_varbins_samples[i][kQELike_DIS]->Integral()<<" dis event"<<endl;
    cout<<i<<" has "<<h_fit_varbins_samples[i][kQELike_2p2h]->Integral()<<" 2p2h event"<<endl;
    //for( vector<int>::iterator itType = all_universe_histos.begin(); itType != all_universe_histos.end(); itType++) PrintMnvH1D( h_fit_varbins_samples[i][*itType] );
  }


  cout << "Done with Tejin's track" << endl;
  fin->cd();
  fin->Close();
  cout << "Done booking histograms" << endl;


  cout<< " Done with tracks "<<endl;
  //vector<MnvH2D*> h_vect_dalphat_scales;
  //vector<MnvH1D*> h_vect_yvar_scales;
  //vector<MnvH2D*> h_vect_dthetaP_scales;
  //vector<MnvH2D*> h_vect_dthetaR_scales;
  //vector<MnvH2D*> h_vect_dalphat_scales;
  //vector<MnvH2D*> h_vect_pn_scales;
  //vector<MnvH2D*> h_vect_dpt_scales;
  //vector<MnvH2D*> h_vect_dthetaPq2qe_scales;
  //vector<MnvH2D*> h_vect_dthetaRq2qe_scales;
  //for( int i = 0; i<n_samples; i++ ) 
  //{
  //  string hname;
  //  axis_binning xbins;

  //  hname = Form("h_weights_q2qe_qelikenot_%02d", Samples[i] );
  //  h_vect_yvar_scales.push_back( new MnvH1D( hname.c_str(), hname.c_str(), ybins.nbins, &(ybins.bin_edges[0]) ) ) ;

  //}

  ////vector<MnvH2D*> h_dalphat_scales, h_dthetaR_scales, h_dthetaP_scales;
  //getAllUniverseScales1D(h_fit_varbins_samples, h_vect_yvar_scales);
  
  //ExpandHistos( h_vect_yvar_scales, h_vect_dthetaP_scales, 1);
  cout << "*********************************" <<endl;
  cout << "GetSingleType" << endl;
  cout << "*********************************" << endl;
  vector<MnvH1D*> hs_vect_yvar_scales;
  vector<MnvH2D*> hs_vect_dthetaP_scales;
  vector<MnvH2D*> hs_vect_dthetaR_scales;
  vector<MnvH2D*> hs_vect_dalphat_scales;
  vector<MnvH2D*> hs_vect_pn_scales;
  vector<MnvH2D*> hs_vect_dpt_scales;
  //vector<MnvH2D*> hs_vect_dthetaPq2qe_scales;
  //vector<MnvH2D*> hs_vect_dthetaRq2qe_scales;
  int nTypes = bg_categories_to_scale.size();
  for( int i = 0; i<nTypes; i++ ) 
  {
    string typeName = bg_categories_to_scale_names[i];
    string hname;
    axis_binning xbins;

    hname = Form("hs_weights_yvar_bgType_%s", typeName.c_str() );
    hs_vect_yvar_scales.push_back( new MnvH1D( hname.c_str(), hname.c_str(), ybins.nbins, &(ybins.bin_edges[0]) ) ) ;
    
    //xbins = dthetaReactbins;
    //hname = Form( "hs_weights_dthetaR_yvarbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dthetaR_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    //xbins = dthetaPerpbins;
    //hname = Form( "hs_weights_dthetaP_yvarbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dthetaP_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );
    //
    //xbins = dalphatbins;
    //hname = Form( "hs_weights_dalphat_yvarbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dalphat_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    //xbins = pnbins;
    //hname = Form( "hs_weights_pn_yvarbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_pn_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    //xbins = dptbins;
    //hname = Form( "hs_weights_dpt_yvarbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dpt_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), xbins.nbins, &(xbins.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0]) ) );

    //axis_binning hdlx = dthetaPerpbins; 
    //axis_binning hdlz = Q2bins;
    //axis_binning global_x;
    //GetGlobal_x( hdlx, ybins, hdlz, global_x );
    //hname = Form( "hs_weights_dthetaPq2qeptbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dthetaPq2qe_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), global_x.nbins, &(global_x.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] ) ) );


    //hdlx = dthetaReactbins; 
    //hdlz = Q2bins;
    //GetGlobal_x( hdlx, ybins, hdlz, global_x );
    //hname = Form( "hs_weights_dthetaRq2qeptbins_bgType_%s", typeName.c_str() ) ;
    //hs_vect_dthetaRq2qe_scales.push_back( new MnvH2D( hname.c_str(), hname.c_str(), global_x.nbins, &(global_x.bin_edges[0]), ybins.nbins, &( ybins.bin_edges[0] )) ) ;


  }

  //vector<MnvH2D*> hs_dalphat_scales, hs_dthetaR_scales, hs_dthetaP_scales;

  getAllUniverseScalesSingleType1D(h_fit_varbins_samples, hs_vect_yvar_scales);

  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaP_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaR_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dalphat_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_pn_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dpt_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaPyvar_scales, 1);
  //ExpandHistos( hs_vect_yvar_scales, hs_vect_dthetaRyvar_scales, 1);



  cout << "*********************************" << endl;
  cout << "Done Saving Ouput" << endl;
  cout << "*********************************" << endl;
  TFile *tmp = new TFile(filenames[3].c_str(),"RECREATE"); //micblob
  tmp->cd();
  for(int i=0;i<n_samples;i++){
    cout<<"Saving... "<<Samples[i]<<endl;
    //h_vect_yvar_scales[i]->Write();
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
    //hs_vect_dthetaP_scales[i]->Write();
    //hs_vect_dthetaR_scales[i]->Write();
    //hs_vect_dalphat_scales[i]->Write();
    //hs_vect_pn_scales[i]->Write();
    //hs_vect_dpt_scales[i]->Write();
    //hs_vect_dthetaPq2qe_scales[i]->Write();
    //hs_vect_dthetaRq2qe_scales[i]->Write();
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
  par.push_back( Form("%s/ana/rootfiles/MuonEventSignal_2track.root",getenv("CCQENUROOT") ) );
  par.push_back( Form("%s/ana/rootfiles/SideBandFitWeights.root",getenv("CCQENUROOT") ) ); //qelikenot subtraction
  par.push_back( Form("%s/ana/rootfiles/SignalFitWeightsAngle.root",getenv("CCQENUROOT") ) ); //output file


  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
  
  return SideBandFit(par);
}

