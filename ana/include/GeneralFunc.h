#ifndef __GEN_FUNC_H
#define __GEN_FUNC_H
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"


using namespace CCQENU_ANA;

double minQ2 = 0.05;
double maxQ2 = 4;

MnvH2D* ProjectionHDL( vector<MnvH2D*> &XYinZ, string name,HyperDimLinearizer* hdl, axis_binning &xbin, axis_binning& ybin, axis_binning &zbin,  string dims="xy", int projmin=1, int projmax=-1);//inclusive proj
MnvH2D* ProjectionXZ( vector<MnvH2D*> &hists, string name,  int ymin = 1, int ymax = -1 );

void GetGlobal_x( axis_binning &xbins, axis_binning&ybins, axis_binning&zbins, axis_binning &global_x, int mode = 0 );
HyperDimLinearizer* GetHDL( axis_binning &xbins, axis_binning &ybins, axis_binning &zbins );
HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x, bool declareErrBand = true);

void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuTruth* evt, double wgt , bool runSys = true, bool pass_cv = true);
void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt , bool runSys = true, bool runLat = true, bool pass_cv = true);
void FillHDL2D( std::string name, CCQENuUtils *utils, string var_x, string var_y, string var_z, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt , bool runSys = true, bool runLat = true, bool pass_cv = true);

void FillHDL2DTruth( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuTruth* evt, double wgt , bool runSys = true);

void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, TH2D *h2d, TH3D *h3d );
void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, MnvH2D *h2d, MnvH3D *h3d );
void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, MnvH2D **h2d, MnvH3D **h3d );

void ExpandHisto( TH1D* h1, TH2D* h2, int xaxis=0, bool err=true );// 0:x->expand into y, 1:y->expand into x
void ExpandHisto1D( TH1D* h1, TH1D* h2 );// 0:x->expand into y, 1:y->expand into x
void ExpandHisto( MnvH1D* h1, MnvH2D*h2, int xaxis = 0, bool err=true );//MnvH2D doesn't have any errorbands defined
void ExpandHisto1D( MnvH1D* h1, MnvH1D*h2);//MnvH2D doesn't have any errorbands defined
void ExpandHistos( vector<MnvH1D*> h1vec, vector<MnvH2D*> h2vec, int xaxis = 0, bool err=true );//MnvH2D doesn't have any errorbands defined
void ExpandHistos1D( vector<MnvH1D*> h1vec, vector<MnvH1D*> h2vec );//MnvH2D doesn't have any errorbands defined

void SetHists1bin( TH1D* h1, TH2D* h2, int ibin, bool fillx = true );//set to the x bin or ybin
void SetHists1bin( MnvH1D* h1, MnvH2D* h2, int ibin, bool fillx = true );
void SetHists( vector<MnvH1D*> h1vec, MnvH2D* h2, bool fillx = true );

XYZVector Convert( TVector3 v ){ return XYZVector( v.X(), v.Y(),v.Z()); };
XYZTVector Convert( TLorentzVector v ){ return XYZTVector( v.X(), v.Y(),v.Z(), v.T()); };
TVector3 Convert( XYZVector  v ) { return TVector3(v.X(), v.Y(),v.Z()); };
TLorentzVector Convert( XYZTVector v) { return TLorentzVector(v.X(), v.Y(),v.Z(), v.T() ); };

template<class A,class B>
vector<A> Convert( vector<B> &vec )
{
  vector<A> ret;
  for( auto v: vec ) ret.push_back( Convert(v) );
  return ret;
}
template vector<XYZTVector> Convert( vector<TLorentzVector> &vec );
template vector<XYZVector> Convert( vector<TVector3> &vec );
template vector<TLorentzVector> Convert( vector<XYZTVector> &vec );
template vector<TVector3>  Convert( vector<XYZVector> &vec );


//Projections
//MnvH1D* ProjectionY

template< class T > TLorentzVector GetParticleMomentum( T* event, int i );
template< class T > std::vector< TLorentzVector > GetAllParticlesMomenta( T*event, int pdg = 2212, double Ethreshold = 0 );


string GetQuadrantCode( XYZVector &beam, XYZVector &muon, XYZVector &proton, XYZVector &neutron );
int getRegions(double dthetaP, double dthetaR);

//----- HyderDim Related ----------
void GetGlobal_x( axis_binning &xbins, axis_binning&ybins, axis_binning&zbins, axis_binning &global_x, int mode/* = 0*/ )
{
  std::vector< std::vector<double> > bins3D({xbins.bin_edges, ybins.bin_edges, zbins.bin_edges});
  HyperDimLinearizer hdl(bins3D,mode);
  int n_bins = (xbins.bin_edges.size()+1)*(zbins.bin_edges.size()+1 );
  vector<double> bins;
  for( int i = 0; i<n_bins; ++i ) bins.push_back( i );
  global_x.bin_edges = bins;
  global_x.nbins = bins.size()-1;
  global_x.min = bins.front();
  global_x.max = bins.back();
}




MnvH2D* ProjectionHDL( vector<MnvH2D*> &XYinZ, string name,HyperDimLinearizer* hdl, axis_binning &xbin, axis_binning& ybin, axis_binning &zbin,  string dims/*="xy"*/, int projmin/*=1*/, int projmax/*=-1*/)//inclusive proj
{
  if(!( dims=="xy" || dims == "zy" || dims == "xz") ) throw "dims must be xy, zy or xz";

  MnvH2D* ret = nullptr;
  cout<<ret<<endl;
  if(dims == "xy")
  {
     ret= ProjectionXZ( XYinZ, name, projmin, projmax );
  }
  else
  {
     if(dims=="xz") ret=new MnvH2D( name.c_str(), name.c_str(), xbin.nbins, &(xbin.bin_edges[0]),
                   zbin.nbins, &(zbin.bin_edges[0]));
     if(dims=="zy") ret=new MnvH2D( name.c_str(), name.c_str(), zbin.nbins, &(zbin.bin_edges[0]),
                   ybin.nbins, &(ybin.bin_edges[0]));
     cout<<"add missing errorbands..";
     ret->AddMissingErrorBandsAndFillWithCV( (*XYinZ[0]) );
     cout<<"done"<<endl;
     bool projectX = (dims=="xz");
     for(unsigned int iz = 0; iz< XYinZ.size();iz++)
     {
       MnvH1D* XorY =(MnvH1D*) XYinZ[iz]->Projection(name.c_str(), projectX,projmin, projmax,"" );
       if(projectX)//xz
       {//fill into hist y axis, since histy is z 
           SetHists1bin( XorY, ret, iz+1, false); 
       } else //zy
       {//fill into hist x axis, since histx is z 
           SetHists1bin( XorY, ret, iz+1, true) ;
       }
       delete XorY;
     }
  }
  //for( vector<MnvH2D*>::iterator it = XYinZ.begin(); it!= XYinZ.end();++it) delete *it;
  cout<<ret<<endl;
  return ret;
}



MnvH2D* ProjectionXZ( vector<MnvH2D*> &hists, string name,  int ymin /*= 1*/, int ymax /*= -1*/ )
{
  MnvH2D* ret = nullptr;
  int nbins = hists.size();
  if( ymin > nbins ) return ret;


  ret = (MnvH2D*) hists[ymin-1]->Clone( name.c_str() );
  int maxbin = nbins;
  if (ymax > ymin ) maxbin = ymax;
  for( int ibin = ymin-1; ibin < maxbin; ++ibin ) ret->Add( hists[ibin] );
  return ret;
}


HyperDimLinearizer* GetHDL( axis_binning &xbins, axis_binning &ybins, axis_binning &zbins )
{
  std::vector< std::vector<double> > bins3D({xbins.bin_edges, ybins.bin_edges, zbins.bin_edges});
  HyperDimLinearizer *hdl = new HyperDimLinearizer(bins3D,0);
  return hdl;
}


HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x, bool declareErrBand/* = true*/)
{
  //std::vector< std::vector<double> > bins3D({xbins.bin_edges, ybins.bin_edges, zbins.bin_edges});

  HyperDimLinearizer *hdl = GetHDL( xbins, ybins, zbins );

  int n_bins = (xbins.bin_edges.size()+1)*(zbins.bin_edges.size()+1 );
  vector<double> bins;
  for( int i = 0; i<n_bins; ++i ) bins.push_back( i );

  global_x.bin_edges = bins;
  global_x.nbins = bins.size()-1;
  global_x.min = bins.front();
  global_x.max = bins.back();
  utils->bookHistos( h, name.c_str(), title.c_str(), global_x, ybins );
  if(declareErrBand)
  {
    utils->addLatErrorBands( h );
    utils->addVertErrorBands( h );
  }
  utils->HDLMap[ name.c_str() ] = hdl;
  return hdl;
}

void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuTruth* evt, double wgt , bool runSys /*= true*/, bool pass_cv)
{
  //std::vector<double> values({x,y,z});
  HyperDimLinearizer* hdl = utils->HDLMap[ name.c_str() ];
  double globalx = hdl->GetBin( std::vector<double>({x,y,z}) ).first+0.0001;
  utils->fillHistosV3( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );

  if( runSys && !isData)
  {
	  utils->fillVertErrorBands( utils->histos2D[ name.c_str() ],globalx, y,evt);
  }
}
void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt , bool runSys /*= true*/, bool runLat, bool pass_cv )
{
  //std::vector<double> values({x,y,z});
    HyperDimLinearizer* hdl = utils->HDLMap[ name.c_str() ];
    double globalx = hdl->GetBin( std::vector<double>({x,y,z}) ).first+0.0001;

  if(pass_cv)
  {
    utils->fillHistosV3( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );
    //utils->fillBlobHist( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );
 
    if( runSys && !isData) utils->fillVertErrorBands( utils->histos2D[ name.c_str() ],globalx, y,evt);
  }
  if( runSys && !isData && runLat)
  {
	  utils->fillLatErrorBands( utils->histos2D[ name.c_str() ],"dthetaPdthetaR","q2",globalx, y,evt,pass_cv);
  }
}

void FillHDL2D( std::string name, CCQENuUtils *utils, string var_x, string var_y, string var_z,double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt , bool runSys /*= true*/, bool runLat, bool pass_cv )
{
  HyperDimLinearizer* hdl = utils->HDLMap[ name.c_str() ];
  double globalx = hdl->GetBin( std::vector<double>({x,y,z}) ).first+0.0001;
  if(pass_cv)
  {
    //std::vector<double> values({x,y,z});
    utils->fillHistosV3( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );
    //utils->fillBlobHist( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );

    if( runSys && !isData) utils->fillVertErrorBands( utils->histos2D[ name.c_str() ],globalx, y,evt);
  }
  if( runSys && runLat && !isData ) utils->fillLatErrorBands( utils->histos2D[ name.c_str() ], var_x, var_y, var_z, x, y, z, evt, hdl, pass_cv, true  );
}


void FillHDL2DTruth( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuTruth* evt, double wgt , bool runSys /*= true*/)
{
  //std::vector<double> values({x,y,z});
  HyperDimLinearizer* hdl = utils->HDLMap[ name.c_str() ];
  double globalx = hdl->GetBin( std::vector<double>({x,y,z}) ).first+0.0001;
  utils->fillHistosV3( utils->histos2D[ name.c_str() ], globalx, y, isData, evt, wgt );

  if( runSys && !isData)
  {
	  utils->fillVertErrorBands( utils->histos2D[ name.c_str() ],globalx, y,evt);
  }
}



void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, TH2D* h2d, TH3D *h3d )
{
  //cout<<"start Convert TH2D/TH3D"<<endl;
  for(int i = 0;i<h2d->GetNbinsX(); ++i )
  {
    int global_x = i+1;
    vector<int> xzbins = hdl->GetValues( global_x );
    cout<<xzbins.size()<<endl;
    int xbin = xzbins[0]-1;
    int zbin = xzbins[2];
    cout<<Form("Bin Conversion: %i-> %i %i", global_x+1, xbin,zbin)<<endl;
    
    //int nbinsy = h2d->GetNbinsY();
    //cout<<nbinsy<<endl;
    for( int ii = 0; ii< h2d->GetNbinsY();ii++ )
    {
      int ybin = ii+1;
      //cout<<"xyzbins: "<<xbin<<", "<<ybin<<", "<<zbin<<endl;
      double v = h2d->GetBinContent( global_x, ybin );
      double e = h2d->GetBinError( global_x, ybin );
      //cout<<v<<", "<<e<<endl;
      h3d->SetBinContent( xbin,ybin,zbin, v );
      h3d->SetBinError( xbin,ybin,zbin, e );
      //cout<<"error set"<<endl;
    }
  }
  //cout<<Form("conversion-- %.2f %.2f", h2d->Integral(), h3d->Integral() )<<endl;
}


void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, MnvH2D* h2d, MnvH3D *h3d )
{
  //cv
  ConvertHDL2DTo3D( hdl, (TH2D*) h2d, (TH3D*) h3d );

  /*
  //vert err band
  vector<string> vertNames = h2d->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    //cout<<vertNames[k]<<",";
    int nUniverses = h2d->GetVertErrorBand( vertNames[k] )->GetNHists();
    h3d->AddVertErrorBand(vertNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* ph1 = h2d->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      TH3D* ph2 = h3d->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      ConvertHDL2DTo3D(hdl, ph1, ph2);
    }
  }
  //lat err band
  vector<string> latNames = h2d->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    //cout<<latNames[k]<<",";
    int nUniverses = h2d->GetLatErrorBand( latNames[k] )->GetNHists();
    h3d->AddLatErrorBand(latNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* ph1 = h2d->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      TH3D* ph2 = h3d->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      ConvertHDL2DTo3D(hdl, ph1, ph2);
    }
  }
  */
}

void ConvertHDL2DTo3D( HyperDimLinearizer* hdl, MnvH2D** h2d, MnvH3D **h3d )
{
  for(unsigned int i = 0; i< nHistos;i++ ) ConvertHDL2DTo3D( hdl, h2d[i], h3d[i] ) ;
}


//expanding TH1D/MnvH1D to 2D:
void ExpandHisto( TH1D* h1, TH2D* h2, int xaxis/*=0*/ , bool includeError)// 0:x->expand into y, 1:y->expand into x
{
  int nX = h2->GetNbinsX();
  int nY = h2->GetNbinsY();
  for( int ix = 0; ix <= nX+1; ix ++ ) {
    for( int iy = 0; iy <= nY+1;iy ++ ){
      if(xaxis == 0 )
      {
        h2->SetBinContent(ix,iy, h1->GetBinContent(ix) );
        h2->SetBinError(ix,iy, (includeError)? h1->GetBinError(ix):0 );
      }
      else
      {
        h2->SetBinContent(ix,iy, h1->GetBinContent(iy) );
        h2->SetBinError(ix,iy, (includeError)? h1->GetBinError(iy):0 );
      }
    }
  }
}
void ExpandHisto1D( TH1D* h1, TH1D* h2 )// 0:x->expand into y, 1:y->expand into x
{
  int nX = h2->GetNbinsX();
  for( int ix = 0; ix <= nX+1; ix ++ ) {

    //double ic = h1->GetBinContent(ix);
    //double ie = h1->GetBinError(ix);
    //cout<<"bin content/error: "<<ic<<", "<<ie<<endl;
    h2->SetBinContent(ix, h1->GetBinContent(ix) );
    h2->SetBinError(ix, h1->GetBinError(ix) );
    //cout<<"h2 bin content/error: "<<h2->GetBinContent(ix)<<", "<<h2->GetBinError(ix)<<endl;
  }
}

void ExpandHisto( MnvH1D* h1, MnvH2D*h2, int xaxis/* = 0*/, bool includeError )//MnvH2D doesn't have any errorbands defined
{
  h2->ClearAllErrorBands();
  h2->Reset();
  //AddLat and VertErrBands
  //AddErrBands( h1, h2 );
  //CV:
  //cout<<"Expanding "<<h2->GetName()<<endl;
  ExpandHisto( (TH1D*) h1, (TH2D*)h2, xaxis, includeError );
  //LatErrBand
  //cout<<"Expanding VertErrBand"<<endl;
  vector<string> vertNames = h1->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    //cout<<vertNames[k]<<",";
    int nUniverses = h1->GetVertErrorBand( vertNames[k] )->GetNHists();
    h2->AddVertErrorBand(vertNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH1D* ph1 = h1->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      TH2D* ph2 = h2->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      ExpandHisto( ph1, ph2, xaxis, includeError );
    }
  }

  //cout<<"Expanding LatErrBand"<<endl;
  vector<string> latNames = h1->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    int nUniverses  = h1->GetLatErrorBand( latNames[k] )->GetNHists();
    h2->AddLatErrorBand(latNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) 
    {
      TH1D* ph1 = h1->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      TH2D* ph2 = h2->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      ExpandHisto( ph1, ph2, xaxis, includeError );
    }
  }
}

void ExpandHisto1D( MnvH1D* h1, MnvH1D*h2)//MnvH2D doesn't have any errorbands defined
{
  //AddLat and VertErrBands
  //AddErrBands( h1, h2 );
  //CV:
  ExpandHisto1D( (TH1D*) h1, (TH1D*)h2 );
  //VertErrBand
  vector<string> vertNames = h1->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    int nUniverses = h1->GetVertErrorBand( vertNames[k] )->GetNHists();
    h2->AddVertErrorBand(vertNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++)
    {
      TH1D* ph2 = h2->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      TH1D* ph1 = h1->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      ExpandHisto1D( ph1,ph2);
    }
      
  }

  //LatErrBand
  vector<string> latNames = h1->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    int nUniverses  = h1->GetLatErrorBand( latNames[k] )->GetNHists();
    h2->AddLatErrorBand(latNames[k],nUniverses);
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) 
    {
      TH1D* ph2 = h2->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      TH1D* ph1 = h1->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      ExpandHisto1D(ph1, ph2);
    }
  }
}


void ExpandHistos( vector<MnvH1D*> h1vec, vector<MnvH2D*> h2vec, int xaxis/* = 0*/, bool includeError)//MnvH2D doesn't have any errorbands defined
{
  for( unsigned int i = 0; i< h1vec.size() ; i++ ) ExpandHisto( h1vec[i], h2vec[i], xaxis, includeError );
}
void ExpandHistos1D( vector<MnvH1D*> h1vec, vector<MnvH1D*> h2vec )//MnvH2D doesn't have any errorbands defined
{
  for( unsigned int i = 0; i< h1vec.size() ; i++ ) ExpandHisto1D( h1vec[i], h2vec[i]);
}


void SetHists1bin( TH1D* h1, TH2D* h2, int ibin, bool fillx /*= true*/ )//set to the x bin or ybin
{

  if(fillx)
  {
    for( int j = 0; j<h1->GetNbinsX();j++) 
    {
      int jbin = j+1;
      h2->SetBinContent( ibin,jbin, h1->GetBinContent(jbin) );
      h2->SetBinError( ibin,jbin, h1->GetBinError(jbin) );
    }
  }
  else
  {
    for( int j = 0; j<h1->GetNbinsX();j++) 
    {
      int jbin = j+1;
      h2->SetBinContent( jbin,ibin, h1->GetBinContent(jbin) );
      h2->SetBinError( jbin,ibin, h1->GetBinError(jbin) );
    }
  }
}

void SetHists1bin( MnvH1D* h1, MnvH2D* h2, int ibin, bool fillx /*(= true */)
{
  //AddLat and VertErrBands
  //AddErrBands( h1, h2 );
  //CV:
  //cout<<"Expanding "<<h2->GetName()<<endl;
  SetHists1bin( (TH1D*) h1, (TH2D*)h2, ibin, fillx );
  //VertErrBand
  vector<string> vertNames = h1->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    int nUniverses = h1->GetVertErrorBand( vertNames[k] )->GetNHists();
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++)
    {
      TH2D* ph2 = h2->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      TH1D* ph1 = h1->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      SetHists1bin(ph1,ph2, ibin, fillx);
    }
  }

  //LatErrBand
  vector<string> latNames = h1->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    int nUniverses  = h1->GetLatErrorBand( latNames[k] )->GetNHists();
    for (int iUniv = 0 ;iUniv< nUniverses; iUniv++) 
    {
      TH2D* ph2 = h2->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      TH1D* ph1 = h1->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      SetHists1bin(ph1,ph2, ibin, fillx);
    }
  }
}

void SetHists( vector<MnvH1D*> h1vec, MnvH2D* h2, bool fillx /*= true*/ )
{
  for(unsigned int i = 0;i< h1vec.size();i++)
  {
    int ibin = i+1;
    SetHists1bin( h1vec[i], h2, ibin, fillx );
  }
}

//search particle
template< class T >
TLorentzVector GetParticleMomentum( T* event, int i )
{
  TLorentzVector vec(-999,-999,-999,-999);
  if ( i > event->mc_nFSPart ) return vec;
  vec.SetPxPyPzE( event->mc_FSPartPx[i], event->mc_FSPartPy[i], event->mc_FSPartPz[i], event->mc_FSPartE[i]); 
  return vec;
}

//TLorentzVector GetParticleMomentum( CCQENuEvent *event, int i );
//TLorentzVector GetParticleMomentum( CCQENuTruth *event, int i );

template< class T >
std::vector< TLorentzVector > GetAllParticlesMomenta( T*event, int pdg/* = 2212*/, double Ethreshold/* = 0*/ )
{
  std::vector< TLorentzVector> ret;
  for( int i = 0; i< event->mc_nFSPart; ++i )
  {
    if( pdg != event->mc_FSPartPDG[i] ) continue;
    if( event->mc_FSPartE[i] < Ethreshold ) continue;
    ret.push_back( GetParticleMomentum( event, i )) ;
  }
  std::sort( ret.begin(), ret.end(), [](TLorentzVector a, TLorentzVector b){ return a.E() > b.E(); } );
  return ret;
}

//std::vector< TLorentzVector > GetAllParticlesMomenta( CCQENuEvent*event, int pdg = 2212, double Ethreshold = 0 );
//std::vector< TLorentzVector > GetAllParticlesMomenta( CCQENuTruth*event, int pdg = 2212, double Ethreshold = 0 );

string GetQuadrantCode( XYZVector &beam, XYZVector &muon, XYZVector &proton, XYZVector &neutron )
{
  XYZVector xcoord = beam.Cross(muon).Unit();
  string p = ( xcoord.Dot(proton) <= 0 )? "l":"r";
  string n = ( xcoord.Dot(neutron) <= 0 )? "l":"r";
  return (p+n);//so right should be positive..
}

string GetDuoCode( XYZVector &beam, XYZVector &muon, XYZVector &proton, XYZVector &neutron )
{
  double pnAngle = TMath::ACos( proton.Unit().Dot( neutron.Unit() ) ) *180/TMath::Pi() ;
  string duoname="";
  if(pnAngle>40) duoname="g40";
  else duoname="l40";
  return duoname;
}

int getRegions(double dthetaP, double dthetaR)
{
  if( abs(dthetaP) > 55 || dthetaR> 10 || -55 > dthetaR )          return -1;
  if( abs(dthetaP) < 10 && 10 > dthetaR && dthetaR > -10 )           return 0; //hydrogen
  else if ( abs(dthetaP) < 20 && abs(dthetaR)<20 )               return 1; //QE-C 1
  else if ( abs(dthetaP) < 30 && abs(dthetaR)<30 )              return 2; //QE-C 3
  else if ( abs(dthetaP) < 20 && dthetaR<-30 && dthetaR >=-40  ) return 3; //2p2h/RES 4
  else if ( abs(dthetaP) > 20 && dthetaR<-30 && dthetaR >=-40  ) return 4; //2p2h/RES 1
  else if ( abs(dthetaP) < 55 && dthetaR>-55 && dthetaR < -40  ) return 5; //2p2h/RES 2
  else return -1;
}

int parsePartType( int pdg )
{
  pdg = abs(pdg);
  if ( pdg == 22 ) return 0;
  if ( pdg == 11 ) return 1;
  if ( pdg == 111 ) return 2;
  if ( pdg == 211 ) return 3;
  if ( pdg == 2212 ) return 4;
  if ( pdg == 2112 ) return 5;
  if ( pdg == 13 ) return 6;
  return -1;
}

string returnPartName( unsigned int type )
{
  if(type == 0 ) return "#gamma";
  if(type == 1 ) return "e";
  if(type == 2 ) return "#pi0";
  if(type == 3 ) return "#pi+/-";
  if(type == 4 ) return "p";
  if(type == 5 ) return "n";
  if(type == 6 ) return "#mu";
  return "";
}

TH2D* Normalize( TH2D* h, bool column = true )
{
  int NX = h->GetNbinsX();
  int NY = h->GetNbinsY();
  TH2D* ret = (TH2D*) h->Clone( Form( "%s-norm-column_%d", h->GetName(), column ) );
  if( column )
  {
    for( int i = 0; i< NX; i++ )
    {
      int ibin = i+1;
      double sum = h->Integral( ibin, ibin, 1, NY );
      if( sum == 0 ) continue;
      for( int j = 0; j< NY; j++ )
      {
        int jbin = j+1;
        ret->SetBinContent( ibin, jbin, h->GetBinContent( ibin, jbin )/sum );
      }
    }
  }
  else
  {
    for( int i = 0; i< NY; i++ )
    {
      int ibin = i+1;
      double sum = h->Integral( 1, NX, ibin, ibin );
      if( sum == 0 ) continue;
      for( int j = 0; j< NX; j++ )
      {
        int jbin = j+1;
        ret->SetBinContent( jbin, ibin, h->GetBinContent( jbin, ibin )/sum );
      }
    }
  }
  return ret;
}


bool passQ2( double q2 )
{
  return (q2>minQ2 && q2<maxQ2);
}


#endif
