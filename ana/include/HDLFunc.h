#ifndef __HDL_FUNC_H
#define __HDL_FUNC_H
#include "include/CCQENuUtils.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH1D.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"


using namespace CCQENU_ANA;


class HDLHistos
{
  public:
  HDLHistos( HyperDimLinearizer *hdl, axis_binning& xbin, axis_binning& ybin, axis_binning &zbin, MnvH2D** results, bool includeSys );
  ~HDLHistos();

  MnvH1D* ProjectionX(vector<MnvH2D*> &slices, string name = "_py", int ymin = 0, int ymax = -1, int zmin = 0, int zmax = -1 );
  MnvH1D* ProjectionY(vector<MnvH2D*> &slices, string name = "_px", int xmin = 0, int xmax = -1, int zmin = 0, int zmax = -1 );
  MnvH1D* ProjectionZ(vector<MnvH2D*> &slices, string name = "_pz", int xmin = 0, int xmax = -1, int ymin = 0, int ymax = -1 );

  vector<MnvH1D*> ProjectionsX( vector<vector<MnvH2D*>> &slices, string name, int ymin=0, int ymax=-1, int zmin=0, int zmax=-1 );
  vector<MnvH1D*> ProjectionsY( vector<vector<MnvH2D*>> &slices, string name, int xmin=0, int xmax=-1, int zmin=0, int zmax=-1 );
  vector<MnvH1D*> ProjectionsZ( vector<vector<MnvH2D*>> &slices, string name, int xmin=0, int xmax=-1, int ymin=0, int ymax=-1 );



  MnvH2D* ProjectionXZ (vector<MnvH2D*> &slices, string name="_pxz", int ymin=0, int ymax=-1 );
  MnvH2D* ProjectionXY (vector<MnvH2D*> &slices, string name="_pxy", int zmin=0, int zmax=-1 );
  MnvH2D* ProjectionYZ (vector<MnvH2D*> &slices, string name="_pyz", int xmin=0, int xmax=-1 );

  vector<MnvH2D*> ProjectionsXZ (vector<vector<MnvH2D*>> &slices, string name="_pxz", int ymin=0, int ymax=-1 );
  vector<MnvH2D*> ProjectionsXY (vector<vector<MnvH2D*>> &slices, string name="_pxy", int zmin=0, int zmax=-1 );
  vector<MnvH2D*> ProjectionsYZ (vector<vector<MnvH2D*>> &slices, string name="_pyz", int xmin=0, int xmax=-1 );

  void Add( vector<MnvH1D*> &h1, vector<MnvH1D*> &h2 ); //add h2 to h1
  void Add( vector<MnvH2D*> &h1, vector<MnvH2D*> &h2 ); //add h2 to h1

  vector<vector<MnvH2D*>> GetHistos(){ return this->histos; }

  TAxis* GetXaxis() { return this->axes->GetXaxis(); }
  TAxis* GetYaxis() { return this->axes->GetYaxis(); }
  TAxis* GetZaxis() { return this->axes->GetZaxis(); }

  TH3D* GetAxisHisto() { return this->axes; }
  private:
  vector<vector<MnvH2D*>> histos;
  axis_binning xbin, ybin, zbin;
  TH3D* axes;
};

HDLHistos::HDLHistos(HyperDimLinearizer *hdl, axis_binning& xbin, axis_binning& ybin, axis_binning &zbin, MnvH2D** results, bool includeSys )
{
  //Fill Histos
  histos.clear();
  for( unsigned int i = 0; i< nHistos; i++ )
  {
    histos.push_back( hdl->Get2DMnvHistos( results[i], includeSys ));
  }

  //Set axis
  this->axes = new TH3D( Form("%s_axes",results[0]->GetName() ),Form("%s_axes",results[0]->GetName() ), xbin.nbins, &(xbin.bin_edges[0]),  ybin.nbins, &(ybin.bin_edges[0]),  zbin.nbins, &(zbin.bin_edges[0])); 
  this->xbin = xbin;
  this->ybin = ybin;
  this->zbin = zbin;

}

HDLHistos::~HDLHistos()
{
  delete this->axes;
  for(uint i = 0; i<histos.size();i++)
  {
    for(uint j = 0; j< histos[i].size();j++)
      delete histos[i][j];

  }
}



MnvH1D* HDLHistos::ProjectionX( vector<MnvH2D*> &slices, string name, int ymin, int ymax, int zmin, int zmax )
{
  //vector<MnvH2D*> slices = this->histos[ihisto]; // i.e.  q2qe vs dthetaP in vectors of dthetaR
  int Zmin = (zmax >= zmin )? zmin : 1;
  int Zmax = (zmax >= zmin )? zmax : axes->GetNbinsZ();
  if ( Zmax > axes->GetNbinsZ() ) Zmax = axes->GetNbinsZ();
  MnvH1D* ret = (MnvH1D*) slices[Zmin]->ProjectionX(name.c_str(), ymin, ymax);
  for( int z = Zmin+1; z<=Zmax;z++ ) 
  {
    MnvH1D* proj = (MnvH1D*)slices[z]->ProjectionX(name.c_str(), ymin, ymax) ;
    ret->Add( proj );
    delete proj;
  }
  return ret;
}

MnvH1D* HDLHistos::ProjectionY( vector<MnvH2D*> &slices, string name, int xmin, int xmax, int zmin, int zmax )
{
  //if( ihisto < 0 || ihisto > nHistos ){ return vector<MnvH1D*>(); }
  //vector<MnvH2D*> slices = this->histos[ihisto]; // i.e.  q2qe vs dthetaP in vectors of dthetaR
  int Zmin = (zmax >= zmin )? zmin : 1;
  int Zmax = (zmax >= zmin )? zmax : axes->GetNbinsZ();
  if ( Zmax > axes->GetNbinsZ() ) Zmax = axes->GetNbinsZ();

  MnvH1D* ret = new MnvH1D(name.c_str(), name.c_str(), this->ybin.nbins, &(this->ybin.bin_edges[0]) );
  ret->AddMissingErrorBandsAndFillWithCV( *(slices[0]) );

  for( int z = Zmin; z<=Zmax;z++ ) 
  {
    MnvH1D* proj = (MnvH1D*)slices[z]->ProjectionY(name.c_str(), xmin, xmax) ;
    ret->Add( proj );
    delete proj;

  }
  return ret;
}

MnvH1D* HDLHistos::ProjectionZ( vector<MnvH2D*> &slices, string name, int xmin, int xmax, int ymin, int ymax )
{
  //if( ihisto < 0 || ihisto > nHistos ){ return vector<MnvH1D*>(); }
  //vector<MnvH2D*> slices = this->histos[ihisto]; // i.e.  q2qe vs dthetaP in vectors of dthetaR
  MnvH1D* ret = new MnvH1D(name.c_str(), name.c_str(), this->zbin.nbins, &(this->zbin.bin_edges[0]) );
  ret->AddMissingErrorBandsAndFillWithCV( *(slices[0]) );

  double integral, error;
  for( int z = 1; z <= axes->GetNbinsZ(); z++ )
  {
    //cv
    integral = slices[z]->IntegralAndError(xmin,xmax, ymin, ymax, error );
    ret->SetBinContent( z, integral);
    ret->SetBinError( z, error);
  }

  vector<string> vertNames = ret->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    //cout<<vertNames[k]<<",";
    unsigned int nUniverses = ret->GetVertErrorBand( vertNames[k] )->GetNHists();
    ret->AddVertErrorBand(vertNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH1D* r = ret->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        TH2D* s = slices[z]->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
        integral = s->IntegralAndError(xmin,xmax, ymin,ymax, error);
        r->SetBinContent( z, integral );
        r->SetBinError( z, error);
      }
    }
  }
  vector<string> latNames = ret->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    //cout<<latNames[k]<<",";
    unsigned int nUniverses = ret->GetLatErrorBand( latNames[k] )->GetNHists();
    ret->AddLatErrorBand(latNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH1D* r = ret->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        TH2D* s = slices[z]->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
        integral = s->IntegralAndError(xmin,xmax, ymin,ymax, error);
        r->SetBinContent( z, integral );
        r->SetBinError( z, error);
      }
    }
  }
  return ret;
}


MnvH2D* HDLHistos::ProjectionXZ (vector<MnvH2D*> &slices, string name, int ymin, int ymax )
{
  MnvH2D* ret = new MnvH2D(name.c_str(), name.c_str(), this->xbin.nbins, &(this->xbin.bin_edges[0]), this->zbin.nbins, &(this->zbin.bin_edges[0]) );
  ret->AddMissingErrorBandsAndFillWithCV( *(slices[0]) );

  double integral, error;
  for( int z = 1; z <= axes->GetNbinsZ(); z++ )
  {
    for( int x = 1; x <= axes->GetNbinsX(); x++ )
    {
      //cv
      integral = slices[z]->IntegralAndError(x,x, -1, 0, error );
      ret->SetBinContent( x,z, integral);
      ret->SetBinError( x,z, error);
    }
  }

  vector<string> vertNames = ret->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    //cout<<vertNames[k]<<",";
    unsigned int nUniverses = ret->GetVertErrorBand( vertNames[k] )->GetNHists();
    ret->AddVertErrorBand(vertNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* r = ret->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        TH2D* s = slices[z]->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
        for( int x = 1; x<axes->GetNbinsX(); x++ )
        {
          integral = s->IntegralAndError(x,x, 0,-1, error);
          r->SetBinContent( x,z, integral );
          r->SetBinError( x,z, error);
        }
      }
    }
  }
  vector<string> latNames = ret->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    //cout<<latNames[k]<<",";
    unsigned int nUniverses = ret->GetLatErrorBand( latNames[k] )->GetNHists();
    ret->AddLatErrorBand(latNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* r = ret->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        for( int x = 1; x<=axes->GetNbinsX();x++)
        {
          TH2D* s = slices[z]->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
          integral = s->IntegralAndError(x,x,0,-1, error);
          r->SetBinContent( x,z, integral );
          r->SetBinError( x,z, error);
        }
      }
    }
  }
  return ret;
}

MnvH2D* HDLHistos::ProjectionYZ (vector<MnvH2D*> &slices, string name, int xmin, int xmax )
{
  MnvH2D* ret = new MnvH2D(name.c_str(), name.c_str(), this->ybin.nbins, &(this->ybin.bin_edges[0]), this->zbin.nbins, &(this->zbin.bin_edges[0]) );
  ret->AddMissingErrorBandsAndFillWithCV( *(slices[0]) );

  double integral, error;
  for( int z = 1; z <= axes->GetNbinsZ(); z++ )
  {
    for( int y = 1; y <= axes->GetNbinsX(); y++ )
    {
      //cv
      integral = slices[z]->IntegralAndError(-1, 0,y,y,error );
      ret->SetBinContent( y,z, integral);
      ret->SetBinError( y,z, error);
    }
  }

  vector<string> vertNames = ret->GetVertErrorBandNames();
  for( unsigned int k = 0; k < vertNames.size(); ++k )
  {
    //cout<<vertNames[k]<<",";
    unsigned int nUniverses = ret->GetVertErrorBand( vertNames[k] )->GetNHists();
    ret->AddVertErrorBand(vertNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* r = ret->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        TH2D* s = slices[z]->GetVertErrorBand( vertNames[k] )->GetHist(iUniv);
        for( int y = 1; y<axes->GetNbinsX(); y++ )
        {
          integral = s->IntegralAndError(0,-1, y,y, error);
          r->SetBinContent( y,z, integral );
          r->SetBinError( y,z, error);
        }
      }
    }
  }
  vector<string> latNames = ret->GetLatErrorBandNames();
  for( unsigned int k = 0; k < latNames.size(); ++k )
  {
    //cout<<latNames[k]<<",";
    unsigned int nUniverses = ret->GetLatErrorBand( latNames[k] )->GetNHists();
    ret->AddLatErrorBand(latNames[k],nUniverses);
    for (unsigned int iUniv = 0 ;iUniv< nUniverses; iUniv++) {
      TH2D* r = ret->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
      for( int z = 1; z<=axes->GetNbinsZ(); z++ )
      {
        for( int y = 1; y<=axes->GetNbinsX();y++)
        {
          TH2D* s = slices[z]->GetLatErrorBand( latNames[k] )->GetHist(iUniv);
          integral = s->IntegralAndError(0,-1,y,y,error);
          r->SetBinContent( y,z, integral );
          r->SetBinError( y,z, error);
        }
      }
    }
  }
  return ret;
}
MnvH2D* HDLHistos::ProjectionXY (vector<MnvH2D*> &slices, string name, int zmin, int zmax )
{
  int Zmin = (zmax >= zmin )? zmin : 1;
  int Zmax = (zmax >= zmin )? zmax : axes->GetNbinsZ();
  Zmax = min(axes->GetNbinsZ(), Zmax);
  Zmin = max(Zmin,0);
  MnvH2D* ret = (MnvH2D*) slices[Zmin]->Clone(name.c_str());
  for( int z = Zmin+1; z <= Zmax; z++ )
  {
    ret->Add( slices[z] );
  }
  return ret;
}




vector<MnvH1D*> HDLHistos::ProjectionsX( vector<vector<MnvH2D*>> &slices, string name, int ymin, int ymax, int zmin, int zmax )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionX( slices[i],Form("%s_%s",name.c_str(), names[i].c_str()), ymin, ymax, zmin, zmax ) ) ;
  return ret;

}
vector<MnvH1D*> HDLHistos::ProjectionsY( vector<vector<MnvH2D*>> &slices, string name, int xmin, int xmax, int zmin, int zmax )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionY( slices[i], Form("%s_%s",name.c_str(), names[i].c_str()), xmin, xmax, zmin, zmax ) ) ;
  return ret;

}
vector<MnvH1D*> HDLHistos::ProjectionsZ( vector<vector<MnvH2D*>> &slices, string name, int xmin, int xmax, int ymin, int ymax )
{
  vector<MnvH1D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionZ( slices[i],Form("%s_%s",name.c_str(), names[i].c_str()), xmin, xmax, ymin, ymax ) ) ;
  return ret;

}


vector<MnvH2D*> HDLHistos::ProjectionsXY( vector<vector<MnvH2D*>> &slices, string name, int zmin, int zmax )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionXY( slices[i],Form("%s_%s",name.c_str(), names[i].c_str()), zmin, zmax ) ) ;
  return ret;
}


vector<MnvH2D*> HDLHistos::ProjectionsXZ( vector<vector<MnvH2D*>> &slices, string name, int ymin, int ymax )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionXZ( slices[i],Form("%s_%s",name.c_str(), names[i].c_str()), ymin, ymax ) ) ;
  return ret;
}

vector<MnvH2D*> HDLHistos::ProjectionsYZ( vector<vector<MnvH2D*>> &slices, string name, int xmin, int xmax )
{
  vector<MnvH2D*> ret;
  for( unsigned int i = 0; i< slices.size(); i++ ) ret.push_back( ProjectionYZ( slices[i],Form("%s_%s",name.c_str(), names[i].c_str()), xmin, xmax ) ) ;
  return ret;
}

void HDLHistos::Add( vector<MnvH1D*> &h1, vector<MnvH1D*> &h2 )
{
  for(unsigned int i = 0; i< h1.size(); i++ ) h1[i]->Add( h2[i] );
}


void HDLHistos::Add( vector<MnvH2D*> &h1, vector<MnvH2D*> &h2 )
{
  for(unsigned int i = 0; i< h1.size(); i++ ) h1[i]->Add( h2[i] );
}
#endif
