#include "include/GlobalParameters.h"
#include "include/CCQENuUtils.h"
#include "TParameter.h" 
#include "PlotUtils/HyperDimLinearizer.h"

#include "include/GeneralFunc.h" //in CCQENuNeutron/include/
#include "include/EventCuts.h"
#include "include/CommonBins.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
using namespace CCQENU_ANA;


//forward declaration
//HyperDimLinearizer* DeclareHDLHistos2D( CCQENuUtils* utils, MnvH2D **h, std::string name, std::string title, axis_binning &xbins, axis_binning &ybins, axis_binning &zbins, axis_binning &global_x );
//void FillHDL2D( std::string name, CCQENuUtils *utils, double x, double y, double z, bool isData, CCQENuEvent* evt, double wgt );
//================================================

//int npMode = 1;//0:neutron, 1:proton, 10:neutron+proton




int main( int argc, char *argv[])
{
  XYZVector vtx(0,0,0);
  XYZVector coneAxis(1,0,0);
  vector<XYZVector> Positions({ XYZVector(10,0,0), XYZVector(10,5,0) } );

  double angle = 10;
  Cone cone(vtx, coneAxis, angle, 99999999 );
  cout<<"=========== cone angle: "<<angle<<endl;
  for( auto pos: Positions )
  {
    cout<<ACos( pos.Unit().Dot(cone.GetAxis().Unit() ) )*180/TMath::Pi()<<", IsInsideCone "<<cone.InsideCone( pos, -1 )<<endl;
    cout<<"Simulate"<<endl;
    cout<<"Default angle: "<<cone.GetDefaultAngleRadian()<<endl;

    angle = cone.GetDefaultAngleRadian();
    //else angle*=TMath::Pi()/180;
    XYZVector posWRTvtx = pos - cone.GetVtx();
    cout<<"Vtx: "<<cone.GetVtx().X()<<", "<<cone.GetVtx().Y()<<", "<<cone.GetVtx().Z()<<endl;
    double z = posWRTvtx.Dot( cone.GetAxis() );
    if ( z < 0 || z > 9999999999999 ) cout<<"Too far"<<endl;
    double cosTheta = posWRTvtx.Unit().Dot( cone.GetAxis().Unit() );
    double theta = TMath::ACos( cosTheta );
    cout<<"Theta: "<<theta<<", "<<theta*180/TMath::Pi()<<" angle: "<<angle<<endl;
    cout<< (abs(theta) < angle)<<endl;
  }
  return 0;
}
