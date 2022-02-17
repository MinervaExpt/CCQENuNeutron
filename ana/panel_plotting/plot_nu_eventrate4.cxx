#include "plot_nu_eventrate_regions_neutron_proton_func.cxx"

int main(int argc, char* argv[])
{
  string f_orig = argv[1];
  string f_region = argv[2];
  string f_signal_weighted = argv[3];
  string sideband = argv[4];
  nu = atoi( argv[5] );
  cout<<f_orig<<endl;
  if (argc == 7 ) nPlaylists = atoi(argv[6]);
  //draw scales with errors
  //
  //Get Files
  TFile* file_orig = new TFile( f_orig.c_str(),"read");
  TFile* file_region = new TFile( f_region.c_str(),"read");
  TFile* file_signal_weighted = new TFile( f_signal_weighted.c_str(),"read");


  //Draw 1D
  //DrawFits( file_signal_weighted );
  ////
  //string model = "";
  //draw3DChunk( f_orig, f_signal_weighted, sideband, true, model );//do bckfitting
  //draw3DChunk( f_orig, f_signal_weighted, sideband, false , model);//not do bckfitting

  //Draw regions
  //drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband); //do bck fitting
  //drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband ); //no bck fitting

  //drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband, "vtx"); //do bck fitting
  //drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband, "vtx" ); //no bck fitting

  //drawRegionAnglePlots( file_region, file_signal_weighted, true , sideband, "nonvtx"); //do bck fitting
  //drawRegionAnglePlots( file_region, file_signal_weighted, false, sideband, "nonvtx" ); //no bck fitting
  //

  PlotRecoil( file_orig, file_signal_weighted );
  //draw1DDistros( file_orig, file_signal_weighted, "" );

  //Diagnostics
  //NeutronDiagnostics( file_orig, file_signal_weighted );

  return 0;
}
