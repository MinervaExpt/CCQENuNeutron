#include "include/CCQENuPlotUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HistogramUtils.h"
#include "TParameter.h" 

using namespace CCQENU_ANA;

double GetPOTScale(TFile *file){
  TVector2 *m1 = (TVector2*)file->Get("pot");
  double pot_data = m1->X();
  double pot_mc = m1->Y();
  double scale = pot_data/pot_mc;
  delete m1;
  if(scale > 1) cout << "I expected more POT in MC than data. Make sure you ran your jobs correctly or know this is true" << endl;
  //  cout <<"Scale of the standard playlists " << scale << endl;
  return scale;
}

double GetPOTScale(TFile *final_sample, TFile *mconly_sample){
  //Check if pot in final_sample is already normalized
  TVector2 *m1 = (TVector2*)final_sample->Get("pot");
  double pot_data_final_sample = m1->X();
  double pot_mc_final_sample = m1->Y();
  if(pot_data_final_sample!=pot_mc_final_sample){
    cout << "You are attempting to add a mconly sample (such as a 2p2h sample) to a final sample which is NOT pot normalized" << endl;
    cout << "Final sample POT (data,mc): ("<<pot_data_final_sample << "," << pot_mc_final_sample << ")" << endl;
    exit (EXIT_FAILURE);
  }
  //
 TVector2 *m2 = (TVector2*)mconly_sample->Get("pot");
  double pot_mc_mconly_sample = m2->Y();
  double scale = pot_mc_final_sample/pot_mc_mconly_sample;
  cout <<"Scale of the special playlists " << scale<<"\t" <<  pot_mc_final_sample << "\t" <<pot_mc_mconly_sample  << endl;
  delete m1;
  delete m2;

  if(scale > 1) cout << "I expected more POT in MC than data. Make sure you ran your jobs correctly or know this is true" << endl;
  return scale;

}


int PlaylistAdder(string basefilename, string stage, std::vector<string> playlists, bool mconly) 
{

  //Make sure input stage is valid
  bool validStage = false;
  if(stage=="MuonEventSelection"||stage=="EffPurity"||stage=="BackgroundTemplates"||stage=="Migration"||stage=="ProtonEfficiencyStudy") validStage=true;
  if(!validStage){
    cout << "You have provided an invalid stage to add. Valid stages are:MuonEventSelection, EffPurity, BackgroundTemples, Migration" << endl;
    return 1;
  }

  //Need to open files
  const int num_playlists = playlists.size()-4;//Full par input vector less 3 header values
  if(mconly and num_playlists !=2){
    cout << "You are passing a mconly sample. There should just be two playlists the final_sample and the mconly_sample. Passed in that order" << endl;
    exit (EXIT_FAILURE);
  }
  TFile *f_input[num_playlists];
  // read files to get histogram
  for(int i=0;i<num_playlists;i++){
    string filename = basefilename+playlists[i+4]+".root";
    cout << "Opening " << filename << endl;
    f_input[i] = new TFile(filename.c_str());
  
    if (f_input[i]->IsZombie() || f_input[i]->GetListOfKeys()->IsEmpty()){
      Error("Playlistadder","Could not get histogram ROOT file or it was empty.");
      return 1;
    }
  }
  string outname = "";
  if(!mconly) outname = basefilename+"CombinedPlaylists.root";
  else outname = basefilename+"SpecialSampleIncluded.root";
  TFile *f_output = new TFile(outname.c_str(),"RECREATE");
  
  //Do I have to carry extra variables, pot normalize or other specific items

  //Do I POT normalize? Might be that I always do, but keep option open to not renormalize.
  bool POTNormalize = false;
  if(stage=="MuonEventSelection") POTNormalize = true;
  else if(stage=="EffPurity") POTNormalize = true;//Also need to recalc eff from raw distributions....no data POT
  else if(stage=="BackgroundTemplates") POTNormalize=true;
  else if(stage=="Migration") POTNormalize = true;
  else if (stage=="ProtonEfficiencyStudy") POTNormalize = true;
  //Do I need to carry n_events
  bool CombineNEvents = false;
  if(stage=="MuonEventSelection") CombineNEvents = true;
  else if(stage=="EffPurity") CombineNEvents = true;
  else if(stage=="BackgroundTemplates") CombineNEvents=true;
  else if(stage=="Migration") CombineNEvents = true;
  else if (stage=="ProtonEfficiencyStudy") CombineNEvents = false;
  
  //Do I need to carry sample
  bool SampleObject = false;
  if(stage=="MuonEventSelection") SampleObject = true;
  else if(stage=="EffPurity") SampleObject = false;
  else if(stage=="BackgroundTemplates") SampleObject=true;
  else if(stage=="Migration") SampleObject = true;
  else if(stage=="ProtonEfficiencyStudy") SampleObject = false;

  //Do I need to carry flux constraint
  bool fluxConstraint = false;
  if(stage=="MuonEventSelection") fluxConstraint = true;
  else if(stage=="EffPurity") fluxConstraint = true;
  else if(stage=="BackgroundTemplates") fluxConstraint=true;
  else if(stage=="Migration") fluxConstraint = true;
  else if(stage=="ProtonEfficiencyStudy") fluxConstraint = false;

  //Addition stuff
  //loop over objects in first input file. Grab the name. figure out if it is a special case or the typical mnvh1d/2d. 
  //Get objects of same name from other input files and add to the original object. MAYBE doing operations like POT normalization

//vector stuff
  float  pot_x = 0;
  float  pot_y = 0;
  float  events_x = 0;
  float  events_y = 0;
  //Add hist
  TList *file0List = f_input[0]->GetListOfKeys();
  TIter next(file0List);
  TObject *ob = 0;
  while((ob=next())){
    string name = ob->GetName();
    string classname = ((TKey*)ob)->GetClassName();
    if(name == "pot"){
      for(int i=0;i<num_playlists;i++){
	TVector2 *m1 = (TVector2*)f_input[i]->Get("pot");
	if(!mconly||i==0){//These are for samples which are additional bkg/signal which need to scaled to the final sample pot, but do not contribute to pot
	  pot_x+=m1->X();//data
	  pot_y+=m1->Y();//mc
	}
	delete m1;
      }  
    }
    else if(name == "n_events" && CombineNEvents){
      for(int i=0;i<num_playlists;i++){
	TVector2 *m1 = (TVector2*)f_input[i]->Get("n_events");
	if(!mconly||i==0)events_x+=m1->X();//data
	events_y+=m1->Y();//mc
	delete m1;
      }
    }
    else{
      
      if(classname=="PlotUtils::MnvH2D"){
	//	if(!mconly)cout << "2D Histogram\t" << name << endl;
	PlotUtils::MnvH2D *hist2D = (PlotUtils::MnvH2D*)f_input[0]->Get(name.c_str());
	if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)hist2D->Scale(GetPOTScale(f_input[0]));//not special sample and don't find data
	for(int i=1;i<num_playlists;i++){
	  PlotUtils::MnvH2D *tmp2D = (PlotUtils::MnvH2D*)f_input[i]->Get(name.c_str());
	  //	  cout << f_input[i]->GetName() << "\t" << tmp2D << endl;
	  if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp2D->Scale(GetPOTScale(f_input[i]));//not special sample and don't find data
	  if(mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp2D->Scale(GetPOTScale(f_input[0],f_input[i])); // scale the mconly_sample by correct pot
	  if(mconly&&POTNormalize&&name.find("data") != std::string::npos) continue;//skip data histograms for mconly sample
	  hist2D->Add(tmp2D);
	  delete tmp2D;
	}
	f_output->cd();
	hist2D->Write(name.c_str());
	delete hist2D;
      }
      else if(classname=="PlotUtils::MnvH1D"){
	//	if(!mconly)cout << "1D Histogram\t" << name <<"\t" <<name.find("data") << endl;
	PlotUtils::MnvH1D *hist1D = (PlotUtils::MnvH1D*)f_input[0]->Get(name.c_str());
	if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)hist1D->Scale(GetPOTScale(f_input[0]));

	for(int i=1;i<num_playlists;i++){
	  PlotUtils::MnvH1D *tmp1D = (PlotUtils::MnvH1D*)f_input[i]->Get(name.c_str());
	  if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp1D->Scale(GetPOTScale(f_input[i]));
	  if(mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp1D->Scale(GetPOTScale(f_input[0],f_input[i])); // scale the mconly_sample by correct pot
	  if(mconly&&POTNormalize&&name.find("data") != std::string::npos) continue;//skip data histograms for mconly sample
	  hist1D->Add(tmp1D);
	  delete tmp1D;
	}
	f_output->cd();
	hist1D->Write();
	delete hist1D;
      }
      else if(classname=="PlotUtils::MnvH3D"){
	//	if(!mconly)cout << "3D Histogram\t" << name << endl;
	PlotUtils::MnvH3D *hist3D = (PlotUtils::MnvH3D*)f_input[0]->Get(name.c_str());
	if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)hist3D->Scale(GetPOTScale(f_input[0]));//not special sample and don't find data
	for(int i=1;i<num_playlists;i++){
	  PlotUtils::MnvH3D *tmp3D = (PlotUtils::MnvH3D*)f_input[i]->Get(name.c_str());
	  //	  cout << f_input[i]->GetName() << "\t" << tmp3D << endl;
	  if(!mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp3D->Scale(GetPOTScale(f_input[i]));//not special sample and don't find data
	  if(mconly&&POTNormalize&&name.find("data") == std::string::npos)tmp3D->Scale(GetPOTScale(f_input[0],f_input[i])); // scale the mconly_sample by correct pot
	  if(mconly&&POTNormalize&&name.find("data") != std::string::npos) continue;//skip data histograms for mconly sample
	  hist3D->Add(tmp3D);
	  delete tmp3D;
	}
	f_output->cd();
	hist3D->Write(name.c_str());
	delete hist3D;
      }

      else{
	cout << "I don't know what this is: "<< classname << "\tskipping" << endl;
      }
    }//else
  }//obj next
  ob = NULL;
  file0List = NULL;
  //  delete ob;
  //  delete file0List;
  cout << pot_x <<"\t" << pot_y << endl;
  cout << "Single one off writes" << endl;
  //Do single one off operations
  //copy sample if needed
  //copy fluxconstraint object as necessary
  if(SampleObject){
    cout << "Writing sample object" << endl;
    TMap *map = (TMap*)f_input[0]->Get("sample");
    map->Print();
    f_output->cd();
    f_output->WriteTObject( map, "sample" );
    delete map;
  }
  /*
  if(fluxConstraint){
    cout << "Writing fluxconstrain object" << endl;
    string tmpName = "appliedFluxConstraint";
    TParameter<bool> *fluxConstrain = (TParameter<bool>*)f_input[0]->Get(tmpName.c_str());
    if(fluxConstrain){
      f_output->cd();
      f_output->WriteTObject( fluxConstrain );
    }
    else{
      tmpName = "fluxConstraintHisto";
      fluxConstrain = (TParameter<bool>*)f_input[0]->Get(tmpName.c_str());
      if(fluxConstrain){
	f_output->cd();
	f_output->WriteTObject( fluxConstrain );
      }
      else{
	cout << "Your flux constraint object has an unexpected name. Please check" << endl;
      }
    }

    delete fluxConstrain;
  }
  */
  if(POTNormalize){
    cout << "Writing pot object" << endl;
    cout << "Sample is POT normalized. The final data POT is: " << pot_x << endl;
    TVector2 pot(pot_x,pot_x);
    f_output->cd();
    pot.Write("pot");
  }
  else{ 
    cout << "Sample is NOT POT normalized. The final data POT is: " << pot_x << " and the final MC POT is " << pot_y<< endl;
    TVector2 pot(pot_x,pot_y);
    f_output->cd();
    pot.Write("pot");
  }
  if(CombineNEvents){
    cout << "Writing n_event object" << endl;
    TVector2 n_events(events_x,events_y);
    cout << "Final event count for data is: " << events_x << " and for MC: " << events_y << endl;
    f_output->cd();
    n_events.Write("n_events");
  }
  
  f_output->Close();
  return 0;
};

  int main( int argc, char *argv[])
{
  ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);
  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "PlaylistAdder outfile_root  file_genie_variation_xsection\n\n"<<
      "\t-Base string of files. Assumes files are of the form XXXXXX_minervaY.root where Y is any playlist(s).\n"<<
      "\t-Type of file we are adding. Valid choices: MuonEventSelection, EffPurity, BackgroundTemples, Migration\n"<<
      "\t-MCOnly sample. Valid choices: 0 = no special mconly sample (2p2h for instance), 1 = yes\n"<<
      "\t-List of playlists. Example would be minerva1 minerva7 minerva 9. This will add 1,7,9 together.\n"
      "\t*********************************************\n"<<

    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;

    return 0;
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("PlaylistAdder");
  par.push_back( "MuonSelection_");
  par.push_back( "MuonSelection");
  //! Set user parameters
  for( int i=0; i<argc; ++i){    
    cout << i << endl;
    if(i>2) par.push_back(argv[i]);
    else par.at(i) = argv[i];
   
  }
  bool mconly = false;
  if(atoi(par[3].c_str())!=0) mconly=true;
  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

    
  return PlaylistAdder(par[1], par[2], par, mconly);
}
