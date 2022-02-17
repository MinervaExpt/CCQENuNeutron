#include "include/CCQENuPlotUtils.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HistogramUtils.h"
#include "TParameter.h" 

using namespace CCQENU_ANA;


int CombineSamples(std::vector<string> input_files, string output_file) 
{

  //Need to open files
  const int num_input_files = input_files.size();
  TFile *f_input[num_input_files];
  // read files to get histogram
  for(int i=0;i<num_input_files;i++){
    string filename = input_files[i];
    cout << "Opening " << filename << endl;
    f_input[i] = new TFile(filename.c_str());
    if (f_input[i]->IsZombie() || f_input[i]->GetListOfKeys()->IsEmpty()){
      Error("Playlistadder","Could not get histogram ROOT file or it was empty.");
      return 1;
    }
  }
  cout << "Creating output file "<< output_file << endl;
  string outname = output_file;
  TFile *f_output = new TFile(outname.c_str(),"RECREATE");

  //need to loop over all keys in input_files. Only save the special pot/events etc from the first file.
  for(int i=0;i<num_input_files;i++){
    cout << "Adding in objects from file " << i << endl;
    f_input[i]->cd();
    TList *file0List = f_input[i]->GetListOfKeys();
    TIter next(file0List);
    TObject *ob = 0;
    while((ob=next())){
      string name = ob->GetName();
      string classname = ((TKey*)ob)->GetClassName();
      //      cout << i << "\t" << name << "\t" << classname << endl;
      if(classname=="PlotUtils::MnvH1D" || classname=="PlotUtils::MnvH2D"){
	f_output->cd();
	((TKey*)ob)->ReadObj()->Write();
      }
      else if(i==0){//not mnvh*d, is it first file (get other objects like sample, pot, n_events etc.)
	f_output->cd();
	((TKey*)ob)->ReadObj()->Write(name.c_str());
      }
    }//while obj
  }//for playlists

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
      "CombineSamples outfile_root  file_genie_variation_xsection\n\n"<<
      "\t-Output file name.\n"
      "\t-List of input files.\n"
      "\t*********************************************\n"<<

    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;

    return 0;
  }

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("CombineSamples");
  par.push_back( "test_q2.root");
  par.push_back( "test.root");
  std::vector<std::string> in;
  //! Set user parameters
  for( int i=0; i<argc; ++i){    
    cout << i << endl;
    if(i>1) in.push_back(argv[i]);
    else par.at(i) = argv[i];
   
  }

  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

    
  return CombineSamples(in, par[1] );
}
