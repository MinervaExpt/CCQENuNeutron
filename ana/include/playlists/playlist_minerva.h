#include "TChain.h"

namespace CCQENU_ANA{

  void get_mc_files( TChain* tree , string playlist){
    ifstream f;
    string mybase = getenv("MY_CCQENU");
    if(playlist == "minervame1A")f.open(Form("%s/include/playlists/CCQENu_minervame1A_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1B")f.open(Form("%s/include/playlists/CCQENu_minervame1B_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1C")f.open(Form("%s/include/playlists/CCQENu_minervame1C_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1D")f.open(Form("%s/include/playlists/CCQENu_minervame1D_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1E")f.open(Form("%s/include/playlists/CCQENu_minervame1E_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1F")f.open(Form("%s/include/playlists/CCQENu_minervame1F_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1G")f.open(Form("%s/include/playlists/CCQENu_minervame1G_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1L")f.open(Form("%s/include/playlists/CCQENu_minervame1L_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1M")f.open(Form("%s/include/playlists/CCQENu_minervame1M_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1N")f.open(Form("%s/include/playlists/CCQENu_minervame1N_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1O")f.open(Form("%s/include/playlists/CCQENu_minervame1O_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1P")f.open(Form("%s/include/playlists/CCQENu_minervame1P_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame5A")f.open(Form("%s/include/playlists/CCQENu_minervame5A_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6A")f.open(Form("%s/include/playlists/CCQENu_minervame6A_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6B")f.open(Form("%s/include/playlists/CCQENu_minervame6B_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6C")f.open(Form("%s/include/playlists/CCQENu_minervame6C_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6D")f.open(Form("%s/include/playlists/CCQENu_minervame6D_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6E")f.open(Form("%s/include/playlists/CCQENu_minervame6E_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6F")f.open(Form("%s/include/playlists/CCQENu_minervame6F_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6G")f.open(Form("%s/include/playlists/CCQENu_minervame6G_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6H")f.open(Form("%s/include/playlists/CCQENu_minervame6H_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6I")f.open(Form("%s/include/playlists/CCQENu_minervame6I_MC_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6J")f.open(Form("%s/include/playlists/CCQENu_minervame6J_MC_Inextinguishable_merged.txt",mybase.c_str()));

    //else if(playlist == "minervame6H")f.open(Form("%s/include/playlists/modified_6H_list.txt",mybase.c_str()));
    //else if(playlist == "minervame6I")f.open(Form("%s/include/playlists/modified_6I_list.txt",mybase.c_str()));
    //else if(playlist == "minervame6J")f.open(Form("%s/include/playlists/modified_6J_list.txt",mybase.c_str()));
    else cout << "The specified playlist is not in the list for playlist_minerva.h" << endl;

    while( f.good() ){
      
      // read line
      string line;
      getline( f, line );

      // skip empty and comment lines
      if( line.empty() ) continue;
      if( line.find( "#" )!=string::npos ) continue;

      // add file into tree
      tree->Add( line.c_str() );

    }
  }

  void get_data_files( TChain* tree, string playlist ){

    ifstream f;
    string mybase = getenv("MY_CCQENU");

    if(playlist == "minervame1A")f.open(Form("%s/include/playlists/CCQENu_minervame1A_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1B")f.open(Form("%s/include/playlists/CCQENu_minervame1B_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1C")f.open(Form("%s/include/playlists/CCQENu_minervame1C_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1D")f.open(Form("%s/include/playlists/CCQENu_minervame1D_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1E")f.open(Form("%s/include/playlists/CCQENu_minervame1E_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1F")f.open(Form("%s/include/playlists/CCQENu_minervame1F_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1G")f.open(Form("%s/include/playlists/CCQENu_minervame1G_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1L")f.open(Form("%s/include/playlists/CCQENu_minervame1L_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1M")f.open(Form("%s/include/playlists/CCQENu_minervame1M_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1N")f.open(Form("%s/include/playlists/CCQENu_minervame1N_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1O")f.open(Form("%s/include/playlists/CCQENu_minervame1O_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame1P")f.open(Form("%s/include/playlists/CCQENu_minervame1P_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame5A")f.open(Form("%s/include/playlists/CCQENu_minervame5A_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6A")f.open(Form("%s/include/playlists/CCQENu_minervame6A_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6B")f.open(Form("%s/include/playlists/CCQENu_minervame6B_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6C")f.open(Form("%s/include/playlists/CCQENu_minervame6C_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6D")f.open(Form("%s/include/playlists/CCQENu_minervame6D_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6E")f.open(Form("%s/include/playlists/CCQENu_minervame6E_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6F")f.open(Form("%s/include/playlists/CCQENu_minervame6F_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6G")f.open(Form("%s/include/playlists/CCQENu_minervame6G_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6H")f.open(Form("%s/include/playlists/CCQENu_minervame6H_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6I")f.open(Form("%s/include/playlists/CCQENu_minervame6I_DATA_Inextinguishable_merged.txt",mybase.c_str()));
    else if(playlist == "minervame6J")f.open(Form("%s/include/playlists/CCQENu_minervame6J_DATA_Inextinguishable_merged.txt",mybase.c_str()));

    else cout << "The specified playlist is not in the list for playlist_minerva.h" << endl; 
    
    while( f.good() ){
      
      // read line
      string line;
      getline( f, line );
      
      // skip empty and comment lines
      if( line.empty() ) continue;
      if( line.find( "#" )!=string::npos ) continue;

      // add file into tree
      tree->Add( line.c_str() );

    }

  }

  void get_truth_files( TChain* tree, string playlist ){
    get_mc_files( tree , playlist ); 
  }

  void get_mc_files_merged( TChain* tree ){
    tree->Add("/minerva/data/users/kenai/minmoddepccqe/merged/MergedMCTree.root");
  }

  void get_data_files_merged( TChain* tree ){
    tree->Add("/minerva/data/users/kenai/minmoddepccqe/merged/MergedDataTree.root");
  }

  void get_truth_files_merged( TChain* tree ){
    tree->Add("/minerva/data/users/kenai/minmoddepccqe/merged/MergedTruthTree.root");
  }

};

