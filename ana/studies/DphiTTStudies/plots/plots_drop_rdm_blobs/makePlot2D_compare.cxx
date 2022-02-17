{
#include <vector>
#include <pair>
#include <utility>
#include <map>

  using namespace std;
  bool doRatio = false;
  using namespace PlotUtils;

  //double dropFrac[]={0., 0.1,0.2,0.3}; int nFrac = sizeof(dropFrac)/sizeof(dropFrac[0]);
  //double dropEnergy[]={0.,5., 10., 15.,20.};int nEnergy = sizeof(dropEnergy)/sizeof(dropEnergy[0]);
  double dropFrac[]={0.,0.2,0.3}; int nFrac = sizeof(dropFrac)/sizeof(dropFrac[0]);
  double dropEnergy[]={0.,10.};int nEnergy = sizeof(dropEnergy)/sizeof(dropEnergy[0]);

  string tag = "none";

  vector<double> DropParF;
  vector<double> DropParE;
  for( int i = 0; i<nFrac; i++ )
  {
    for( int j = 0; j< nEnergy; j++ )
    {
      double F = dropFrac[i];
      double E = dropEnergy[j];
      DropParF.push_back(F);
      DropParE.push_back(E);
    }
  }

  

  TFile* f = new TFile(Form("../../root_files/dbs_CombinedPlaylists.root",tag.c_str()),"read");
  TFile* f2 = new TFile(Form("../../root_files/dbs_rw-1_CombinedPlaylists.root",tag.c_str()),"read");
  string type = "mc";


  map<pair<double,double>, TH2D*> h_blobE, h_blobDist, h_blobN;
  map<pair<double,double>, TH2D*> h_blobE_comp, h_blobDist_comp, h_blobN_comp;
  map<pair<double,double>, TH2D*> h_blobE_data, h_blobDist_data, h_blobN_data;
  map<pair<double,double>, TH1D*> h_blobE1D, h_blobDist1D, h_blobN1D;
  TVector2 pot = (TVector2) f->Get("pot");
  int rebinQ2 = 3;
  double scale = pot->X()/pot->Y();



  for( int i = 0 ; i< DropParF.size(); i++ )
  {
    double drop_fra = DropParF[i];
    double drop_ene = DropParE[i];
    pair<double, double> w(drop_fra, drop_ene);
    cout<<i<<", "<<w.first<<", "<<w.second<<endl;


    string comp = "mc";
    h_blobE_comp[w] = (TH2D*) f2->Get( Form("h_blobE_q2qe_%0.2f_%0.1f_%s", w.first, w.second, comp.c_str() ) );
    h_blobE_comp[w]->GetXaxis()->SetTitle("Energy (MeV)");
    h_blobE_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobE_comp[w]->Rebin2D(1, rebinQ2);
    h_blobE_comp[w]->SetLineColor(i+1);
    htmp = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.2f_%0.1f_%s", w.first, w.second, comp.c_str() ) ); htmp->Rebin2D(1,rebinQ2); h_blobE_comp[w]->Divide(htmp);

    h_blobDist_comp[w] = (TH2D*) f2->Get( Form("h_blobDist_q2qe_%0.2f_%0.1f_%s",  w.first, w.second, comp.c_str() ) );
    h_blobDist_comp[w]->GetXaxis()->SetTitle("Dist (mm)");
    h_blobDist_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobDist_comp[w]->Rebin2D(5,rebinQ2);
    h_blobDist_comp[w]->SetLineColor(i+1);
    htmp = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.2f_%0.1f_%s", w.first, w.second, comp.c_str() ) ); htmp->Rebin2D(5,rebinQ2); h_blobDist_comp[w]->Divide(htmp);

    h_blobN_comp[w] = (TH2D*) f2->Get( Form("h_blobN3D_q2qe_%0.2f_%0.1f_%s",  w.first, w.second, comp.c_str() ) );
    h_blobN_comp[w]->GetXaxis()->SetTitle("N Blobs");
    h_blobN_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobN_comp[w]->Rebin2D(1,rebinQ2);
    h_blobN_comp[w]->SetLineColor(i+1);
    htmp = (TH2D*) f->Get( Form("h_blobN3D_q2qe_%0.2f_%0.1f_%s", w.first, w.second, comp.c_str() ) ); htmp->Rebin2D(1,rebinQ2); h_blobN_comp[w]->Divide(htmp);
  }


  TCanvas* c = new TCanvas("c","c");
  c->cd();

  string opt = "histl";
  TLegend *l = new TLegend( 0, 0, 1,1 );

  TLatex* text = new TLatex();
  text->SetTextColor(kRed);
  text->SetTextFont(82);
  text->SetTextSize(0.1);
  text->SetTextAlign(22);

  TLine * line  = new TLine(0,1,2500,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1);
  line->SetLineStyle(2);



  pair<double, double> DropPars0(0,0);

  int imin = 0;
  int NPlots = h_blobE_comp[DropPars0]->GetNbinsY();
  int nX = int(TMath::Sqrt(NPlots))+1;
  int n = DropParF.size();
  cout<<"Dividing canvas: "<<nX<<", "<<nX-1<<endl;
  c->Divide(nX,nX-1,0,0);

  TH2D* h_blobE_base = (TH2D*) h_blobE_comp[DropPars0]->Clone("h_blobE_base");
  for( int i = 0; i< n; i++ )
  {
    pair<double, double> w(DropParF[i], DropParE[i]);
    if(i!=0 && (w.first == 0 || w.second == 0 )) continue;
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_blobE_comp[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first,w.second,j+1), j+1,j+1);
      h->SetLineColor(i+1);
      h->GetYaxis()->SetRangeUser(.8,1.2);
      c->cd(j+1);
      gPad->SetGrid();
      cout<<j+1<<endl;
      h->Draw(opt.c_str() );
    }

    l->AddEntry(h_blobE_comp[w], Form("%.2f_%.1f",w.first, w.second ) );
    if ( i==0 ) opt+="same";
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  c->Print(Form( "Compare_blobE2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );

  TH2D* h_blobDist_base = (TH2D*) h_blobDist_comp[DropPars0]->Clone("h_blobDist_base");
  opt="histl";
  for( int i = 0; i< n; i++ )
  {

    pair<double, double> w(DropParF[i], DropParE[i]);
    if(i!=0 && (w.first == 0 || w.second == 0) ) continue;
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_blobDist_comp[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first,w.second,j+1), j+1,j+1);
      h->SetLineColor(i+1);
      h->GetYaxis()->SetRangeUser(.95,1.05);
      h->GetXaxis()->SetRangeUser(0,2500);
      c->cd(j+1);
      gPad->SetGrid();
      cout<<j+1<<endl;
      h->Draw(opt.c_str() );
    }

    if ( i==0) opt+="same";
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  c->Print(Form( "Compare_blobDist2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );

  TH2D* h_blobN_base = (TH2D*) h_blobN_comp[DropPars0]->Clone("h_blobN_base");
  opt="hist";
  for( int i = 0; i< n; i++ )
  {
    pair<double, double> w(DropParF[i], DropParE[i]);
    if(i!=0 &&(w.first == 0 || w.second == 0 ) ) continue;
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_blobN_comp[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first,w.second,j+1), j+1,j+1);
      h->SetLineColor(i+1);
      h->GetYaxis()->SetRangeUser(0,2);
      h->GetXaxis()->SetRangeUser(0,7);
      h->SetMarkerSize(0);
      c->cd(j+1);
      gPad->SetGrid();
      cout<<j+1<<endl;
      h->Draw(opt.c_str() );
    }

    if ( i==0 ) opt+="same";
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  c->Print(Form( "Compare_blobN2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );





}



