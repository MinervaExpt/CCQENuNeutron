{
#include <vector>
#include <pair>
#include <utility>
#include <map>

  using namespace std;
  bool doRatio = false;
  using namespace PlotUtils;

  double dropFrac[]={0., 0.1,0.2,0.3}; int nFrac = sizeof(dropFrac)/sizeof(dropFrac[0]);
  double dropEnergy[]={0.,5., 10., 15.,20.};int nEnergy = sizeof(dropEnergy)/sizeof(dropEnergy[0]);

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
  string type = "mc";


  map<pair<double,double>, TH2D*> h_blobE, h_blobDist, h_nBlobs;
  map<pair<double,double>, TH2D*> h_blobE_comp, h_blobDist_comp, h_nBlobs_comp;
  map<pair<double,double>, TH2D*> h_blobE_data, h_blobDist_data, h_nBlobs_data;
  map<pair<double,double>, TH1D*> h_blobE1D, h_blobDist1D, h_nBlobs1D;
  TVector2 pot = (TVector2) f->Get("pot");
  int rebinQ2 = 3;
  double scale = pot->X()/pot->Y();


  for( int i = 0 ; i< DropParF.size(); i++ )
  {
    double drop_fra = DropParF[i];
    double drop_ene = DropParE[i];
    pair<double, double> w(drop_fra, drop_ene);
    cout<<i<<", "<<w.first<<", "<<w.second<<endl;

    h_blobE_data[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.2f_%.1f_data", w.first, w.second ) );
    h_blobE_data[w]->GetXaxis()->SetTitle("Energy (MeV)");
    h_blobE_data[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobE_data[w]->Rebin2D(1, rebinQ2);

    h_blobDist_data[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.2f_%0.1f_data",  w.first, w.second ) );
    h_blobDist_data[w]->GetXaxis()->SetTitle("Dist (mm)");
    h_blobDist_data[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobDist_data[w]->Rebin2D(5,rebinQ2);

    h_nBlobs_data[w] = (TH2D*) f->Get( Form("h_nBlobs_q2qe_%0.2f_%0.1f_data",  w.first, w.second ) );
    h_nBlobs_data[w]->GetXaxis()->SetTitle("N Blobs");
    h_nBlobs_data[w]->GetYaxis()->SetTitle("Events/bin");
    h_nBlobs_data[w]->Rebin2D(1,rebinQ2);

    string comp = "qelike";
    h_blobE_comp[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.2f_%0.1f_%s", w.first, w.second, comp.c_str() ) );
    h_blobE_comp[w]->GetXaxis()->SetTitle("Energy (MeV)");
    h_blobE_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobE_comp[w]->Rebin2D(1, rebinQ2);

    h_blobDist_comp[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.2f_%0.1f_%s",  w.first, w.second, comp.c_str() ) );
    h_blobDist_comp[w]->GetXaxis()->SetTitle("Dist (mm)");
    h_blobDist_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobDist_comp[w]->Rebin2D(5,rebinQ2);

    h_nBlobs_comp[w] = (TH2D*) f->Get( Form("h_nBlobs_q2qe_%0.2f_%0.1f_%s",  w.first, w.second, comp.c_str() ) );
    h_nBlobs_comp[w]->GetXaxis()->SetTitle("N Blobs");
    h_nBlobs_comp[w]->GetYaxis()->SetTitle("Events/bin");
    h_nBlobs_comp[w]->Rebin2D(1,rebinQ2);






    h_blobE[w] = (TH2D*) f->Get( Form("h_blobE_q2qe_%0.2f_%0.1f_%s", w.first, w.second, type.c_str() ) );
    h_blobE[w]->GetXaxis()->SetTitle("Energy (MeV)");
    h_blobE[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobE[w]->Rebin2D(1, rebinQ2);
    h_blobDist[w] = (TH2D*) f->Get( Form("h_blobDist_q2qe_%0.2f_%0.1f_%s", w.first, w.second, type.c_str() ));
    h_blobDist[w]->GetXaxis()->SetTitle("Dist (mm)");
    h_blobDist[w]->GetYaxis()->SetTitle("Events/bin");
    h_blobDist[w]->Rebin2D(5,rebinQ2);
    h_nBlobs[w] = (TH2D*) f->Get( Form("h_nBlobs_q2qe_%0.2f_%0.1f_%s", w.first, w.second, type.c_str() ));
    h_nBlobs[w]->GetXaxis()->SetTitle("# of Blobs");
    h_nBlobs[w]->GetYaxis()->SetTitle("Events/bin");
    h_nBlobs[w]->Rebin2D(1,rebinQ2);



    scale = h_blobE_data[w]->Integral()/h_blobE[w]->Integral();
    h_blobE[w]->Scale(scale);
    h_blobDist[w]->Scale(scale);
    h_nBlobs[w]->Scale(scale);

    h_blobE_comp[w]->Scale(scale);
    h_blobDist_comp[w]->Scale(scale);
    h_nBlobs_comp[w]->Scale(scale);


    h_blobE1D[w] = (TH1D*) h_blobE[w]->ProjectionX();
    h_blobDist1D[w] = (TH1D*) h_blobDist[w]->ProjectionX();
    h_blobDist1D[w]->Rebin(2);
    h_nBlobs1D[w] = (TH1D*) h_nBlobs[w]->ProjectionX();


    h_blobE1D[w]->SetLineColor(i+1);
    h_blobDist1D[w]->SetLineColor(i+1);
    h_nBlobs1D[w]->SetLineColor(i+1);
    cout<<h_blobDist1D[w]<<endl;
  }

/*

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

  int imin = 0;
  int NPlots = h_blobE[1]->GetNbinsY();
  int nX = int(TMath::Sqrt(NPlots))+1;
  cout<<"Dividing canvas: "<<nX<<", "<<nX-1<<endl;
  c->Divide(nX,nX-1,0,0);

  TH2D* h_blobE_base = (TH2D*) h_blobE[DropPars[0]]->Clone("h_blobE_base");
  for( int i = 0; i< n; i++ )
  {
    pair<double,double> w = DropPars[i];
    if(doRatio) h_blobE[w]->Divide( h_blobE_base);
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_blobE[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first,w.second,j+1), j+1,j+1);
      h->SetLineColor(i+1);
      if(doRatio) h->GetYaxis()->SetRangeUser(.5,1.5);
      else h->GetYaxis()->SetRangeUser(0, h_blobE_data[w]->GetMaximum() * 1.5 );
      c->cd(j+1);
      cout<<j+1<<endl;
      h->Draw(opt.c_str() );
    }

    l->AddEntry(h_blobE1D[w], Form("%.2f_%.1f",w.first, w.second ) );
    if ( i==imin ) opt+="same";
  }

  if(doRatio) h_blobE_data[DropPars[0]]->Divide( h_blobE_base );
  if(doRatio) h_blobE_comp[DropPars[0]]->Divide( h_blobE_base );
  for( int j = 0; j< NPlots; j++ )
  {
    c->cd(j+1);
    if(type == "mc")
    {
      TH1D* hdata = (TH1D*) h_blobE_data[DropPars[0]]->ProjectionX(Form("data_%.2f_%.1f_%d",w.first, w.second,j+1), j+1,j+1);
      hdata->SetMarkerStyle(0);
      hdata->Draw("esame");
      //hdata->Draw("esame");
    }
      //TH1D* h_comp = (TH1D*) h_blobE_comp[1.0]->ProjectionX(Form("qelike_%.1f_%d",w,j+1), j+1,j+1);
      //h_comp->SetLineColor(kBlue);
      //h_comp->SetLineStyle(2);
      //h_comp->Draw( "histlsame" );

    double q2Low = h_blobE[1]->GetYaxis()->GetBinLowEdge(j+1);
    double q2High = h_blobE[1]->GetYaxis()->GetBinUpEdge(j+1);
    text->DrawLatexNDC(0.5,0.9, Form("(%.3f, %.3f)",q2Low,q2High));
    //gPad->SetGrid();
    line->Draw();
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  //h_blobE[0]->Draw("colz");
  //c->SetGrid();
  c->Print(Form( "BlobE2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );


  opt="histl";
  TH2D* h_blobDist_base = (TH2D*) h_blobDist[DropPars[0]]->Clone( "h_blobDist_base" );
  for( int i = imin; i< n; i++ )
  {
    pair<double,double> w = DropPars[i];
    if(doRatio) h_blobDist[w]->Divide( h_blobDist_base );
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_blobDist[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first, w.second,j+1), j+1,j+1);
      h->GetXaxis()->SetRangeUser(0,2500);
      h->SetLineColor(i+1);
      if(doRatio) h->GetYaxis()->SetRangeUser(0.5,1.5);
      else h->GetYaxis()->SetRangeUser(0, h_blobDist_data[w]->GetMaximum() * 1.5 );
      c->cd(j+1);
      h->Draw(opt.c_str() );
    }

    if ( i==imin ) opt+="same";
  }


  if(doRatio) h_blobDist_data[DropPars[0]]->Divide( h_blobDist_base);
  for( int j = 0; j< NPlots; j++ )
  {
    c->cd(j+1);
    if(type == "mc")
    {
      TH1D* hdata = (TH1D*) h_blobDist_data[DropPars[0]]->ProjectionX(Form("data_%.2f_%.1f_%d",w.first, w.second,j+1), j+1,j+1);
      hdata->SetMarkerStyle(0);
      hdata->Draw("esame");
      hdata->Draw("esame");
    }

    double q2Low = h_blobE[DropPars[0]]->GetYaxis()->GetBinLowEdge(j+1);
    double q2High = h_blobE[DropPars[0]]->GetYaxis()->GetBinUpEdge(j+1);
    text->DrawLatexNDC(0.5,0.9, Form("(%.3f, %.3f)",q2Low,q2High));
    //gPad->SetGrid();
    line->Draw();
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  //h_blobDist[0]->Draw("colz");
  //c->SetGrid();
  c->Print(Form( "BlobDist2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );

  opt="histl";
  TH2D* h_nBlobs_base = (TH2D*) h_nBlobs[DropPars[0]]->Clone( "h_nBlobs_base" );
  for( int i = imin; i< n; i++ )
  {
    pair<double,double> w = DropPars[i];
    if(doRatio) h_nBlobs[w]->Divide( h_nBlobs_base );
    for( int j = 0; j< NPlots; j++ )
    {
      TH1D* h = (TH1D*) h_nBlobs[w]->ProjectionX(Form("%.2f_%.1f_%d",w.first, w.second,j+1), j+1,j+1);
      h->GetXaxis()->SetRangeUser(0,2500);
      h->SetLineColor(i+1);
      if(doRatio) h->GetYaxis()->SetRangeUser(0.0,1.5);
      else h->GetYaxis()->SetRangeUser(0, h_nBlobs_data[w]->GetMaximum() * 1.5 );
      c->cd(j+1);
      h->Draw(opt.c_str() );
    }

    if ( i==imin ) opt+="same";
  }


  if(doRatio) h_nBlobs_data[DropPars[0]]->Divide( h_nBlobs_base);
  for( int j = 0; j< NPlots; j++ )
  {
    c->cd(j+1);
    if(type == "mc")
    {
      TH1D* hdata = (TH1D*) h_nBlobs_data[DropPars[0]]->ProjectionX(Form("data_%.2f_%.1f_%d",w.first,w.second,j+1), j+1,j+1);
      hdata->SetMarkerStyle(0);
      hdata->Draw("esame");
      hdata->Draw("esame");
    }

    double q2Low = h_blobE[1]->GetYaxis()->GetBinLowEdge(j+1);
    double q2High = h_blobE[1]->GetYaxis()->GetBinUpEdge(j+1);
    text->DrawLatexNDC(0.5,0.9, Form("(%.3f, %.3f)",q2Low,q2High));
    //gPad->SetGrid();
    line->Draw();
  }

  c->cd(NPlots+1);
  l->Draw();
  c->cd(NPlots+2);
  text->DrawLatexNDC(.5,.5,type.c_str() );
  //h_nBlobs[0]->Draw("colz");
  //c->SetGrid();
  
  c->Print(Form( "BlobN2D_%s_%s_ratio-%d.pdf", type.c_str(), tag.c_str(), doRatio) );
*/

}



