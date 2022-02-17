{
  string tune = "CV2ElasticNuwroSF";
  int nW = 10;
  int chi2weights[nW]= {1,5,10,15,20,25,30,35,40,45,50,60,70,100,200,300,500,700,1000,1500,1800,2000,2500,3000,3500,4000,4500,5000};
  //int nW = 11;
  //int chi2weights[nW]={1,3,5,10,15,20,50,100,200,500,1000,};
  int colors[11]={kBlack, kGray, kRed, kPink, kOrange, kYellow, kSpring, kGreen,kTeal, kCyan, kAzure};
  string opt = "histl";
  TCanvas* c = new TCanvas("c","c");
  vector<double>X,Y,X0,Y0;

  TGraph* g = new TGraph();
  TGraph* g0 = new TGraph();
  TGraph* gsum = new TGraph();
  TGraph* gall = new TGraph();

  int I = 0;
  for( int i = 0; i<nW; i++ )
  {
    cout<<i<<endl;
    if(i!=0) opt = "histlsame";
    int w = chi2weights[i];
    TFile *f = new TFile(Form("../../rootfiles/minervameNu/%s/test_signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_%d_%d_noRES.root", tune.c_str(), w,1000 ) , "read" );

    double chi2=0;
    double h0Chi2 = 0;
    for( int h = 0; h<24;h++)
    {

      TH1D* fit = (TH1D*) f->Get(Form("h_q2qe_%02d_mc_fit",h));
      TH1D* orig = (TH1D*)    f->Get(Form("h_q2qe_%02d_mc_orig",h));
      TH1D* data = (TH1D*)  f->Get(Form("h_q2qe_%02d_data",h));

      double dchi2 = 0;
      for( int j=0;j<fit->GetNbinsX();j++ )
      {
        double v= fit->GetBinContent(j+1);
        double d = data->GetBinContent(j+1);
        //cout<<h<<", "<<j<<": "<<d<<endl;
        dchi2+=(d>0)? (v-d)*(v-d)/d:0;
      }
      cout<<dchi2<<", ";
      if( h!=0 ) chi2+= (dchi2/23.);
      if( h == 0 ) h0Chi2 = dchi2;//g0->SetPoint(i, w, dchi2 );
      gall->SetPoint(I,w, dchi2);
      I++;

    }
    cout<<i<<", "<<h0Chi2<<", "<<chi2/23.<<endl;
    g->SetPoint(i,w,chi2);
    g0->SetPoint(i,w,h0Chi2);
    gsum->SetPoint( i,w, chi2+h0Chi2 );
    f->Close();
  }
  cout<<"Drawing Legend"<<endl;
  for(int i = 0; i<X.size();i++) cout<<X[i]<<", "<<Y[i]<<endl;
  TH2D* weight_chi2 = new TH2D("weight_ch2", "weight_chi2",10,0,150,100,0,1000);
  weight_chi2->GetXaxis()->SetTitle("signal weight");
  weight_chi2->GetYaxis()->SetTitle("#chi^{2}, NDF = 18");
  weight_chi2->Draw();


  g->SetFillColor(0);
  g0->SetFillColor(0);
  gsum->SetFillColor(0);
  gall->SetFillColor(0);
  g->SetFillStyle(0);
  g0->SetFillStyle(0);
  gsum->SetFillStyle(0);
  gall->SetFillStyle(0);
  g->SetLineColor(kBlue);
  g0->SetLineColor(kRed);
  gall->SetLineColor(kBlack);
  gsum->SetLineColor(kGray);

  g->Draw("C*");
  g0->Draw("CP*");
  gsum->Draw("CP");
  gall->Draw("P*");

  TLegend *leg= new TLegend(.55,.75,.85,.9 );
  leg->AddEntry(g0,"signal region");
  leg->AddEntry(g,"sidebands average");
  leg->AddEntry(gsum,"sum");
  leg->AddEntry(gall,"individual #chi^{2}");
  leg->Draw();
  c->SetGrid();
  //c->SetLogx();
  //c->SetLogy();
  c->Print( Form("summary/%s_chi2.pdf", tune.c_str(), j ) );
  c->Close();
}
