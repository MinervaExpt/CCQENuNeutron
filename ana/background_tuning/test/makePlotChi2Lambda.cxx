{
  using namespace PlotUtils;
  string tune = "CV2ElasticNuwroSF";
  int nW = 11;
  int nLam = 15;
  //string lambdas[nLam] ={ "0","0.1","0.2","0.5","1","3","5","10","20","50","100","500","1000" };
  //float lambdas_v[nLam] ={  0 , 0.1 , 0.2 , 0.5 , 1 , 3 , 5 , 10 , 20 , 50 , 100 , 500 , 1000  };

  //string lambdas[nLam] ={ "0","0.001","0.005","0.01","0.05","0.1","0.2", "0.5","0.7",  "1","2","3","5", "7", "9"};
  //float lambdas_v[nLam] ={0,0.001,0.005,0.01,0.05,0.1,0.2, 0.5,0.7,  1,2,3,5, 7, 9};
  //nLam = 15;
  //string lambdas[nLam] ={ "0","1","5","10","20","30","50","100","200","500","1000","2000","5000","10000","20000"};
  //float lambdas_v[nLam]={ 0,1,5,10,20,30,50,100,200,500,1000, 2000, 5000, 10000, 20000};
  //
  nLam = 8;
  string lambdas[nLam] ={ "0","1","10","100","1000","10000", "50000","100000"};
  float lambdas_v[nLam]={ 0,1,10,100,1000,10000,50000, 100000};
  int rw = 30;

  


  int chi2weights[nW]={1,3,5,10,15,20,50,100,200,500,1000,};
  int colors[nW]={kBlack, kGray, kRed, kPink, kOrange, kYellow, kSpring, kGreen,kTeal, kCyan, kAzure};
  string opt = "histl";
  TCanvas* c = new TCanvas("c","c");
  vector<double>X,Y,X0,Y0;

  TGraph* g = new TGraph();
  TGraph* g0 = new TGraph();
  TGraph* gsum = new TGraph();
  TGraph* gall = new TGraph();

  int I = 0;
  for( int i = 0; i<nLam; i++ )
  {
    cout<<i<<endl;
    if(i!=0) opt = "histlsame";
    string w = lambdas[i];
    TFile *f = new TFile(Form("../../rootfiles/minervameNu/%s/test_signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_%d_%s_noRES.root", tune.c_str(), rw,w.c_str() ) , "read" );
    cout<< "file is open? "<<f->IsOpen()<<endl;


    //TH1D* h = (TH1D*) f-Get(
    MnvH1D*h0 = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qe_oth");
    MnvH1D*h1 = (MnvH1D*) f->Get("hs_weights_yvar_bgType_2p2h");
    MnvH1D*h2 = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_scp");
    MnvH1D*h3 = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_snp");
    MnvH1D*h4 = (MnvH1D*) f->Get("hs_weights_yvar_bgType_qelikenot_mp");

    cout<<"Integral: "<<h0->Integral()<<endl;
    
    vector<TH1D> hists;
    hists.push_back(  h0->GetCVHistoWithError() );
    hists.push_back(  h1->GetCVHistoWithError() );
    hists.push_back(  h2->GetCVHistoWithError() );
    hists.push_back(  h3->GetCVHistoWithError() );
    hists.push_back(  h4->GetCVHistoWithError() );
    
    cout<<"created hists holder"<<endl;
    cout<<hists[0].GetName()<<endl;


    double reg = 0;
    double lam = lambdas_v[i];

    cout<<"==================lambda = "<<lam<<endl;
    if(lam>0)
    {
      for(int h = 0; h<hists.size();h++)
      {
        cout<<"at hist "<<h<<" : ";
        for(int l = 0; l<hists[0].GetNbinsX()-2;l++)
        {
          cout<<"bin "<<l+1<<endl;
          double v0 = hists[h].GetBinContent(l+1);
          double v1 = hists[h].GetBinContent(l+2);
          double v2 = hists[h].GetBinContent(l+3);
          double dreg = (v0+v2-2*v1);
          dreg*=dreg;
          cout<<"( "<<v0<<" "<<v1<<" "<<v2<<" )";
          reg+= dreg;
        }
      }
    }
    cout<<"reg: "<<reg<<", lam: "<<lam<<endl;
    reg*=lam;

    reg=lam;

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
      gall->SetPoint(I,reg, dchi2);
      I++;
    }

    cout<<i<<", "<<h0Chi2<<", "<<chi2/23.<<", "<<reg<<endl;
    cout<<"reg: "<<reg<<endl;
    g->SetPoint(i,reg,chi2);
    g0->SetPoint(i,reg,h0Chi2);
    gsum->SetPoint( i,reg, chi2+h0Chi2 );
    f->Close();
    cout<<"File Closed"<<endl;
  }
  cout<<"Drawing Legend"<<endl;
  for(int i = 0; i<X.size();i++) cout<<X[i]<<", "<<Y[i]<<endl;
  //TH2D* weight_chi2 = new TH2D("weight_ch2", "weight_chi2",10,0.001,15,10,0,.15);
  TH2D* weight_chi2 = new TH2D("weight_ch2", "weight_chi2",10,0.001,120000,10,0,1200);
  weight_chi2->GetXaxis()->SetTitle("reg");
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
  c->SetLogx();
  //c->SetLogy();
  c->Print( Form("summary/%s_chi2_lam_regularized_w%d.pdf", tune.c_str(),rw, j ) );
  c->Close();
}
