{
  string tune = "CV2Elastic";
  int W = 15;

  int nW = 11;
  int nLam = 13;
  int chi2weights[nW]={1,3,5,10,15,20,50,100,200,500,1000,};
  string lambdas[nLam] ={ "0","0.1","0.2","0.5","1","3","5","10","20","50","100","500","1000" };
  float lambdas_v[nLam] ={  0 , 0.1 , 0.2 , 0.5 , 1 , 3 , 5 , 10 , 20 , 50 , 100 , 500 , 1000  };
  int colors[nLam]={kBlack, kGray, kRed, kPink, kOrange, kYellow, kSpring, kGreen,kTeal, kCyan, kAzure, kRed+1, kPink+2};
  string opt = "histl";
  for( int j = 0; j<24;j++)
  {
    TCanvas* c = new TCanvas("c","c");
    TLegend *leg= new TLegend(.1,.1,.9,.3 );
    leg->SetNColumns(3);
    for( int i = 0; i<nLam; i++ )
    {
      if(i>0 && i < 8) continue;
      cout<<i<<endl;
      if(i!=0) opt = "histlsame";
      //int w = chi2weights[i];
      TFile *f = new TFile(Form("../../rootfiles/minervameNu/%s/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_noRes_chi2_w_%d_lambda_%s.root", tune.c_str(), W, lambdas[i].c_str() ) , "read" );
      TH1D* ratio = (TH1D*) f->Get(Form("h_q2qe_%02d_mc_fit",j));
      TH1D* mc = (TH1D*)    f->Get(Form("h_q2qe_%02d_mc_orig",j));
      TH1D* data = (TH1D*)  f->Get(Form("h_q2qe_%02d_data",j));
      //TH1D* mc_ratio = (TH1D*) f->Get("hs_weights_yvar_bgType_qe_oth");
      //TH1D* mc_ratio = (TH1D*) f->Get(Form("h_weights_q2qe_mc_%02d",j));
      TH1D* mc_w_2p2h = (TH1D*) f->Get(Form("hs_weights_yvar_bgType_2p2h"));
      TH1D* mc_w_res = (TH1D*) f->Get(Form("hs_weights_yvar_bgType_res"));
      TH1D* mc_w_qe = (TH1D*) f->Get(Form("hs_weights_yvar_bgType_qe_oth"));
      cout<<"mc: "<<mc->Integral()<<endl;;
      cout<<"fit: "<<ratio->Integral()<<endl;;
      cout<<"data:"<<data->Integral()<<endl;;
      ratio->Divide(data);
      ratio->SetLineColor( colors[i] );
      ratio->GetYaxis()->SetRangeUser(-4,5);
      ratio->Draw( opt.c_str() );
      mc_w_2p2h->SetLineColor( colors[i] );
      mc_w_2p2h->SetLineStyle( 2 );
      mc_w_2p2h->SetLineWidth( 1 );


      mc_w_res->SetLineColor( colors[i] );
      mc_w_res->SetLineStyle( 6 );
      mc_w_res->SetLineWidth( 1 );

      mc_w_qe->SetLineColor( colors[i] );
      mc_w_qe->SetLineStyle( 1 );
      mc_w_qe->SetLineWidth( 1 );

      mc_w_qe->Draw("histlsame");
      mc_w_2p2h->Draw("histlsame");
      mc_w_res->Draw("histlsame");

      //leg->AddEntry( ratio, Form("w = %d", w) );
      leg->AddEntry( ratio, Form("w = %d, l = %.1f", W, lambdas_v[i]) );
      f->Close();
    }
    cout<<"Drawing Legend"<<endl;
    leg->Draw();
    c->SetGrid();
    c->SetLogx();
    c->Print( Form("spectrum/%s_w_%d_lambda_%02d.pdf", tune.c_str(),W, j ) );
    c->Close();
  }
}
