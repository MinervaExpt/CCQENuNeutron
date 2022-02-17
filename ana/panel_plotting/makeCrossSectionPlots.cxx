{
  bool doPvalues=false;
  //doPvalues=true;
  //
  bool doRatio = false;
  doRatio = true;


  using namespace PlotUtils;

  vector<string> tunes;
  tunes.push_back("CV2");
  tunes.push_back("CV2Elastic");
  tunes.push_back("CV2ElasticNuwroSF");
  tunes.push_back("CV2ElasticNuwroLFG");

  vector<string> names;
  names.push_back("mnvGENIEv2");
  names.push_back("mnvGENIEv2+Elastic Fix");
  names.push_back("mnvGENIEv2+Elastic Fix+Nuwro SF");
  names.push_back("mnvGENIEv2+Elastic Fix+Nuwro LFG");

  bool pvalues[2]={0,1};
  bool ratio[2]={0,1};

  for( int ip = 0; ip<2; ip++ )
  {
    for( int ir = 0; ir<2;ir++)
    {
      doPvalues = pvalues[ip];
      doRatio = ratio[ir];
      vector<TH1D*> histos;
      TCanvas *c = new TCanvas("c","c");
      c->Divide(2,2);
      TLatex text;
      for( int i = 0; i< tunes.size(); i++ )
      {
        c->cd(i+1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetMargin(.15,.05,.15,.1);
        TFile *f = new TFile(Form("../rootfiles/minervameAntiNu/%s/CrossSectionHist_hydrogen.root", tunes[i].c_str()));
        MnvH1D* mc  = (MnvH1D*) f->Get("h_q2qe_region_00_qelike_qe_h_nobck_hydrogen_unfold_effcor_cross_section");
        MnvH1D* data= (MnvH1D*) f->Get("h_q2qe_region_00_data_nobck_hydrogen_unfold_effcor_cross_section");

        mc->SetLineColor(kRed);

        TH1D* hdataStat= &data->GetCVHistoWithStatError();
        TH1D* hdataTotal=&data->GetCVHistoWithError();
        TH1D* hmc = &mc->GetCVHistoWithStatError();

        if(!doRatio)
        {
          hmc->Scale(1,"width");
          hdataStat->Scale(1,"width");
          hdataTotal->Scale(1,"width");
        }



        hmc->SetMaximum( hmc->GetMaximum()*1.2);
        hmc->GetXaxis()->SetNdivisions(4);
        hmc->GetXaxis()->SetLabelSize(0.05);
        hmc->GetYaxis()->SetLabelSize(0.05);
        hmc->GetXaxis()->SetTitleSize(0.05);
        hmc->GetYaxis()->SetTitleSize(0.05);


        double markerSize=0.05;
        hdataStat->SetMarkerSize(markerSize);
        hdataTotal->SetMarkerSize(markerSize);

        histos.push_back(hdataStat);
        histos.push_back(hdataTotal);
        histos.push_back(hmc);

        hmc->GetXaxis()->SetRangeUser(0.05,5);
        hmc->SetMinimum(0);
        hmc->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
        hmc->GetYaxis()->SetTitle("d#sigma/dQ^{2} (cm^{2}/(GeV/c)^{2}/proton)");
        if(doRatio) 
        {
          hmc->GetYaxis()->SetTitle("Ratio to MC");
          MnvH1D* hmcc = (MnvH1D*) hmc->Clone("hmcc");
          hmc->Divide( hmcc );
          hdataStat->Divide( hmcc );
          hdataTotal->Divide( hmcc );
          hmc->GetYaxis()->SetRangeUser(0,2);
        }
        hmc->Draw("histl");
        hdataStat->Draw("pe1same");
        hdataTotal->Draw("pe1same");
        text.SetTextAngle(0);
        text.SetTextFont(102);
        text.SetTextSize(0.04);
        text.SetTextColor( kBlue );
        text.DrawLatexNDC(0.2,0.85, names[i].c_str() );
        //calculate simple p-value
        if(doPvalues)
        {
          int bin_min = hmc->GetXaxis()->FindBin(0.05);
          int bin_max = hmc->GetXaxis()->FindBin(5);
          cout<<tunes[i]<<"=============="<<endl;
          for ( int ibin = bin_min ; ibin <= bin_max; ibin ++ )
          {
            double center = hmc->GetXaxis()->GetBinCenter(ibin);
            double scale = (doRatio)? 1.: 1.e40;
            double err = hdataTotal->GetBinError( ibin )*scale;
            double vdata = hdataTotal->GetBinContent( ibin )*scale;
            double vmc = hmc->GetBinContent( ibin )*scale;
            double delta = TMath::Abs( vmc - vdata );
            double nSigma = delta/err;

            //cout<<"err: "<<err<<"\t"<<"vdata: "<<vdata<<"\t"<<"vmc: "<<vmc<<"\t"<<"delta: "<<delta<<"\t"<<endl;



            double sf = 1-ROOT::Math::normal_cdf( nSigma );
            cout<<nSigma<<"\t"<<sf<<endl;
            text.SetTextFont(52);
            text.SetTextColor(kRed);
            text.SetTextAngle(45);
            text.SetTextSize(0.05);
            if(!doRatio)text.DrawLatex( center, 18e-39, Form("%.4f", sf) );
            if(doRatio) text.DrawLatex( center, 0.2, Form("%.4f", sf) );
          }
        }

        

      }
      c->Print(Form("plots/minervameAntiNu/CrossSections-ratio_%d-pvlaue_%d.pdf",doRatio, doPvalues));
    }
  }
}
