{
  TFile* f = new TFile("../test_muon_ela_SF_AllBlobs.root","read");
  TH2D* h2p2h = (TH2D*) f->Get("h_n2DBlobs_q2qe_qelike_2p2h");
  TH2D* hmc = (TH2D*) f->Get("h_n2DBlobs_q2qe_mc");
  TH2D* hdata = (TH2D*) f->Get("h_n2DBlobs_q2qe_data");
  TH2D* hqelikeh = (TH2D*) f->Get("h_n2DBlobs_q2qe_qelike_qe_h");
  TH2D* hqelikeoth = (TH2D*) f->Get("h_n2DBlobs_q2qe_qelike_qe_oth");
  TVector2* pot = (TVector2*) f->Get("pot");
  hdata->Scale(pot->Y()/pot->X() );

  int ymax = 18;
  h2p2h1D = (TH1D*) h2p2h->ProjectionY("h2p2h1D",1,ymax);
  hmc1D = (TH1D*) hmc->ProjectionY("hmc1D",1,ymax);
  hdata1D = (TH1D*) hdata->ProjectionY("hdata1D",1,ymax);
  hqelikeh1D = (TH1D*) hqelikeh->ProjectionY("hqelikeh1D",1,ymax);
  hqelikeoth1D = (TH1D*) hqelikeoth->ProjectionY("hqelikeoth1D",1,ymax);

  //h2p2h1D->Divide(hmc1D);
  //hqelikeh1D->Divide(hmc1D);
  //hqelikeoth1D->Divide(hmc1D);

  //hqelikeoth1D->GetYaxis()->SetRangeUser(0,2);
  hdata1D->Draw("p");
  hmc1D->Draw("histsame");
  hqelikeoth1D->Draw("histlsame");
  hqelikeh1D->Draw("histsame");
  h2p2h1D->Draw("histcsame");
  c1->SetLogx();
  c1->Print("test.pdf");

}
