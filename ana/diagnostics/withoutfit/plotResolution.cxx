{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TChain *CCQENu = new TChain("CCQENu");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/00/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/01/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/02/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/03/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/04/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/05/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/06/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/07/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/08/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/09/*.root");
  CCQENu->Add("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr//minerva/persistent/users/tejinc/CCQENu/Processing_NotUpdateNeutronBlobPos/minervame5A_mc/grid/central_value/minerva/ana/v21r1p1/00/12/30/10/*.root");
  CCQENu->Draw("CCQENu_BlobX-CCQENu_BlobMCTrackX>>h(25,-100,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");

  CCQENu->Draw("CCQENu_BlobX-CCQENu_BlobMCTrackX>>h(25,-100,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");
  c1 = new TCanvas("c1","c1");
  TH1D* hx = (TH1D*) h;
  hx->SetMarkerColor(kBlack);
  hx->SetLineColor(kBlack);
  hx->SetTitle("#delta x");
  hx->GetXaxis()->SetTitle("Distance (mm)");
  hx->Draw("pe");
  hx->Fit("gaus","","",-20,20 );
  c1->Print("dx.pdf");
  
  CCQENu->Draw("CCQENu_BlobY-CCQENu_BlobMCTrackY>>h1(25,-100,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");
  TH1D* hy =(TH1D*)  h1;
  hy->SetMarkerColor(kBlack);
  hy->SetLineColor(kBlack);
  hy->SetTitle("#delta y");
  hy->GetXaxis()->SetTitle("Distance (mm)");
  hy->Draw("pe");
  hy->Fit("gaus","","",-40,40 );
  c1->Print("dy.pdf");
  
  CCQENu->Draw("CCQENu_BlobZ-CCQENu_BlobMCTrackZ>>h2(25,-100,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");
  TH1D* hz =(TH1D*)  h2;
  hz->SetMarkerColor(kBlack);
  hz->SetLineColor(kBlack);
  hz->SetTitle("#delta z");
  hz->GetXaxis()->SetTitle("Distance (mm)");
  hz->Draw("pe");
  hz->Fit("gaus","","",-20,30 );
  c1->Print("dz.pdf");
  
  
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1);
  CCQENu->Draw("((CCQENu_BlobX-CCQENu_BlobMCTrackX)^2+(CCQENu_BlobY-CCQENu_BlobMCTrackY)^2)^.5>>h3(25,0,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");
  TH1D* hxy =(TH1D*)  h3;
  hxy->SetMarkerColor(kBlack);
  hxy->SetLineColor(kBlack);
  hxy->SetTitle("#sqrt{#deltax^{2}+#deltay^{2}}");
  hxy->GetXaxis()->SetTitle("Distance (mm)");
  hxy->Draw("pe");
  //hxy->Fit("gaus","","",-20,30 );
  c1->Print("dxy.pdf");

  gStyle->SetOptFit(1);
  CCQENu->Draw("((CCQENu_BlobX-CCQENu_BlobMCTrackX)^2+(CCQENu_BlobY-CCQENu_BlobMCTrackY)^2+(CCQENu_BlobZ-CCQENu_BlobMCTrackZ)^2)^.5>>h4(10,0,100)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","pe");
  TH1D* hr = (TH1D*) h4;
  hr->SetMarkerColor(kBlack);
  hr->SetLineColor(kBlack);
  hr->SetTitle("#delta R");
  hr->GetXaxis()->SetTitle("Distance (mm)");
  hr->Draw("pe");
  hr->Fit("landau","","",0,100);
  //hr->Fit("gaus","","",-20,30 );
  c1->Print("dr.pdf");

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  CCQENu->Draw("(CCQENu_BlobY-CCQENu_BlobMCTrackY):(CCQENu_BlobX-CCQENu_BlobMCTrackX)>>h2d(10,-50,50,10,-50,50)","mc_intType==1 && CCQENu_BlobIs3D && CCQENu_BlobMCPID==2212","colz");
  TH1D* hxy2d = (TH1D*) h2d;
  hxy2d->SetMarkerColor(kBlack);
  hxy2d->SetLineColor(kBlack);
  hxy2d->SetTitle("y:x");
  hxy2d->GetXaxis()->SetTitle("#deltax (mm)");
  hxy2d->GetYaxis()->SetTitle("#deltay (mm)");
  hxy2d->Draw("colz");
  //hxy->Fit("gaus","","",-20,30 );
  c1->Print("dxy2d.pdf");


}


