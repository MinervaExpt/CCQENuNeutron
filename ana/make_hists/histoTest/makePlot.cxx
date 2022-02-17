#include <vector>
#include <TH1D>
const int nbins = 2;
int binslow[nbins] = { 1,5 };
int binshigh[nbins] = { 4,8 };

vector<TH1D*> separateHist( TH2D* h, string name)
{
  vector<TH1D*> hists(2);
  int binlow = h->GetYaxis()->FindBin(0);
  int binmid = h->GetYaxis()->FindBin(0.4);
  int binhigh = h->GetYaxis()->FindBin(1);
  TH1D* hlow = (TH1D*) h->ProjectionX(Form("%s_low",name.c_str()), binlow, binmid );
  TH1D* hhigh = (TH1D*) h->ProjectionX(Form("%s_high",name.c_str()), binmid+1, binhigh );
  cout<<"got hlow/hhigh: "<<hlow<<", "<<hhigh<<endl;
  hists[0]=hlow;
  hists[1]=hhigh;
  cout<<"hists in vector: "<<hists[0]<<", "<<hists[1]<<endl;
  hists[0]->GetYaxis()->SetRangeUser(0, 3500);
  hists[1]->GetYaxis()->SetRangeUser(0, 800);
  cout<<"max set"<<hists[0]->GetMaximum()<<endl;

  return hists;
}

void makePlot()
{
  vector<string> tunes;
  tunes.push_back("CV2");
  tunes.push_back("CV2Elastic");
  tunes.push_back("CV2ElasticNuwroSF");
  tunes.push_back("CV2ElasticNuwroLFG");


  for( int i =  0; i< tunes.size(); i++ )
  {
    cout<<i<<endl;
    TFile *f = new TFile( Form("../../rootfiles/minervameAntiNu/%s/MuonEventSelection_MakeFlux-1_Multiplicity-1_Sample-Signal_CombinedPlaylists.root", tunes[i].c_str() ) );
    TFile *fbck = new TFile( Form("../../rootfiles/minervameAntiNu/%s/signal_weights_minervameNu_q2qe_flux_2D_optimized_Type_weighted_w_lambda_20_1000_noRES.root", tunes[i].c_str() ) );
    TVector2 *pot = (TVector2*) f->Get("pot");
    double scale = pot->X()/pot->Y();
    TH2D* hmc = (TH2D*) f->Get("h_recoil_inc_q2qe_mc");
    TH2D* h_qelike_qe_h = (TH2D*) f->Get("h_recoil_inc_q2qe_qelike_qe_h");
    TH2D* h_qelike_qe_oth = (TH2D*) f->Get("h_recoil_inc_q2qe_qelike_qe_oth");
    TH2D* h_qelike_2p2h = (TH2D*) f->Get("h_recoil_inc_q2qe_qelike_2p2h");
    TH2D* h_qelike_res = (TH2D*) f->Get("h_recoil_inc_q2qe_qelike_res");
    TH2D* h_qelike_dis = (TH2D*) f->Get("h_recoil_inc_q2qe_qelike_dis");
    TH2D* h_qelikenot = (TH2D*) f->Get("h_recoil_inc_q2qe_qelikenot");
    TH2D* h_qelikenot_singlechargedpion = (TH2D*) f->Get("h_recoil_inc_q2qe_qelikenot_singlechargedpion");
    TH2D* h_qelikenot_singleneutralpion = (TH2D*) f->Get("h_recoil_inc_q2qe_qelikenot_singleneutralpion");
    TH2D* h_qelikenot_multipion = (TH2D*) f->Get("h_recoil_inc_q2qe_qelikenot_multipion");
    TH2D* h_qelikenot_nopions = (TH2D*) f->Get("h_recoil_inc_q2qe_qelikenot_nopions");

    TH1D* bck_qelike_qe_oth = (TH1D*) fbck->Get("




    hmc->Scale(scale);
    h_qelike_qe_oth->Scale(scale);

    hmc->SetLineColor(kRed);
    h_qelike_qe_oth->SetLineColor(kGreen);
    h_qelike_qe_h->SetLineColor(kMagenta);
    h_qelikenot->SetLineColor(kGray);




    TH2D* hdata = (TH2D*) f->Get("h_recoil_inc_q2qe_data");
    cout<<"got histos"<<endl;

    vector<TH1D*> hists_mc = separateHist(hmc,"mc");
    vector<TH1D*> hists_data = separateHist(hdata,"data");
    vector<TH1D*> hists_carbon = separateHist(h_qelike_qe_oth,"carbon");
    vector<TH1D*> hists_hydrogen = separateHist(h_qelike_qe_h,"hydrogen");
    vector<TH1D*> hists_qelikenot = separateHist(h_qelikenot,"qelikenot");
    hists_data[0]->GetXaxis()->SetTitle("Recoil (MeV)");
    hists_data[1]->GetXaxis()->SetTitle("Recoil (MeV)");
    cout<<"converted histos"<<endl;

    TCanvas* c = new TCanvas("c1","c1", 800,400);
    c->Divide(2,1);

    c->cd(1);
    hists_data[0]->Draw();
    hists_mc[0]->Draw("histsame");
    hists_carbon[0]->Draw("histlsame");
    hists_hydrogen[0]->Draw("histlsame");
    hists_qelikenot[0]->Draw("histlsame");

    c->cd(2);
    hists_data[1]->Draw();
    hists_mc[1]->Draw("histsame");
    hists_carbon[1]->Draw("histlsame");
    hists_hydrogen[1]->Draw("histlsame");
    hists_qelikenot[1]->Draw("histlsame");

    TLegend* leg = new TLegend(.7,.6,.9,.9);
    leg->AddEntry(hists_data[1],"data");
    leg->AddEntry(hists_mc[1],"mc","l");
    leg->AddEntry(hists_carbon[1],"QELike C","l");
    leg->AddEntry(hists_hydrogen[1],"QELike H","l");
    leg->AddEntry(hists_qelikenot[1],"QELikeNot","l");
    leg->Draw();
    c->Print(Form("recoil_inc_%s.pdf", tunes[i].c_str()));
  }
}
