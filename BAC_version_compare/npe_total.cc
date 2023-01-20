void npe_total(){

  gStyle->SetOptStat(0);
  int N = 5;
  TFile* file[N];
  TH1* hist[N];
  TEfficiency *eff[N];
  for(int i=0;i<N;i++){
    file[i] = new TFile(Form("simul_%d.root",i+4));
    hist[i] = (TH1*)file[i]->Get("hist_total");
    eff[i] = (TEfficiency*)file[i]->Get("eff_total");
    hist[i]->SetMarkerColor(i+1);
    hist[i]->SetLineColor(i+1);
    eff[i]->SetMarkerColor(i+1);
    eff[i]->SetLineColor(i+1);
    eff[i]->SetMarkerStyle(22);
    hist[i]->SetTitle(";N_{p.e.};a.u.");
  }
  TLegend *le_e = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<N;i++){
    le_e->AddEntry(hist[i],Form("version %d",i+1));
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  hist[0]->Draw();
  hist[1]->Draw("same");
  hist[2]->Draw("same");
  hist[3]->Draw("same");
  hist[4]->Draw("same");
  le_e->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  eff[0]->Draw();
  eff[1]->Draw("same");
  eff[2]->Draw("same");
  eff[3]->Draw("same");
  eff[4]->Draw("same");
  le_e->Draw();

  
}
