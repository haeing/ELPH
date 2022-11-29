void npe_picture(){

  Int_t F = 3; //# of files
  Int_t y_pos[F];
  //y_pos[0] = 0;
  //y_pos[1] = -10;
  y_pos[0] = -11;
  y_pos[1] = -23;
  y_pos[2] = -35;
  
  TFile* file[F];
  TTree *tree[F];

  Int_t N;
  Int_t x_position[F][11];
  Double_t mean[F][11];
  Double_t sigma[F][11];
  Double_t simul_mean[F][11];
  Double_t simul_sigma[F][11];
  Double_t efficiency[F][11];
  Double_t effi_error[F][11];

  TGraphErrors *npe[F];
  TGraphErrors *npe_simul[F];
  TGraphErrors *eff[F];
  auto mg_npe = new TMultiGraph();
  auto mg_eff = new TMultiGraph();

  for(int i=0;i<F;i++){
    file[i] = new TFile(Form("parameter_%d.root",y_pos[i]),"read");
    tree[i] = (TTree*)file[i]->Get("tree");
    tree[i]->SetBranchAddress("N",&N);
    tree[i]->SetBranchAddress("x_position",x_position[i]);
    tree[i]->SetBranchAddress("mean",mean[i]);
    tree[i]->SetBranchAddress("sigma",sigma[i]);
    tree[i]->SetBranchAddress("simul_mean",simul_mean[i]);
    tree[i]->SetBranchAddress("simul_sigma",simul_sigma[i]);
    tree[i]->SetBranchAddress("efficiency",efficiency[i]);
    tree[i]->SetBranchAddress("effi_error",effi_error[i]);
    tree[i]->GetEntry(0);
    Double_t x_pos[N];
    Double_t pa_mean[N];
    Double_t pa_sigma[N];
    Double_t pa_si_mean[N];
    Double_t pa_si_sigma[N];
    Double_t pa_effi[N];
    Double_t pa_effi_err[N];
    Double_t x_err[N];
    for(int j=0;j<N;j++){
      x_pos[j] = 10.0*x_position[i][j];
      x_err[j] = 0.05;
      pa_mean[j] = mean[i][j];
      pa_sigma[j] = sigma[i][j];
      pa_si_mean[j] = simul_mean[i][j];
      pa_si_sigma[j] = simul_sigma[i][j];
      pa_effi[j] = efficiency[i][j];
      pa_effi_err[j] = effi_error[i][j];
    }
    npe[i] = new TGraphErrors(N,x_pos,pa_mean,x_err,pa_sigma);
    npe[i]->SetMarkerStyle(24);
    npe[i]->SetMarkerSize(2);
    npe[i]->SetMarkerColor(kBlack);
    npe[i]->SetLineColor(kBlack);

    npe_simul[i] = new TGraphErrors(N,x_pos,pa_si_mean,x_err,pa_si_sigma);
    npe_simul[i]->SetMarkerStyle(22);
    npe_simul[i]->SetMarkerSize(2);
    npe_simul[i]->SetMarkerColor(kBlue);
    npe_simul[i]->SetLineColor(kBlue);

    eff[i]= new TGraphErrors(N,x_pos,pa_effi,x_err,pa_effi_err);
    eff[i]->SetMarkerStyle(24);
    eff[i]->SetMarkerSize(1.5);
    
    
    
    
    mg_npe->Add(npe[i]);
    mg_npe->Add(npe_simul[i]);

    mg_eff->Add(eff[i]);
    
  }

  eff[0]->SetMarkerColor(1);
  eff[1]->SetMarkerColor(2);
  eff[2]->SetMarkerColor(4);
  //eff[3]->SetMarkerColor(8);
  //eff[4]->SetMarkerColor(28);
  
  eff[0]->SetLineColor(1);
  eff[1]->SetLineColor(2);
  eff[2]->SetLineColor(4);
  //eff[3]->SetLineColor(8);
  //eff[4]->SetLineColor(28);

  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<F;i++){
    le_s->AddEntry(npe[i],"Experiment");
    le_s->AddEntry(npe_simul[i],"Simulation");
  }

  TLegend *le_e = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<F;i++){
    le_e->AddEntry(eff[i],Form("y = %d mm",y_pos[i]));
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  mg_npe->SetTitle("Npe;X [mm];N_{p.e.}");
  mg_npe->Draw("AP");
  le_s->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  c2->cd();
  mg_eff->SetTitle("Efficiency;X [mm];Efficiency");
  mg_eff->Draw("AP");
  le_e->Draw();
  
    
  


  

  
	      
}