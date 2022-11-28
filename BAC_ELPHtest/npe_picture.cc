void npe_picture(){

  Int_t F = 1; //# of files
  Int_t y_pos[F];
  y_pos[0] = 0;
  y_pos[1] = -10;
  y_pos[2] = -11;
  y_pos[3] = -23;
  y_pos[4] = -35;
  
  TFile* file[F];
  TTree *tree[F];

  Int_t N;
  Int_t x_position[F][11];
  Double_t mean[F][11];
  Double_t sigma[F][11];
  Double_t simul_mean[F][11];
  Double_t simul_sigma[F][11];

  TGraphErrors *npe[F];
  TGraphErrors *npe_simul[F];
  auto mg_npe = new TMultiGraph();

  for(int i=0;i<F;i++){
    file[i] = new TFile(Form("parameter_%d.root",y_pos[i]),"read");
    tree[i] = (TTree*)file[i]->Get("tree");
    tree[i]->SetBranchAddress("N",&N);
    tree[i]->SetBranchAddress("x_position",x_position[i]);
    tree[i]->SetBranchAddress("mean",mean[i]);
    tree[i]->SetBranchAddress("sigma",sigma[i]);
    tree[i]->SetBranchAddress("simul_mean",simul_mean[i]);
    tree[i]->SetBranchAddress("simul_sigma",simul_sigma[i]);
    tree[i]->GetEntry(0);
    Double_t x_pos[N];
    Double_t pa_mean[N];
    Double_t pa_sigma[N];
    Double_t pa_si_mean[N];
    Double_t pa_si_sigma[N];
    Double_t x_err[N];
    for(int j=0;j<N;j++){
      x_pos[j] = 10.0*x_position[i][j];
      x_err[j] = 0.05;
      pa_mean[j] = mean[i][j];
      pa_sigma[j] = sigma[i][j];
      pa_si_mean[j] = simul_mean[i][j];
      pa_si_sigma[j] = simul_sigma[i][j];
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
    
    mg_npe->Add(npe[i]);
    mg_npe->Add(npe_simul[i]);
    
  }

  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<F;i++){
    le_s->AddEntry(npe[i],"Experiment");
    le_s->AddEntry(npe_simul[i],"Simulation");
  }
  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  mg_npe->SetTitle("Npe;X [mm];N_{p.e.}");
  mg_npe->Draw("AP");
  le_s->Draw();
  


  

  
	      
}
