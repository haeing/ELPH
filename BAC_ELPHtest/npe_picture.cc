void npe_picture(){

  Int_t F = 1; //# of files
  Int_t y_pos[F];
  y_pos[0] = -10;
  //y_pos[1] = -10;
  //y_pos[0] = -11;
  //y_pos[1] = -23;
  //y_pos[2] = -35;
  
  TFile* file[F];
  TTree *tree[F];

  Int_t N;
  Int_t x_position[F][11];
  Double_t mean[F][11];
  Double_t sigma[F][11];
  Double_t error[F][11];
  Double_t simul_mean[F][11];
  Double_t simul_sigma[F][11];
  Double_t simul_error[F][11];
  Int_t total_event[F][11];
  Int_t selected_event[F][11];

  TGraphErrors *npe[F];
  TGraphErrors *npe_simul[F];
  auto mg_npe = new TMultiGraph();
    
  auto mg_eff = new TMultiGraph();

  TEfficiency* eff[F];

  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  //c2->cd();
  //eff[0]->SetTitle("Efficiency;X [mm];Efficiency");
  //eff[0]->Draw("");
  //eff[1]->Draw("same");
  //eff[2]->Draw("same");
  //le_e->Draw();

  TLegend *le_e = new TLegend(0.8,0.5,0.48,0.6);

  TCanvas *c3[F];


  for(int i=0;i<F;i++){
    //eff[i] = new TEfficiency(Form("eff%d",i),"Efficiency;X [mm];Efficiency",100,-60,70);
    file[i] = new TFile(Form("parameter_%d.root",y_pos[i]),"read");
    eff[i] = (TEfficiency*)file[i]->Get("eff");
    tree[i] = (TTree*)file[i]->Get("tree");
    tree[i]->SetBranchAddress("N",&N);
    tree[i]->SetBranchAddress("x_position",x_position[i]);
    tree[i]->SetBranchAddress("mean",mean[i]);
    tree[i]->SetBranchAddress("sigma",sigma[i]);
    tree[i]->SetBranchAddress("error",error[i]);
    tree[i]->SetBranchAddress("simul_mean",simul_mean[i]);
    tree[i]->SetBranchAddress("simul_sigma",simul_sigma[i]);
    tree[i]->SetBranchAddress("simul_error",simul_error[i]);
    tree[i]->SetBranchAddress("total_event",total_event[i]);
    tree[i]->SetBranchAddress("selected_event",selected_event[i]);
    tree[i]->GetEntry(0);
    Double_t x_pos[N];
    Double_t pa_mean[N];
    Double_t pa_sigma[N];
    Double_t pa_error[N];
    Double_t pa_si_mean[N];
    Double_t pa_si_sigma[N];
    Double_t pa_si_error[N];
    Double_t pa_tot[N];
    Double_t pa_selected[N];
    
    Double_t x_err[N];
    for(int j=0;j<N;j++){
      x_pos[j] = 10.0*x_position[i][j];
      x_err[j] = 0.05;
      pa_mean[j] = mean[i][j];
      pa_sigma[j] = sigma[i][j];
      pa_error[j] = error[i][j];
      pa_si_mean[j] = simul_mean[i][j];
      pa_si_sigma[j] = simul_sigma[i][j];
      pa_si_error[j] = simul_error[i][j];
      pa_tot[j] = total_event[i][j];
      pa_selected[j] = selected_event[i][j];
    }
    npe[i] = new TGraphErrors(N,x_pos,pa_mean,x_err,pa_error);
    npe[i]->SetMarkerStyle(24);
    npe[i]->SetMarkerSize(2);
    npe[i]->SetMarkerColor(kBlack);
    npe[i]->SetLineColor(kBlack);

    npe_simul[i] = new TGraphErrors(N,x_pos,pa_si_mean,x_err,pa_si_error);
    npe_simul[i]->SetMarkerStyle(22);
    npe_simul[i]->SetMarkerSize(2);
    npe_simul[i]->SetMarkerColor(kBlue);
    npe_simul[i]->SetLineColor(kBlue);

    //eff[i]= new TEfficiency("eff","Efficiency;X [mm];Efficiency",10,0,1);
    eff[i]->SetMarkerStyle(24);
    eff[i]->SetMarkerColor(i+1);
    eff[i]->SetLineColor(i+1);
    eff[i]->SetMarkerSize(1.5);
    le_e->AddEntry(eff[i],Form("y = %d mm",y_pos[i]));
    c2->cd();
    if(i==0)eff[i]->Draw("AP");
    else if(i>0)eff[i]->Draw("same");

    c3[i] = new TCanvas("c3","c3",800,650);
    c3[i]->cd();
    eff[i]->Draw("AP");
    
    
    
    
    
    mg_npe->Add(npe[i]);
    mg_npe->Add(npe_simul[i]);

    //mg_eff->Add(eff[i]);
    
  }

  c2->cd();
  le_e->Draw("same");

  /*
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
  */


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
