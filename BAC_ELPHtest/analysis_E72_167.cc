void analysis_E72_167(){
  Int_t X = 6;
  Int_t Y = 6;
  Int_t M = 5;
  Int_t x_pos[X];
  Int_t y_pos[Y];
  Int_t mom[M];
  for(int i=0;i<X;i++)x_pos[i] = 0+10*i;
  for(int i=0;i<Y;i++)y_pos[i] = 0+10*i;

  for(int i=0;i<M;i++)mom[i] = 700+20*i;

  TFile *file[X][Y][M];
  TTree *tree[X][Y][M];

  Int_t nhMppc[X][Y][M];
  Int_t mppcnum[X][Y][M][1000];
  TH1D *hist_simul[X][Y][M];

  TH1D* hist_total = new TH1D("hist_total","hist_total",40,0,40);
  TF1 *fit_simul[X][Y][M];

  Double_t pa_simul[X][Y][M][3];
  Double_t pa_error[X][Y][M];

  TGraphErrors *npe[M][X];

  TCanvas *c1[M];
  Int_t evt[X][Y][M];

  TEfficiency *eff[M][X];
  Bool_t passed;

  
  for(int i=0;i<M;i++){
    c1[i] = new TCanvas(Form("c1%d",i),"c1",800,650);
    c1[i]->Divide(Y,X);
  }

  for(int i=0;i<M;i++){
    for(int j=0;j<X;j++){
      //eff[i][j] = new TEfficiency(Form("eff%d%d",i,j),Form("mppc = %d mm;Y [mm];Efficiency",mom[i]),110,-60,60);
      eff[i][j] = new TEfficiency(Form("eff%d%d",i,j),Form("p_{#pi} = %d MeV/c;Y [mm];Efficiency",mom[i]),110,-60,60);
      eff[i][j]->SetMarkerStyle(22);
      eff[i][j]->SetMarkerSize(1.5);
      eff[i][j]->SetMarkerColor(j+1);
      eff[i][j]->SetLineColor(j+1);
    }
  }
  
  for(int i=0;i<X;i++){
    for(int j=0;j<Y;j++){
      for(int k=0;k<M;k++){
	
	file[i][j][k] = new TFile(Form("~/E72/ELPH_data/simul_167_36/bac_%dmm_%dmm_mom_%d_mppc_36mm.root",x_pos[i],y_pos[j],mom[k]),"read");
	tree[i][j][k] = (TTree*)file[i][j][k]->Get("tree");
	tree[i][j][k]->SetBranchAddress("nhMppc",&nhMppc[i][j][k]);
	tree[i][j][k]->SetBranchAddress("mppcnum",mppcnum[i][j][k]);
	hist_simul[i][j][k] = new TH1D(Form("simul%d%d%d",x_pos[i],y_pos[j],mom[k]),Form("simul%d%d%d",x_pos[i],y_pos[j],mom[k]),50,0,50);
	fit_simul[i][j][k] = new TF1(Form("f_simul%d%d%d",x_pos[i],y_pos[j],mom[k]),"gaus(0)",0,100);
	evt[i][j][k]=0;
	passed = 0;
	
	for(int n=0;n<1000;n++){
	  tree[i][j][k]->GetEntry(n);
	  hist_simul[i][j][k]->Fill(nhMppc[i][j][k]);
	  hist_total->Fill(nhMppc[i][j][k]);

	  if(nhMppc[i][j][k]>7)passed =1;

	  else{passed = 0;}
	  eff[k][i]->Fill(passed,y_pos[j]);
	  
	  
	}
	c1[k]->cd(Y*i+j+1);
	hist_simul[i][j][k]->Fit(fit_simul[i][j][k],"","",0,50);
	fit_simul[i][j][k]->GetParameters(pa_simul[i][j][k]);
	pa_error[i][j][k] = fit_simul[i][j][k]->GetParError(1);

	
      }
    }
  }

  TMultiGraph *mg[M];
  TLegend *le[M];
  TLegend *le_e[M];
  
  for(int i=0;i<M;i++){
    mg[i] = new TMultiGraph();
    le[i] =  new TLegend(0.8,0.5,0.48,0.6);
    le_e[i] =  new TLegend(0.8,0.5,0.48,0.6);
    
  }
  
  Double_t phonum[Y];
  Double_t pho_err[Y];
  Double_t yposition[Y];
  Double_t yposition_err[Y];
  TCanvas *c2[M];
  TCanvas *c3[M];
  
  
  for(int i=0;i<M;i++){
    for(int j=0;j<X;j++){
      for(int k=0;k<Y;k++){
	yposition[k] = y_pos[k]*1.0;
	yposition_err[k] = 0.5;
	phonum[k] = pa_simul[j][k][i][1];
	pho_err[k] = pa_error[j][k][i];
      }
      npe[i][j] = new TGraphErrors(Y,yposition,phonum,yposition_err,pho_err);
      npe[i][j]->SetMarkerStyle(24);
      npe[i][j]->SetMarkerSize(2);
      npe[i][j]->SetMarkerColor(j+1);
      npe[i][j]->SetLineColor(j+1);
      mg[i]->Add(npe[i][j]);
      le[i]->AddEntry(npe[i][j],Form("X = %d mm",x_pos[j]));
      le_e[i]->AddEntry(eff[i][j],Form("X = %d mm",x_pos[j]));
      
      
    }
    mg[i]->SetTitle(Form("p_{#pi} = %d MeV/c;Y [mm];N_{p.e.}",mom[i]));
    //mg[i]->SetTitle(Form("mppc = %d mm;Y [mm];N_{p.e.}",mom[i]));
    c2[i] = new TCanvas(Form("c2%d",i),"c2",800,650);
    
    c2[i]->cd();
    mg[i]->Draw("AP");
    le[i]->Draw();

    c3[i] = new TCanvas(Form("c3%d",i),"c3",800,650);
    c3[i]->cd();
    for(int j=0;j<X;j++){
      if(j==0)eff[i][j]->Draw("");
      else{eff[i][j]->Draw("same");}
      
    }
    le_e[i]->Draw();
      
    
  }

  TCanvas *c_total = new TCanvas("c_total","c_total",800,650);
  hist_total->Draw();
    
    
      

  

	
  
}
