void analysis(){

  //gStyle->SetOptStat(0);

  TFile *file_pe = new TFile("data/run00333.root","read");
  TTree *data_pe = (TTree*)file_pe->Get("tree");
  Double_t ADC_pe[4];
  Double_t ADCs_pe[1];
  data_pe ->SetBranchAddress("E72BACa",ADC_pe);
  data_pe ->SetBranchAddress("E72BACSUMa",ADCs_pe);
  Double_t total_pe = data_pe->GetEntries();

  TH1D* hist_pe[4];
  TH1D* hist_sum = new TH1D("hist_sum","hist_sum",450,50,500);
  TF1* fit_pe_sum = new TF1("fit_sum","gaus(0)",50,500);
  TF1* fit_pe[4];
  for(int j=0;j<4;j++){
    hist_pe[j] = new TH1D(Form("hist_pe_BAC%d",j+1),Form("hist_pe_BAC%d",j+1),450,50,500);
    fit_pe[j] = new TF1(Form("fit_pe_BAC%d",j+1),"gaus(0)",50,500);
  }

  for(int n=0;n<total_pe;n++){
    data_pe->GetEntry(n);
    for(int i=0;i<4;i++){
      hist_pe[i]->Fill(ADC_pe[i]);
    }
    hist_sum->Fill(ADCs_pe[0]);
  }

  Double_t parameter_pe[5][3];

  for(int i=0;i<5;i++){
    if(i<4){
    hist_pe[i]->Fit(fit_pe[i],"Q","",50,500);
    fit_pe[i]->GetParameters(parameter_pe[i]);
    }
    if(i==4){
    hist_sum->Fit(fit_pe_sum,"Q","",50,500);
    fit_pe_sum->GetParameters(parameter_pe[4]);
    }
  }

  //Int_t N=5;
  Int_t N=7;
  TFile *file_po[N];
  TTree *data_po[N];
  Double_t ADCi[N][4];
  Double_t TDCi[N][4][16];
  Double_t ADCs[N][1];
  Double_t TDCs[N][1][16];
  Double_t Ta[N][4][1];
  Double_t Tt[N][4][1][16];

  Double_t tot[N];
  Double_t min = 100000000;


  TH1D *hist_Ta[N][4];
  TH1D *hist_Tt[N][4];

  TH1D *hist_Ta_cut[N][4];
  TH1D *hist_Tt_cut[N][4];
  
  TH1D *hist_inda[N];
  TH1D *hist_indt[N][4];
  TH1D *hist_Suma[N];
  TH1D *hist_Sumt[N];
  TH1D *multi[N];
  TF1 *f_a[N][4];
  TF1 *f_Tt[N][4];
  TF1 *f_inda[N];
  TF1 *f_indt[N][4];
  TF1 *f_sumt[N];
  TF1 *f_suma[N];

  TH2D *ind_sum[N];
  TH2D *adc_tdc_ind[N][4];

  Int_t start_pos = -6;

  file_po[0] = new TFile("data/run00313.root","read");
  file_po[1] = new TFile("data/run00314.root","read");
  file_po[2] = new TFile("data/run00315.root","read");
  file_po[3] = new TFile("data/run00296.root","read");
  file_po[4] = new TFile("data/run00316.root","read");
  file_po[5] = new TFile("data/run00317.root","read");
  file_po[6] = new TFile("data/run00319.root","read");
  
  for(int i=0;i<N;i++){
    data_po[i] = (TTree*)file_po[i]->Get("tree");
    data_po[i] ->SetBranchAddress("E72BACSUMa",ADCs[i]);
    data_po[i] ->SetBranchAddress("E72BACSUMt",TDCs[i]);
    data_po[i] ->SetBranchAddress("E72BACa",ADCi[i]);
    data_po[i] ->SetBranchAddress("E72BACt",TDCi[i]);
    data_po[i] ->SetBranchAddress("T4a",Ta[i][0]);
    data_po[i] ->SetBranchAddress("T4t",Tt[i][0]);
    data_po[i] ->SetBranchAddress("T5a",Ta[i][1]);
    data_po[i] ->SetBranchAddress("T5t",Tt[i][1]);
    data_po[i] ->SetBranchAddress("T6a",Ta[i][2]);
    data_po[i] ->SetBranchAddress("T6t",Tt[i][2]);
    data_po[i] ->SetBranchAddress("T7a",Ta[i][3]);
    data_po[i] ->SetBranchAddress("T7t",Tt[i][3]);
    tot[i] = data_po[i]->GetEntries();
    //if(tot<min)min = tot;
    for(int j=0;j<4;j++){
      hist_Ta[i][j] = new TH1D(Form("hist_T%da_%dcm",j+4,start_pos+2*i),Form("hist_T%da_%dcm",j+4,start_pos+2*i),450,0,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,start_pos+2*i),Form("hist_T%dt_%dcm",j+4,start_pos+2*i),80,600,800);

      hist_Ta_cut[i][j] = new TH1D(Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),450,0,2000);
      hist_Tt_cut[i][j] = new TH1D(Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),80,600,800);
      
      
      hist_indt[i][j] = new TH1D(Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),200,650,850);
      
      f_a[i][j] = new TF1(Form("f_Tt%da_%dcm",j+4,start_pos+2*i),"landau(0)",200,2000);
      f_Tt[i][j] = new TF1(Form("f_Tt%dt_%dcm",j+4,start_pos+2*i),"gaus(0)",650,730);

      f_indt[i][j] = new TF1(Form("f_indt_%dcm_BAC%d",start_pos+2*i,j+1),"gaus(0)",650,850);

      adc_tdc_ind[i][j] = new TH2D(Form("adc_tdc_ind_%dcm_BAC%d",start_pos+2*i,j+1),Form("adc_tdc_ind_%dcm_BAC%d",start_pos+2*i,j+1),100,0,2000,100,670,770);
      
			  
			  
    }
    hist_Suma[i] = new TH1D(Form("hist_Suma_%dcm",start_pos+2*i),Form("hist_Suma_%dcm",start_pos+2*i),100,-10,60);
    hist_Sumt[i] = new TH1D(Form("hist_Sumt_%dcm",start_pos+2*i),Form("hist_Sumt_%dcm",start_pos+2*i),200,650,850);

    hist_inda[i] = new TH1D(Form("hist_inda_%dcm",start_pos+2*i),Form("hist_inda_%dcm",start_pos+2*i),100,-10,60);
    
    f_sumt[i] = new TF1(Form("f_sumt_%dcm",start_pos+2*i),"gaus(0)",650,850);
    f_suma[i] = new TF1(Form("f_suma_%dcm",start_pos+2*i),"landau(0)",-10,60);

    f_inda[i] = new TF1(Form("f_inda_%dcm",start_pos+2*i),"landau(0)",-10,60);
    ind_sum[i] = new TH2D(Form("ind_sum%d",i),Form("ind_sum%d",i),200,0,100,200,0,100);

    multi[i] = new TH1D(Form("multi%d",i),"multiplicity",6,0,6);
    

    
  }


  Double_t entry_Tt[N][4];
  Double_t entry_Ta[N][4];
  Double_t entry_Tt_cut[N][4];
  Double_t entry_Ta_cut[N][4];
  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      for(int j=0;j<4;j++){
	hist_Ta[i][j]->Fill(Ta[i][j][0]);
	
	hist_Tt[i][j]->Fill(Tt[i][j][0][0]);
      }
    }
  }


  Double_t pa_a[N][4][3];
  Double_t pa_t[N][4][3];
  Double_t pa_sumt[N][3];
  Double_t pa_suma[N][3];
  Double_t pa_indt[N][4][3];
  Double_t pa_inda[N][3];

  Int_t pass;
  
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      hist_Ta[i][j]->Fit(f_a[i][j],"Q0","",200,2000);
      f_a[i][j]->GetParameters(pa_a[i][j]);
      std::cout<<"ADC MPV : "<<pa_a[i][j][1]<<" and sigma : "<<pa_a[i][j][2]<<std::endl;
      hist_Tt[i][j]->Fit(f_Tt[i][j],"Q0","",650,730);
      f_Tt[i][j]->GetParameters(pa_t[i][j]);
      std::cout<<"TDC mean : "<<pa_t[i][j][1]<<" and sigma : "<<pa_t[i][j][2]<<std::endl;
    }
  }

  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      pass = 0;
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+50*pa_a[i][j][2]){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-3*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+3*pa_t[i][j][2]){
	    pass+=1;
	    
	  }
	}
	
      }
      if(pass==4){
	for(int j=0;j<4;j++){
	  hist_Ta_cut[i][j]->Fill(Ta[i][j][0]);
	  hist_Tt_cut[i][j]->Fill(Tt[i][j][0][0]);
	}

      }
    }
  }


  
  TCanvas *c1 = new TCanvas("c1","Trigger counters ADC histogram",800,650);
  //c1->Divide(4,N);
  c1->Divide(2,2);
  
  //for(int i=0;i<N;i++){
  for(int i=3;i<4;i++){
    for(int j=0;j<4;j++){
        entry_Ta[i][j] = hist_Ta[i][j]-> GetEntries();
	entry_Tt[i][j] = hist_Tt[i][j]-> GetEntries();
	entry_Ta_cut[i][j] = hist_Ta_cut[i][j]-> GetEntries();
	entry_Tt_cut[i][j] = hist_Ta_cut[i][j]-> GetEntries();
	
      
      //c1->cd(4*i+j+1);
      c1->cd(j+1);
      gPad->SetLogy();
      hist_Ta[i][j]->SetLineColor(kBlack);
      //hist_Ta[i][j]->GetListOfFunctions()->Remove(f_a[i][j]);
      hist_Ta[i][j]->Draw();

      hist_Ta_cut[i][j]->SetLineColor(kRed);
      hist_Ta_cut[i][j]->SetFillColor(kRed);
      hist_Ta_cut[i][j]->SetFillStyle(3001);
      hist_Ta_cut[i][j]->Draw("sames");
      
    }
  }

  TCanvas *c2 = new TCanvas("c2","Trigger counters TDC histogram",800,650);
  //c2->Divide(4,N);
  c2->Divide(2,2);
  //for(int i=0;i<N;i++){
  for(int i=3;i<4;i++){
    for(int j=0;j<4;j++){
      //cc->cd(4*i+j+1);
      c2->cd(j+1);
      gPad->SetLogy();
      hist_Tt[i][j]->SetLineColor(kBlack);
      hist_Tt[i][j]->GetListOfFunctions()->Remove(f_Tt[i][j]);
      hist_Tt[i][j]->Draw();
      hist_Tt_cut[i][j]->SetLineColor(kRed);
      hist_Tt_cut[i][j]->SetFillColor(kRed);
      hist_Tt_cut[i][j]->SetFillStyle(3001);
      hist_Tt_cut[i][j]->Draw("sames");
      
    }
  }
  
  

  /*
  TFile *file3[N];
  TTree *data3[N];
  Double_t ADC_sum[N];
  //gDirectory -> cd("Rint:/");  //Change to Working directory (default)
  
  for(int i=0;i<N;i++){
    file3[i] = new TFile(Form("data/sum_%dcm.root",start_pos+2*i),"recreate");
    data3[i] = new TTree("tree","tree");
    data3[i]->Branch("ADC",&ADC_sum[i],"ADC/D");  //Branch(branchname, &value);
    
  }

  
  

  
  Double_t one_photon= (15.4939+15.753+16.1096+16.0168)*0.945*0.91/4;
  Double_t ind_gain[4];
  ind_gain[0] = 15.4939*0.959;
  ind_gain[1] = 15.753*0.995;
  ind_gain[2] = 16.1096*0.858;
  ind_gain[3] = 16.0168*0.968;
  Double_t numpho;
  Double_t numpho_ind[4];
  Int_t pass_ind;
  Int_t mul_count;


  //BAC TDC histogram
  for(int i=0;i<N;i++){
    file3[i]->cd();
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      pass = 0;
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-1*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+20*pa_a[i][j][2]){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-2*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+2*pa_t[i][j][2]){
	    pass+=1;
	  }
	}
      }
      if(pass==4){
	for(int j=0;j<4;j++)hist_indt[i][j]->Fill(TDCi[i][j][0]);
	hist_Sumt[i]->Fill(TDCs[i][0][0]);
      }
    }
  }


  TCanvas *c3 = new TCanvas("c3","BAC SUM TDC histogram",800,650);
  c3->Divide(N);
  for(int i=0;i<N;i++){
    c3->cd(i+1);
    hist_Sumt[i]->Fit(f_sumt[i],"Q","",650,860);
    f_sumt[i]->GetParameters(pa_sumt[i]);
  }

  TCanvas *c20 = new TCanvas("c20","BAC individual TDC histogram",800,650);
  c20->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c20->cd(4*i+j+1);
      hist_indt[i][j]->Fit(f_indt[i][j],"Q","",700,750);
      f_indt[i][j]->GetParameters(pa_indt[i][j]);
    }
  }

	

  Int_t tot_evt[N];
  Int_t pass_evt[N];
  Int_t pass_mul[N];
  
  
  for(int i=0;i<N;i++){
    file3[i]->cd();
    tot_evt[i] = 0;
    pass_evt[i] = 0;
    pass_mul[i] = 0;
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);

      //Initialize
      pass = 0;
      numpho = 0;
      pass_ind = 0;
      mul_count = 0;
      
	for(int j=0;j<4;j++){
	  numpho_ind[j] = 0;
	  if(Ta[i][j][0]>pa_a[i][j][1]-2*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+20*pa_a[i][j][2]){
	    if(Tt[i][j][0][0]>pa_t[i][j][1]-2*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+2*pa_t[i][j][2]){
	      pass+=1;
	      
	    }
	  }
	}
	if(pass==4){
	  tot_evt[i]+=1;
	  if(TDCs[i][0][0]>pa_sumt[i][1]-3*pa_sumt[i][2]&&TDCs[i][0][0]<pa_sumt[i][1]+3*pa_sumt[i][2]){
	    hist_Suma[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/one_photon);
	    if((ADCs[i][0]-parameter_pe[4][1])/one_photon>0.99)pass_evt[i]+=1;
	    ADC_sum[i] = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	    
	    data3[i]->Fill();

	    for(int j=0;j<4;j++){
	      if(TDCi[i][j][0]>pa_indt[i][j][1]-3*pa_indt[i][j][2]&&TDCi[i][j][0]<pa_indt[i][j][1]+3*pa_indt[i][j][2]){
		numpho += (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
		numpho_ind[j] = (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
	      }
	    }

	    hist_inda[i]->Fill(numpho);
	    ind_sum[i]->Fill(numpho,(ADCs[i][0]-parameter_pe[4][1])/one_photon);
	    for(int j=0;j<4;j++)adc_tdc_ind[i][j]->Fill(ADCi[i][j],TDCi[i][j][0]);

	    //multiplicity
	    for(int j=0;j<4;j++){
	      if(numpho_ind[j]>0.99)mul_count+=1;
	    }
	    multi[i]->Fill(mul_count);
	    if(mul_count>0.99)pass_mul[i]+=1;
	    
	  }
	  
	
	  
	}
	//pass=0;
	  
	
    }
    data3[i]->Write();
    file3[i]->Close();
  }

for(int i=0;i<N;i++){
  std::cout<<"Efficiency(SUM) for position "<<start_pos+i*2<<" cm is "<<1.00*pass_evt[i]/tot_evt[i]<<std::endl;
  std::cout<<"Efficiency(Multiplicity) for position "<<start_pos+i*2<<" cm is "<<1.00*pass_mul[i]/tot_evt[i]<<std::endl;
 }
  
  

  



  TCanvas *c_2d = new TCanvas("c_2d","2D histogram of individual",800,650);
  c_2d->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c_2d->cd(4*i+j+1);
      adc_tdc_ind[i][j]->Draw("colz");
    }
  }
  
  TCanvas *ccom = new TCanvas("ccom","compare ind and sum",800,650);
  ccom->Divide(N);

  for(int i=0;i<N;i++){
    ccom->cd(i+1);
    ind_sum[i]->Draw("colz");
  }

  TCanvas *c_mul = new TCanvas("c_mul","multiplicity",800,650);
  c_mul->Divide(N);
  for(int i=0;i<N;i++){
    c_mul->cd(i+1);
    multi[i]->Draw();
  }



  TCanvas *c4 = new TCanvas("c4","BAC SUM ADC histogram",800,650);
  c4->Divide(N);
  Double_t xpos[N];
  Double_t x_error[N];
  Double_t npe_pos[N];
  Double_t npe_error[N];
  Double_t factor = 1.;
  for(int i=0;i<N;i++){
    c4->cd(i+1);
    //hist_Suma[i]->Scale(factor,"width");
    if(i==0)hist_Suma[i]->Fit(f_suma[i],"Q","",-10,10);
    hist_Suma[i]->Fit(f_suma[i],"Q","",5,100);
    hist_inda[i]->SetLineColor(kBlack);
    hist_inda[i]->Draw("same");
    f_suma[i]->GetParameters(pa_suma[i]);
    
    xpos[i] = start_pos*10+i*20;
    x_error[i] = 0.5;
    npe_pos[i] = pa_suma[i][1];
    npe_error[i] = f_suma[i]->GetParError(1);
  }

  TCanvas *c40 = new TCanvas("c40","BAC indivisual ADC histogram",800,650);
  c40->Divide(N);
  Double_t npe_pos_ind[N];
  Double_t npe_error_ind[N];
  for(int i=0;i<N;i++){
    c40->cd(i+1);
    //hist_Suma[i]->Scale(factor,"width");
    if(i==0)hist_inda[i]->Fit(f_inda[i],"Q","",-10,5);
    hist_inda[i]->Fit(f_inda[i],"Q","",0,100);
    f_inda[i]->GetParameters(pa_inda[i]);
    
    npe_pos_ind[i] = pa_inda[i][1];
    npe_error_ind[i] = f_inda[i]->GetParError(1);
  }

  //lea->Draw();

  //Simulation - position dependence
  TFile *file_simul[N];
  TTree *data_simul[N];
  Int_t simul_npe[N];
  TH1D *hist_simul[N];
  TF1 *fit_simul[N];
  
  file_simul[0] = new TFile("../ELPHsimul/build/elph_221020_m6.5_0_simul.root","read");
  file_simul[1] = new TFile("../ELPHsimul/build/elph_221020_m4.5_0_simul.root","read");
  file_simul[2] = new TFile("../ELPHsimul/build/elph_221020_m2.5_0_simul.root","read");
  file_simul[3] = new TFile("../ELPHsimul/build/elph_221020_m0.5_0_simul.root","read");
  file_simul[4] = new TFile("../ELPHsimul/build/elph_221020_1.5_0_simul.root","read");
  file_simul[5] = new TFile("../ELPHsimul/build/elph_221020_3.5_0_simul.root","read");
  file_simul[6] = new TFile("../ELPHsimul/build/elph_221020_5.5_0_simul.root","read");


  for(int i=0;i<N;i++){
    data_simul[i] = (TTree*)file_simul[i]->Get("tree");
    data_simul[i] ->SetBranchAddress("nhMppc",&simul_npe[i]);
    hist_simul[i] = new TH1D(Form("hist_simul%d",i),Form("hist_simul%d",i),60,0,60);
    fit_simul[i] = new TF1(Form("fit_simul%d",i),"landau(0)",0,60);
  }
    
  for(int n=0;n<14999;n++){
    for(int i=0;i<N;i++){
      data_simul[i]->GetEntry(n);
      hist_simul[i]->Fill(simul_npe[i]);
    }
  }

  Double_t parameter_simul[N][3];
  Double_t npe_simul_pos[N];
  Double_t npe_simul_error[N];

  
  TCanvas *csimul = new TCanvas("csimul","csimul",800,650);
  csimul->Divide(N);
  for(int i=0;i<N;i++){
    csimul->cd(i+1);
    hist_simul[i]->Fit(fit_simul[i],"Q","",0,60);
    fit_simul[i]->GetParameters(parameter_simul[i]);
    npe_simul_pos[i] = parameter_simul[i][1];
    if(i==0||i==N-1)npe_simul_error[i] = 0;
    else{npe_simul_error[i] = fit_simul[i]->GetParError(1);}
  }


  

  TGraphErrors *npe0 = new TGraphErrors(N,xpos,npe_pos,x_error,npe_error);
  TGraphErrors *npe2 = new TGraphErrors(N,xpos,npe_pos_ind,x_error,npe_error_ind);
  
  TGraphErrors *npe1 = new TGraphErrors(N,xpos,npe_simul_pos,x_error,npe_simul_error);
  npe0->SetMarkerStyle(24);
  npe0->SetMarkerColor(1);
  npe0->SetLineColor(1);
  npe0->SetMarkerSize(1);
  npe1->SetMarkerStyle(24);
  npe1->SetMarkerColor(2);
  npe1->SetLineColor(2);
  npe1->SetMarkerSize(1);
  npe2->SetMarkerStyle(24);
  npe2->SetMarkerColor(4);
  npe2->SetLineColor(4);
  npe2->SetMarkerSize(1);

  TMultiGraph *mg = new TMultiGraph();
  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  mg->Add(npe0);
  mg->Add(npe1);
  mg->Add(npe2);
  le_s->AddEntry(npe0,"ELPH test SUM");
  le_s->AddEntry(npe2,"ELPH test Indi.");
  le_s->AddEntry(npe1,"Simulation");
  

  TCanvas *c10 = new TCanvas("c10","c10",800,650);
  c10->cd();
  mg->SetTitle("N_{p.e.} for various X position;X [mm];N_{p.e.}");
  mg->Draw("AP");
  le_s->Draw();


  

  TFile *file_com[2];
  TTree *data_com[2];
  Double_t npe[1];
  Int_t nhMppc[1];
  TH1D *hist_com[2];
  file_com[0] = new TFile("data/sum_0cm.root","read");
  file_com[1] = new TFile("../ELPHsimul/build/elph_221020_0_0_simul.root","read");
  //file_com[1] = new TFile("../../simul_data/221024_ELPH_50_re.root","read");
  //file_com[2] = new TFile("../../BACSimul/build/221024_E72_50.root","read");
  for(int i=0;i<2;i++){
    data_com[i] = (TTree*)file_com[i]->Get("tree");
    hist_com[i] = new TH1D(Form("hist_com%d",i),Form("hist_com%d",i),100,0,100);
  }
  
  data_com[0] ->SetBranchAddress("ADC",&npe[0]);
  data_com[1] ->SetBranchAddress("nhMppc",&nhMppc[0]);


  
  for(int n=0;n<14000;n++){
    for(int i=0;i<2;i++){
      data_com[i]->GetEntry(n);
      
    }
    hist_com[0]->Fill(npe[0]);
    hist_com[1]->Fill(nhMppc[0]);
    hist_com[0]->SetLineColor(kBlack);
    hist_com[1]->SetLineColor(kRed);
    
  }

  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  le->AddEntry(hist_com[0],"ELPH Experiment");
  le->AddEntry(hist_com[1],"Simulation for ELPH");
  

  TCanvas *c5 = new TCanvas("c5","c5",800,650);
  c5->cd();
  hist_com[0]->SetTitle(";N_{p.e.};n");
  hist_com[0]->Draw();
  hist_com[1]->Draw("sames");
  
  le->Draw();

  */

  

  
      
      
      
  
  
  

  
     
}
