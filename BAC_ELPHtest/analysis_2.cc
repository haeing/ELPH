void analysis_2(){

  //Pedestal
  TFile *file_pe = new TFile("../../ELPH_data/exp_data/run00333.root","read");
  //TFile *file_pe = new TFile("../../ELPH_data/exp_data/run00079.root","read");
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

  Double_t pe_sum;
  
  for(int n=0;n<total_pe;n++){
    data_pe->GetEntry(n);
    pe_sum = 0;
    for(int i=0;i<4;i++){
      hist_pe[i]->Fill(ADC_pe[i]);
      pe_sum+=ADC_pe[i];
    }
    hist_sum->Fill(ADCs_pe[0]);
  }

  Double_t parameter_pe[5][3];

  for(int i=0;i<5;i++){
    if(i<4){
    hist_pe[i]->Fit(fit_pe[i],"","",50,500);
    fit_pe[i]->GetParameters(parameter_pe[i]);
    }
    if(i==4){
    hist_sum->Fit(fit_pe_sum,"","",50,500);
    fit_pe_sum->GetParameters(parameter_pe[4]);
    }
  }

  Int_t N=7;
  Int_t start_pos = -6;
  TFile *file_po[N];
  TTree *data_po[N];
  Double_t ADCi[N][4];
  Double_t TDCi[N][4][16];
  Double_t ADCs[N][1];
  Double_t TDCs[N][1][16];
  Double_t Ta[N][4][1];
  Double_t Tt[N][4][1][16];

  Double_t tot[N];

  TH1D *hist_Ta[N][4];
  TH1D *hist_Tt[N][4];

  TH1D *hist_Ta_cut[N][4];
  TH1D *hist_Tt_cut[N][4];


  TH1D *hist_inda[N];
  TH1D *hist_suma[N];

  TH1D *hist_inda_cut[N];
  TH1D *hist_suma_cut[N];

  TH1D *hist_indt[N][4];
  TH1D *hist_sumt[N];

  TH1D *hist_indt_cut[N][4];
  TH1D *hist_sumt_cut[N];

  TH2D *ind_sum[N];

  TF1 *f_a[N][4];
  TF1 *f_t[N][4];
  
  TF1 *f_inda[N];
  TF1 *f_suma[N];

  TF1 *f_indt[N][4];
  TF1 *f_sumt[N];

  //y=0 mm
  file_po[0] = new TFile("../../ELPH_data/exp_data/run00313.root","read");
  file_po[1] = new TFile("../../ELPH_data/exp_data/run00314.root","read");
  file_po[2] = new TFile("../../ELPH_data/exp_data/run00315.root","read");
  file_po[3] = new TFile("../../ELPH_data/exp_data/run00297.root","read");
  file_po[4] = new TFile("../../ELPH_data/exp_data/run00316.root","read");
  file_po[5] = new TFile("../../ELPH_data/exp_data/run00317.root","read");
  file_po[6] = new TFile("../../ELPH_data/exp_data/run00319.root","read");


  //y=-23 mm
  /*
  file_po[0] = new TFile("../../ELPH_data/exp_data/run00048.root","read");
  file_po[1] = new TFile("../../ELPH_data/exp_data/run00052.root","read");
  file_po[2] = new TFile("../../ELPH_data/exp_data/run00053.root","read");
  file_po[3] = new TFile("../../ELPH_data/exp_data/run00045.root","read");
  file_po[4] = new TFile("../../ELPH_data/exp_data/run00054.root","read");
  file_po[5] = new TFile("../../ELPH_data/exp_data/run00055.root","read");
  file_po[6] = new TFile("../../ELPH_data/exp_data/run00056.root","read");
  */

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

    for(int j=0;j<4;j++){
      hist_Ta[i][j] = new TH1D(Form("hist_T%da_%dcm",j+4,start_pos+2*i),Form("hist_T%da_%dcm",j+4,start_pos+2*i),450,0,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,start_pos+2*i),Form("hist_T%dt_%dcm",j+4,start_pos+2*i),250,600,850);

      hist_Ta_cut[i][j] = new TH1D(Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),450,0,2000);
      hist_Tt_cut[i][j] = new TH1D(Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),250,600,850);

      hist_indt[i][j] = new TH1D(Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),500,0,1500);

      hist_indt_cut[i][j] = new TH1D(Form("hist_indt_cut_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_indt_cut_%dcm_BAC%d",start_pos+2*i,j+1),500,0,1500);
      
      f_a[i][j] = new TF1(Form("f_a%da_%dcm",j+4,start_pos+2*i),"landau(0)",200,2000);
      f_t[i][j] = new TF1(Form("f_t%dt_%dcm",j+4,start_pos+2*i),"gaus(0)",650,730);

      f_indt[i][j] = new TF1(Form("f_indt_%dcm_BAC%d",start_pos+2*i,j+1),"gaus(0)",650,850);
    }
    
    hist_inda[i] = new TH1D(Form("hist_inda_%dcm",start_pos+2*i),Form("hist_inda_%dcm",start_pos+2*i),100,-10,60);
    hist_suma[i] = new TH1D(Form("hist_Suma_%dcm",start_pos+2*i),Form("hist_Suma_%dcm",start_pos+2*i),100,-10,60);
    hist_sumt[i] = new TH1D(Form("hist_Sumt_%dcm",start_pos+2*i),Form("hist_Sumt_%dcm",start_pos+2*i),500,0,1500);

    hist_inda_cut[i] = new TH1D(Form("hist_inda_cut_%dcm",start_pos+2*i),Form("hist_inda_cut_%dcm",start_pos+2*i),100,-10,60);
    hist_suma_cut[i] = new TH1D(Form("hist_suma_cut_%dcm",start_pos+2*i),Form("hist_suma_cut_%dcm",start_pos+2*i),100,-10,60);
    hist_sumt_cut[i] = new TH1D(Form("hist_sumt_cut_%dcm",start_pos+2*i),Form("hist_sumt_cut_%dcm",start_pos+2*i),500,0,1500);

    ind_sum[i] = new TH2D(Form("ind_sum%d",i),Form("ind_sum%d",i),260,-100,2500,160,-100,1500);
    

    f_inda[i] = new TF1(Form("f_inda_%dcm",start_pos+2*i),"gaus(0)",-10,60);
    f_sumt[i] = new TF1(Form("f_sumt_%dcm",start_pos+2*i),"gaus(0)",650,850);
    f_suma[i] = new TF1(Form("f_suma_%dcm",start_pos+2*i),"gaus(0)",-10,60);
  }

  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      for(int j=0;j<4;j++){
	hist_Ta[i][j]->Fill(Ta[i][j][0]);
	
	hist_Tt[i][j]->Fill(Tt[i][j][0][0]);
      }
    }
  }

  //Trigger counters Cut condition

  Double_t pa_a[N][4][3];
  Double_t pa_t[N][4][3];
  
  Int_t pass;
  
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      hist_Ta[i][j]->Fit(f_a[i][j],"Q0","",200,2000);
      f_a[i][j]->GetParameters(pa_a[i][j]);
      hist_Tt[i][j]->Fit(f_t[i][j],"Q0","",650,730);
      f_t[i][j]->GetParameters(pa_t[i][j]);
    }
  }

  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      pass = 0;
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
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
  c1->Divide(4,N);
  
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c1->cd(4*i+j+1);
      gPad->SetLogy();
      hist_Ta[i][j]->SetLineColor(kBlack);
      hist_Ta[i][j]->Draw();

      hist_Ta_cut[i][j]->SetLineColor(kRed);
      hist_Ta_cut[i][j]->SetFillColor(kRed);
      hist_Ta_cut[i][j]->SetFillStyle(3001);
      hist_Ta_cut[i][j]->Draw("sames");
      
    }
  }

  TCanvas *c2 = new TCanvas("c2","Trigger counters TDC histogram",800,650);
  c2->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c2->cd(4*i+j+1);
      gPad->SetLogy();
      hist_Tt[i][j]->SetLineColor(kBlack);
      hist_Tt[i][j]->GetListOfFunctions()->Remove(f_t[i][j]);
      hist_Tt[i][j]->Draw();
      hist_Tt_cut[i][j]->SetLineColor(kRed);
      hist_Tt_cut[i][j]->SetFillColor(kRed);
      hist_Tt_cut[i][j]->SetFillStyle(3001);
      hist_Tt_cut[i][j]->Draw("sames");
      
    }
  }

  Double_t one_photon= (15.4939+15.753+16.1096+16.0168)*0.945*0.91*0.5/4;
  Double_t one_photon_ind= (15.4939+15.753+16.1096+16.0168)*0.945/4;
  Double_t ind_gain[4];
  ind_gain[0] = 15.4939*0.959;
  ind_gain[1] = 15.753*0.995;
  ind_gain[2] = 16.1096*0.858;
  ind_gain[3] = 16.0168*0.968;

  Double_t numpho;
  Double_t rawadc;

  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      numpho = 0;
      rawadc = 0;
      pass = 0;

      //Trigger counter cut condition
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
	    pass+=1;
	    
	  }
	}
      }
      if(pass==4){
	hist_suma[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/one_photon);
	hist_sumt[i]->Fill(TDCs[i][0][0]);
	for(int j=0;j<4;j++){
	  numpho += (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
	  rawadc+=ADCi[i][j]-parameter_pe[j][1];
	  hist_indt[i][j]->Fill(TDCi[i][j][0]);
	}
	hist_inda[i]->Fill(numpho);
	ind_sum[i]->Fill(rawadc,ADCs[i][0]-parameter_pe[4][1]);
      }
    }
  }

  //Fitting BAC ADC
  TCanvas *c3 = new TCanvas("c3","BAC ADC of Ind and Sum");
  c3->Divide(N);
  Double_t pa_suma[N][3];
  Double_t pa_inda[N][3];
  for(int i=0;i<N;i++){
    std::cout<<i<<std::endl;
    c3->cd(i+1);
    hist_suma[i]->SetLineColor(kBlack);
    hist_suma[i]->Fit(f_suma[i],"Q","",-10,60);
    f_suma[i]->GetParameters(pa_suma[i]);
    
    hist_inda[i]->SetLineColor(kBlue);
    hist_inda[i]->Draw("sames");
    hist_inda[i]->Fit(f_inda[i],"Q","",-10,60);
    f_inda[i]->GetParameters(pa_inda[i]);
  }

  
  TCanvas *c4 = new TCanvas("c4","ADC of Ind channel and Sum channel",800,650);
  c4->Divide(N,1);
  for(int i=0;i<N;i++){
    ind_sum[i]->SetTitle("Comparision of Ind and Sum ADC;Ind [ADC Ch.];Sum [ADC Ch.]");
    c4->cd(i+1);
    ind_sum[i]->Draw("colz");
  }
 
  
  Double_t numpho_sum;
  Double_t pass_inda;
  Int_t tot_evt[N];
  Int_t evt_suma[N];
  Int_t evt_inda[N];
  
  //Determine ADC cut
  for(int i=0;i<N;i++){
    tot_evt[i] = 0;
    evt_suma[i] = 0;
    evt_inda[i] = 0;
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      numpho = 0;
      pass = 0;
      pass_inda = 0;

      //Trigger counter cut condition
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
	    pass+=1;
	    
	  }
	}
      }
      if(pass==4){
	tot_evt[i]+=1;
	if(ADCs[i][0]<3840){
	  numpho_sum = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	  if(numpho_sum>0.5){
	    evt_suma[i]+=1;
	    hist_suma_cut[i]->Fill(numpho_sum);
	    hist_sumt_cut[i]->Fill(TDCs[i][0][0]);
	  }
	}
	
	for(int j=0;j<4;j++){
	  if(ADCi[i][j]<3840)pass_inda+=1;
	  numpho += (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
	}
	if(pass_inda==4){
	  if(numpho>0.5){
	    evt_inda[i]+=1;
	    hist_inda_cut[i]->Fill(numpho);
	    for(int j=0;j<4;j++){
	      hist_indt_cut[i][j]->Fill(TDCi[i][j][0]);
	    }
	  }
	}
      }
    }
  }


  TCanvas *c5 = new TCanvas("c4","BAC Sum ADC histogram",800,650);
  c5->Divide(N);
  for(int i=0;i<N;i++){
    c5->cd(i+1);
    hist_suma[i]->SetLineColor(kBlack);
    hist_suma[i]->GetListOfFunctions()->Remove(f_suma[i]);
    hist_suma[i]->SetTitle(Form("Sum ADC %d cm;ADC [Ch.];n",start_pos+i*2));
    hist_suma[i]->Draw();
    
    hist_suma_cut[i]->SetLineColor(kRed);
    hist_suma_cut[i]->SetFillColor(kRed);
    hist_suma_cut[i]->SetFillStyle(3001);
    hist_suma_cut[i]->Draw("sames");
  }

  TCanvas *c6 = new TCanvas("c5","BAC Sum TDC histogram",800,650);
  c6->Divide(N);
  for(int i=0;i<N;i++){
    c6->cd(i+1);
    hist_sumt[i]->SetLineColor(kBlack);
    hist_sumt[i]->SetTitle(Form("Sum TDC %d cm;TDC [Ch.];n",start_pos+i*2));
    hist_sumt[i]->Draw();
    
    hist_sumt_cut[i]->SetLineColor(kRed);
    hist_sumt_cut[i]->SetFillColor(kRed);
    hist_sumt_cut[i]->SetFillStyle(3001);
    hist_sumt_cut[i]->Draw("sames");
  }

  TCanvas *c7 = new TCanvas("c6","BAC Ind ADC histogram",800,650);
  c7->Divide(N);
  for(int i=0;i<N;i++){
    c7->cd(i+1);
    hist_inda[i]->SetLineColor(kBlack);
    hist_inda[i]->GetListOfFunctions()->Remove(f_inda[i]);
    hist_inda[i]->SetTitle(Form("Ind ADC %d cm;ADC [Ch.];n",start_pos+i*2));
    hist_inda[i]->Draw();
    
    hist_inda_cut[i]->SetLineColor(kRed);
    hist_inda_cut[i]->SetFillColor(kRed);
    hist_inda_cut[i]->SetFillStyle(3001);
    hist_inda_cut[i]->Draw("sames");
  }

  TCanvas *c8 = new TCanvas("c7","BAC Ind TDC histogram",800,650);
  c8->Divide(N,4);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c8->cd(4*i+j+1);
      hist_indt[i][j]->SetLineColor(kBlack);
      hist_indt[i][j]->SetTitle(Form("Ind TDC %d cm;TDC [Ch.];n",start_pos+i*2));
      hist_indt[i][j]->Draw();
      
      hist_indt_cut[i][j]->SetLineColor(kRed);
      hist_indt_cut[i][j]->SetFillColor(kRed);
      hist_indt_cut[i][j]->SetFillStyle(3001);
      hist_indt_cut[i][j]->Draw("sames");
    }
  }

  
		      
      
    
      
	
    
    
      

  
} 
