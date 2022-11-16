Double_t timewalkfit(Double_t *x, Double_t *par){
  return -par[0]+par[1]/std::sqrt(x[0]);
}

void analysis(){

  //gStyle->SetOptStat(0);


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
  TH1D* hist_pe_indsum = new TH1D("hist_pe_indsum","hist_pe_indsum",250,700,1200);
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
    hist_pe_indsum->Fill(pe_sum);
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


  for(int n=0;n<total_pe;n++){
    data_pe->GetEntry(n);
    pe_sum = 0;
    for(int i=0;i<4;i++){
      //pe_sum+=ADC_pe[i]-parameter_pe[i][1];
      pe_sum+=ADC_pe[i];
    }
    hist_pe_indsum->Fill(pe_sum);
  }


  TCanvas *c_back = new TCanvas("c_back","Adding background of each channel",800,650);
  c_back->cd();
  hist_pe_indsum->Draw();
  






  //Make Histogram, Fitting function
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
  TH1D *hist_inda_se[N][4];
  TH1D *hist_indt[N][4];
  TH1D *hist_Suma[N];
  TH1D *hist_Sumt[N];

  TH1D *hist_inda_raw[N][4];
  TH1D *hist_Suma_raw[N];
  TH1D *hist_inda_com[N];


  
  TH1D *hist_Suma_cut[N];
  TH1D *hist_Sumt_cut[N];
  TH1D *multi[N];
  TF1 *f_a[N][4];
  TF1 *f_Tt[N][4];
  TF1 *f_inda[N];
  TF1 *f_inda_se[N][4];
  TF1 *f_indt[N][4];
  TF1 *f_sumt[N];
  TF1 *f_suma[N];

  TF1 *f_inda_com[N];

  TF1 *f_inda_raw[N][4];
  TF1 *f_suma_raw[N];

  TF1 *f_time_sum[N];

  TH2D *ind_sum[N];
  TH2D *adc_tdc_sum[N];
  TH2D *adc_tdc_sum_c[N];
  TH2D *adc_tdc_ind[N][4];

  Int_t start_pos = -6;

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
    //if(tot<min)min = tot;
    for(int j=0;j<4;j++){
      hist_Ta[i][j] = new TH1D(Form("hist_T%da_%dcm",j+4,start_pos+2*i),Form("hist_T%da_%dcm",j+4,start_pos+2*i),450,0,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,start_pos+2*i),Form("hist_T%dt_%dcm",j+4,start_pos+2*i),250,600,850);

      hist_Ta_cut[i][j] = new TH1D(Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%da_%dcm_cut",j+4,start_pos+2*i),450,0,2000);
      hist_Tt_cut[i][j] = new TH1D(Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),Form("hist_T%dt_%dcm_cut",j+4,start_pos+2*i),250,600,850);
      
      
      hist_indt[i][j] = new TH1D(Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_indt_%dcm_BAC%d",start_pos+2*i,j+1),200,650,850);
      
      hist_inda_se[i][j] = new TH1D(Form("hist_inda_se_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_inda_se_%dcm_BAC%d",start_pos+2*i,j+1),200,-10,60);
      
      hist_inda_raw[i][j] = new TH1D(Form("hist_inda_raw_%dcm_BAC%d",start_pos+2*i,j+1),Form("hist_inda_raw_%dcm_BAC%d",start_pos+2*i,j+1),105,-50,1000);
      
      f_a[i][j] = new TF1(Form("f_Tt%da_%dcm",j+4,start_pos+2*i),"landau(0)",200,2000);
      f_Tt[i][j] = new TF1(Form("f_Tt%dt_%dcm",j+4,start_pos+2*i),"gaus(0)",650,730);

      f_indt[i][j] = new TF1(Form("f_indt_%dcm_BAC%d",start_pos+2*i,j+1),"gaus(0)",650,850);
      
      f_inda_se[i][j] = new TF1(Form("f_inda_se_%dcm_BAC%d",start_pos+2*i,j+1),"landau(0)",-10,60);

      f_inda_raw[i][j] = new TF1(Form("f_inda_raw_%dcm_BAC%d",start_pos+2*i,j+1),"landau(0)",-50,1000);

      adc_tdc_ind[i][j] = new TH2D(Form("adc_tdc_ind_%dcm_BAC%d",start_pos+2*i,j+1),Form("adc_tdc_ind_%dcm_BAC%d",start_pos+2*i,j+1),100,0,2000,100,670,770);
      
			  
			  
    }
    hist_Suma[i] = new TH1D(Form("hist_Suma_%dcm",start_pos+2*i),Form("hist_Suma_%dcm",start_pos+2*i),100,-10,60);
    hist_Sumt[i] = new TH1D(Form("hist_Sumt_%dcm",start_pos+2*i),Form("hist_Sumt_%dcm",start_pos+2*i),200,650,850);

    hist_Suma_raw[i] = new TH1D(Form("hist_Suma_raw_%dcm",start_pos+2*i),Form("hist_Suma_raw_%dcm",start_pos+2*i),105,-50,1000);

    hist_Suma_cut[i] = new TH1D(Form("hist_Suma_%dcm_cut",start_pos+2*i),Form("hist_Suma_%dcm_cut",start_pos+2*i),100,-10,60);
    hist_Sumt_cut[i] = new TH1D(Form("hist_Sumt_%dcm_cut",start_pos+2*i),Form("hist_Sumt_%dcm_cut",start_pos+2*i),200,650,850);

    hist_inda[i] = new TH1D(Form("hist_inda_%dcm",start_pos+2*i),Form("hist_inda_%dcm",start_pos+2*i),100,-10,60);

    hist_inda_com[i] = new TH1D(Form("hist_inda_com_%dcm",start_pos+2*i),Form("hist_inda_com_%dcm",start_pos+2*i),100,-10,60);
    
    f_sumt[i] = new TF1(Form("f_sumt_%dcm",start_pos+2*i),"gaus(0)",-100,100);
    f_suma[i] = new TF1(Form("f_suma_%dcm",start_pos+2*i),"landau(0)",-10,60);

    f_suma_raw[i] = new TF1(Form("f_suma_raw_%dcm",start_pos+2*i),"landau(0)",-50,1000);

    f_inda[i] = new TF1(Form("f_inda_%dcm",start_pos+2*i),"landau(0)",-10,60);
    f_inda_com[i] = new TF1(Form("f_inda_com_%dcm",start_pos+2*i),"landau(0)",-10,60);

    f_time_sum[i] = new TF1(Form("f_time_sum_%dcm",start_pos+2*i),timewalkfit,200,2000,2);
    ind_sum[i] = new TH2D(Form("ind_sum%d",i),Form("ind_sum%d",i),310,-500,2500,200,-500,1500);

    multi[i] = new TH1D(Form("multi%d",i),"multiplicity",6,0,6);

    adc_tdc_sum[i] = new TH2D(Form("adc_tdc_sum_%dcm",start_pos+2*i),Form("adc_tdc_sum_%dcm",start_pos+2*i),100,0,2000,40,670,710);
    adc_tdc_sum_c[i] = new TH2D(Form("adc_tdc_sum_%dcm_c",start_pos+2*i),Form("adc_tdc_sum_%dcm_c",start_pos+2*i),100,0,2000,100,-50,50);
    
    

    
  }


  //Trigger counters Cut condition
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
      //std::cout<<"ADC MPV : "<<pa_a[i][j][1]<<" and sigma : "<<pa_a[i][j][2]<<std::endl;
      hist_Tt[i][j]->Fit(f_Tt[i][j],"Q0","",650,730);
      f_Tt[i][j]->GetParameters(pa_t[i][j]);
      //std::cout<<"TDC mean : "<<pa_t[i][j][1]<<" and sigma : "<<pa_t[i][j][2]<<std::endl;
    }
  }

  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      pass = 0;
      for(int j=0;j<4;j++){
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+50*pa_a[i][j][2]){
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
  //c1->Divide(2,2);
  
  for(int i=0;i<N;i++){
    //for(int i=3;i<4;i++){
    for(int j=0;j<4;j++){
        entry_Ta[i][j] = hist_Ta[i][j]-> GetEntries();
	entry_Tt[i][j] = hist_Tt[i][j]-> GetEntries();
	entry_Ta_cut[i][j] = hist_Ta_cut[i][j]-> GetEntries();
	entry_Tt_cut[i][j] = hist_Ta_cut[i][j]-> GetEntries();
	
      
      c1->cd(4*i+j+1);
      //c1->cd(j+1);
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
  //c2->Divide(2,2);
  for(int i=0;i<N;i++){
    //for(int i=3;i<4;i++){
    for(int j=0;j<4;j++){
      c2->cd(4*i+j+1);
      //c2->cd(j+1);
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
  
  


  TFile *file3[N];
  TTree *data3[N];
  Double_t ADC_sum[N];
  //gDirectory -> cd("Rint:/");  //Change to Working directory (default)
  
  for(int i=0;i<N;i++){
    file3[i] = new TFile(Form("../../ELPH_data/exp_data/sum_%dcm.root",start_pos+2*i),"recreate");
    data3[i] = new TTree("tree","tree");
    data3[i]->Branch("ADC",&ADC_sum[i],"ADC/D");  //Branch(branchname, &value);
    
  }

  
  

  
  Double_t one_photon= (15.4939+15.753+16.1096+16.0168)*0.945*0.91/4;
  Double_t one_photon_ind= (15.4939+15.753+16.1096+16.0168)*0.945/4;
  Double_t ind_gain[4];
  ind_gain[0] = 15.4939*0.959;
  ind_gain[1] = 15.753*0.995;
  ind_gain[2] = 16.1096*0.858;
  //ind_gain[3] = 16.0168*0.968*2;
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
	if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+50*pa_a[i][j][2]){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-3*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+3*pa_t[i][j][2]){
	    pass+=1;
	  }
	}
      }
      if(pass==4){
	adc_tdc_sum[i]->Fill(ADCs[i][0],TDCs[i][0][0]);
      }
    }
  }

  Double_t time_sum[N][2];
  TCanvas *c_2d_sum = new TCanvas("c_2d_sum","2D histogram of SUM",800,650);
  c_2d_sum->Divide(N);
  for(int i=0;i<N;i++){
    c_2d_sum->cd(i+1);
    if(i==0||i==N-1)adc_tdc_sum[i]->Fit(f_time_sum[i],"","",200,1000);
    else{adc_tdc_sum[i]->Fit(f_time_sum[i],"","",200,1500);}
    adc_tdc_sum[i]->Draw("colz");
    f_time_sum[i]->GetParameters(time_sum[i]);
    
  }

	

  
  for(int i=0;i<N;i++){
    file3[i]->cd();
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
	adc_tdc_sum_c[i]->Fill(ADCs[i][0],TDCs[i][0][0]+time_sum[i][0]-(time_sum[i][1]/std::sqrt(ADCs[i][0])));

	for(int j=0;j<4;j++)hist_indt[i][j]->Fill(TDCi[i][j][0]);
	//hist_Sumt[i]->Fill(TDCs[i][0][0]+time_sum[i][0]-(time_sum[i][1]/std::sqrt(ADCs[i][0])));
	hist_Sumt[i]->Fill(TDCs[i][0][0]);
      }
    }
  }


  TCanvas *c3 = new TCanvas("c3","BAC SUM TDC histogram",800,650);
  c3->Divide(N);
  for(int i=0;i<N;i++){
    c3->cd(i+1);
    if(i==0||i==N-1)hist_Sumt[i]->Fit(f_sumt[i],"Q","",650,750);
    else{hist_Sumt[i]->Fit(f_sumt[i],"Q","",650,850);}
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
  Int_t pass_ind_sum[N];
  Double_t adc_ind_to;
  Double_t ind_rawadc;
  

  Double_t tdc_c[N];
  for(int i=0;i<N;i++){
    file3[i]->cd();
    tot_evt[i] = 0;
    pass_evt[i] = 0;
    pass_mul[i] = 0;
    pass_ind_sum[i] = 0;
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);

      //Initialize
      pass = 0;
      numpho = 0;
      pass_ind = 0;
      mul_count = 0;
      adc_ind_to = 0;
      ind_rawadc = 0;
      
	for(int j=0;j<4;j++){
	  numpho_ind[j] = 0;
	  if(Ta[i][j][0]>pa_a[i][j][1]-3*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+50*pa_a[i][j][2]){
	    if(Tt[i][j][0][0]>pa_t[i][j][1]-3*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+3*pa_t[i][j][2]){
	      pass+=1;
	      
	    }
	  }
	}
	if(pass==4){
	  tot_evt[i]+=1;
	  hist_Suma[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/one_photon);
	  hist_Suma_raw[i]->Fill(ADCs[i][0]-parameter_pe[4][1]);
	  for(int j=0;j<4;j++){
	    hist_inda_raw[i][j]->Fill(ADCi[i][j]-parameter_pe[j][1]);
	    hist_inda_se[i][j] ->Fill((ADCi[i][j]-parameter_pe[j][1])/ind_gain[j]);
	  }
	    
	  //tdc_c[i] = TDCs[i][0][0]+time_sum[i][0]-(time_sum[i][1]/std::sqrt(ADCs[i][0]));
	  //if(tdc_c[i]>pa_sumt[i][1]-7*pa_sumt[i][2]&&tdc_c[i]<pa_sumt[i][1]+5*pa_sumt[i][2]){
	  //if(TDCs[i][0][0]>pa_sumt[i][1]-7*pa_sumt[i][2]&&TDCs[i][0][0]<pa_sumt[i][1]+5*pa_sumt[i][2]){
	    hist_Suma_cut[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/one_photon);
	    //hist_Sumt_cut[i]->Fill(TDCs[i][0][0]+time_sum[i][0]-(time_sum[i][1]/std::sqrt(ADCs[i][0])));
	    hist_Sumt_cut[i]->Fill(TDCs[i][0][0]);
	    
	    if((ADCs[i][0]-parameter_pe[4][1])/one_photon>0.99)pass_evt[i]+=1;
	    ADC_sum[i] = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	    
	    data3[i]->Fill();
	    //}

	  for(int j=0;j<4;j++){
	    ind_rawadc += ADCi[i][j]-parameter_pe[j][1];
	    
	    //if(TDCi[i][j][0]>pa_indt[i][j][1]-3*pa_indt[i][j][2]&&TDCi[i][j][0]<pa_indt[i][j][1]+3*pa_indt[i][j][2]){
	      numpho += (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
	      adc_ind_to +=ADCi[i][j]-parameter_pe[j][1];
	      numpho_ind[j] = (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
		//hist_inda_se[i][j] ->Fill(numpho_ind[j]);
		
	      //}
	  }

	  hist_inda_com[i]->Fill(ind_rawadc/one_photon_ind);

	  hist_inda[i]->Fill(numpho);
	  if(numpho>0.99)pass_ind_sum[i]+=1;
	  //ind_sum[i]->Fill(numpho,(ADCs[i][0]-parameter_pe[4][1])/one_photon);
	  ind_sum[i]->Fill(adc_ind_to,ADCs[i][0]-parameter_pe[4][1]);
	  for(int j=0;j<4;j++)adc_tdc_ind[i][j]->Fill(ADCi[i][j],TDCi[i][j][0]);

	    //multiplicity
	  for(int j=0;j<4;j++){
	    if(numpho_ind[j]>0.5)mul_count+=1;
	  }
	  multi[i]->Fill(mul_count);
	  if(mul_count>0.99)pass_mul[i]+=1;
	    
	}
	  
	
	  
	//}
	//pass=0;
	  
	
    }
    data3[i]->Write();
    file3[i]->Close();
  }

for(int i=0;i<N;i++){
  std::cout<<"Efficiency(SUM) for position "<<start_pos+i*2<<" cm is "<<1.00*pass_evt[i]/tot_evt[i]<<std::endl;
  std::cout<<"Efficiency(SUM of Ind) for position "<<start_pos+i*2<<" cm is "<<1.00*pass_ind_sum[i]/tot_evt[i]<<std::endl;
  std::cout<<"Efficiency(Multiplicity) for position "<<start_pos+i*2<<" cm is "<<1.00*pass_mul[i]/tot_evt[i]<<std::endl;
 }
  
  


 TCanvas *c_2d_sum_c = new TCanvas("c_2d_sum_c","2D histogram of SUM after timewalk correction",800,650);
  c_2d_sum_c->Divide(N);
  for(int i=0;i<N;i++){
    c_2d_sum_c->cd(i+1);
    adc_tdc_sum_c[i]->Draw("colz");

  }
  
 TCanvas *ccut = new TCanvas("ccut","BAC TDC cut compare",800,650);
  ccut->Divide(N);
  for(int i=0;i<N;i++){
      ccut->cd(i+1);
      gPad->SetLogy();
      hist_Sumt[i]->SetLineColor(kBlack);
      hist_Sumt[i]->GetListOfFunctions()->Remove(f_sumt[i]);
      hist_Sumt[i]->Draw();
      hist_Sumt_cut[i]->SetLineColor(kRed);
      hist_Sumt_cut[i]->SetFillColor(kRed);
      hist_Sumt_cut[i]->SetFillStyle(3001);
      hist_Sumt_cut[i]->Draw("sames");
      
    
  }
  


  

 


  TCanvas *c_2d = new TCanvas("c_2d","2D histogram of individual",800,650);
  c_2d->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c_2d->cd(4*i+j+1);
      adc_tdc_ind[i][j]->Draw("colz");
    }
  }

  TCanvas *c_ind_se = new TCanvas("c_ind_se","ADC histogram of Individual Channel",800,650);
  c_ind_se->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c_ind_se->cd(4*i+j+1);
      hist_inda_se[i][j]->Draw("colz");
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
    
    //hist_inda[i]->SetLineColor(kBlack);
    //hist_inda[i]->Draw("same");
    hist_Suma[i]->SetLineColor(kBlack);
    hist_Suma_cut[i]->Fit(f_suma[i],"Q","",0,100);
    //hist_Suma[i]->GetListOfFunctions()->Remove(f_suma[i]);
    hist_Suma[i]->Draw();
    
    hist_Suma_cut[i]->SetLineColor(kRed);
    hist_Suma_cut[i]->SetFillColor(kRed);
    hist_Suma_cut[i]->SetFillStyle(3001);
    //if(i==0)hist_Suma_cut[i]->Fit(f_suma[i],"Q","",-10,10);
    //hist_Suma_cut[i]->Fit(f_suma[i],"Q","",0,100);
    hist_Suma_cut[i]->Draw("sames");
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
    hist_inda[i]->Fit(f_inda[i],"Q","",2,100);
    f_inda[i]->GetParameters(pa_inda[i]);
    
    npe_pos_ind[i] = pa_inda[i][1];
    npe_error_ind[i] = f_inda[i]->GetParError(1);
  }

  Double_t pa_inda_com[N][3];
  TCanvas * c50 = new TCanvas("c50","BAC indivisual (divide after summing) ADC histogram",800,650);
  c50->Divide(N);
  Double_t npe_pos_ind_com[N];
  Double_t npe_error_ind_com[N];
  for(int i=0;i<N;i++){
    c50->cd(i+1);
    hist_inda_com[i]->Fit(f_inda_com[i],"Q","",-10,60);
    f_inda_com[i]->GetParameters(pa_inda_com[i]);
    
    npe_pos_ind_com[i] = pa_inda_com[i][1];
    npe_error_ind_com[i] = f_inda_com[i]->GetParError(1);
  }

  //lea->Draw();

  //Simulation - position dependence
  TFile *file_simul[N];
  TTree *data_simul[N];
  Int_t simul_npe[N];
  TH1D *hist_simul[N];
  TF1 *fit_simul[N];
  
  file_simul[0] = new TFile("../../ELPH_data/simul_data/elph_221020_m6_0_simul.root","read");
  file_simul[1] = new TFile("../../ELPH_data/simul_data/elph_221020_m4_0_simul.root","read");
  file_simul[2] = new TFile("../../ELPH_data/simul_data/elph_221020_m2_0_simul.root","read");
  file_simul[3] = new TFile("../../ELPH_data/simul_data/elph_221020_0_0_simul.root","read");
  file_simul[4] = new TFile("../../ELPH_data/simul_data/elph_221020_2_0_simul.root","read");
  file_simul[5] = new TFile("../../ELPH_data/simul_data/elph_221020_4_0_simul.root","read");
  file_simul[6] = new TFile("../../ELPH_data/simul_data/elph_221020_6_0_simul.root","read");


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
  TGraphErrors *npe3 = new TGraphErrors(N,xpos,npe_pos_ind_com,x_error,npe_error_ind_com);
  
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
  npe3->SetMarkerStyle(24);
  npe3->SetMarkerColor(6);
  npe3->SetLineColor(6);
  npe3->SetMarkerSize(1);

  TMultiGraph *mg = new TMultiGraph();
  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  mg->Add(npe0);
  mg->Add(npe1);
  mg->Add(npe2);
  mg->Add(npe3);
  le_s->AddEntry(npe0,"ELPH test SUM");
  le_s->AddEntry(npe2,"ELPH test Indi.");
  le_s->AddEntry(npe3,"ELPH test Indi. 2");
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
  file_com[0] = new TFile("../../ELPH_data/exp_data/sum_0cm.root","read");
  file_com[1] = new TFile("../../ELPH_data/simul_data/elph_221020_m0.5_0_simul.root","read");
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

  Double_t pa_raw[N][5][3];
  TCanvas *craw = new TCanvas("craw","craw",800,650);
  craw->Divide(5,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      craw->cd(5*i+j+1);
      hist_inda_raw[i][j]->Fit(f_inda_raw[i][j],"","",0,1000);
      f_inda_raw[i][j]->GetParameters(pa_raw[i][j]);
    }
    int j = 4;
    craw->cd(5*i+j+1);
    hist_Suma_raw[i]->Fit(f_suma_raw[i],"","",0,1000);
    f_suma_raw[i]->GetParameters(pa_raw[i][j]);
  }

  TLegend *le_ratio = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<4;i++){
    le_ratio->AddEntry(hist_inda_se[4][i],Form("BAC%d",i+1));
  }
  le_ratio->AddEntry(hist_Suma[4],"BAC SUM");

  TCanvas *c_ratio = new TCanvas("c_ratio","compare between ind and SUM",800,650);
  c_ratio->cd();

  hist_Suma[2]->SetLineColor(1);
  hist_Suma[2]->Draw();
  for(int i=0;i<4;i++){
    hist_inda_se[2][i]->SetLineColor(i+6);
    hist_inda_se[2][i]->Draw("sames");
  }
  le_ratio->Draw();
    

  //Inefficiency Study
  TFile *file_pe1 = new TFile("../../ELPH_data/exp_data/run00333.root","read");
  TTree *data_pe1 = (TTree*)file_pe1->Get("tree");
  Double_t ADC_pe1[4];
  Double_t TDC_pe1[4][16];
  Double_t ADCs_pe1[1];
  Double_t TDCs_pe1[1][16];
  data_pe1 ->SetBranchAddress("E72BACa",ADC_pe1);
  data_pe1 ->SetBranchAddress("E72BACt",TDC_pe1);
  data_pe1 ->SetBranchAddress("E72BACSUMa",ADCs_pe1);
  data_pe1 ->SetBranchAddress("E72BACSUMt",TDCs_pe1);
  Double_t total_pe1 = data_pe1->GetEntries();

  TH1D* hist_ADC_pe1[4];
  TH1D* hist_TDC_pe1[4];

  //Range
  Double_t ADC_start = 50;
  Double_t ADC_end = 2000;
  Double_t TDC_start = 0;
  Double_t TDC_end = 1500;
  TH1D* hist_ADCs_pe1 = new TH1D("hist_ADCs_pe1","hist_ADCs_pe1",ADC_end-ADC_start,ADC_start,ADC_end);
  TH1D* hist_TDCs_pe1 = new TH1D("hist_TDCs_pe1","hist_TDCs_pe1",(TDC_end-TDC_start)/10,TDC_start,TDC_end);
  TF1* fit_ADCs_pe1 = new TF1("fit_ADCs_pe1","gaus(0)",ADC_start,ADC_end);
  TF1* fit_TDCs_pe1 = new TF1("fit_TDCs_pe1","gaus(0)",TDC_start,TDC_end);
  TF1* fit_ADC_pe1[4];
  TF1* fit_TDC_pe1[4];
  for(int j=0;j<4;j++){
    hist_ADC_pe1[j] = new TH1D(Form("hist_ADC_pe1_BAC%d",j+1),Form("hist_ADC_pe1_BAC%d",j+1),ADC_end-ADC_start,ADC_start,ADC_end);
    hist_TDC_pe1[j] = new TH1D(Form("hist_TDC_pe1_BAC%d",j+1),Form("hist_TDC_pe1_BAC%d",j+1),(TDC_end-TDC_start)/10,TDC_start,TDC_end);
    fit_ADC_pe1[j] = new TF1(Form("fit_ADC_pe1_BAC%d",j+1),"gaus(0)",ADC_start,ADC_end);
    fit_TDC_pe1[j] = new TF1(Form("fit_TDC_pe1_BAC%d",j+1),"gaus(0)",TDC_start,TDC_end);
  }

  for(int n=0;n<total_pe1;n++){
    data_pe1->GetEntry(n);
    for(int i=0;i<4;i++){
      if(TDC_pe1[i][0]>pa_indt[3][i][1]-3*pa_indt[3][i][2]&&TDC_pe1[i][0]<pa_indt[3][i][1]+3*pa_indt[3][i][2]){
	//if(ADC_pe1[i]>pa_inda[3][i][1]+3*pa_inda[3][i][2]){
	  hist_ADC_pe1[i]->Fill(ADC_pe1[i]);
	  hist_TDC_pe1[i]->Fill(TDC_pe1[i][0]);
	  //}
      }
    }
    if(TDCs_pe1[0][0]>pa_sumt[3][1]-3*pa_sumt[3][2]&&TDCs_pe1[0][0]<pa_sumt[3][1]+3*pa_sumt[3][2]){
      //if(ADCs_pe1[i]>pa_suma[3][1]+3*pa_suma[3][2]){
	hist_ADCs_pe1->Fill(ADCs_pe1[0]);
	hist_TDCs_pe1->Fill(TDCs_pe1[0][0]);
	//}
    }
  }

  Double_t parameter_ADC_pe1[5][3];
  Double_t parameter_TDC_pe1[5][3];

  TCanvas *c_pe1 = new TCanvas("c_pe1","ADC TDC w/o beam",800,650);
  c_pe1->Divide(5,2);
  for(int i=0;i<5;i++){
    if(i<4){
      c_pe1->cd(i+1);
      hist_ADC_pe1[i]->Fit(fit_ADC_pe1[i],"","",ADC_start,ADC_end);
      fit_ADC_pe1[i]->GetParameters(parameter_ADC_pe1[i]);

      c_pe1->cd(5+i+1);
      hist_TDC_pe1[i]->Fit(fit_TDC_pe1[i],"","",TDC_start,TDC_end);
      fit_TDC_pe1[i]->GetParameters(parameter_TDC_pe1[i]);
    }
    if(i==4){
      c_pe1->cd(i+1);
      hist_ADCs_pe1->Fit(fit_ADCs_pe1,"","",ADC_start,ADC_end);
      fit_ADCs_pe1->GetParameters(parameter_ADC_pe1[i]);

      c_pe1->cd(5+i+1);
      hist_TDCs_pe1->Fit(fit_TDCs_pe1,"","",TDC_start,TDC_end);
      fit_TDCs_pe1->GetParameters(parameter_TDC_pe1[i]);
    }
  }



  

  
      
      
      
  
  
  

  
     
}
