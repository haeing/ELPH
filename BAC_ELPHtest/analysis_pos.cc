//Y position : -35 -23 -11 -10 0 12

Int_t y_pos = 9999;

void analysis_pos(){
  //gStyle -> SetOptFit(1);
  gStyle -> SetOptStat(0);
  
  TFile *file_pe;
  if(y_pos==13||y_pos == 12 ||y_pos == 0||y_pos ==-10||y_pos==9999)file_pe = new TFile("../../ELPH_data/exp_data/run00333.root","read");
  else if(y_pos ==-35 || y_pos ==-23 ||y_pos ==-11)file_pe = new TFile("../../ELPH_data/exp_data/run00079.root","read");
  
  Int_t N;  //x position

  if(y_pos == 0)N=7;
  else if(y_pos ==-10)N=3;
  else if(y_pos ==-11)N=7;
  else if(y_pos ==-23)N=7;
  else if(y_pos ==-35)N=11;
  else if(y_pos == 12)N=5;
  else if(y_pos ==13) N=1;
  else if(y_pos ==9999) N=5;

  
  Int_t att = 1;
  

  //Pedestal
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

  Int_t x_pos[N];
  if(y_pos == 0 || y_pos == -11 || y_pos ==-23){
    for(int i=0;i<N;i++)x_pos[i] = -6+2*i;
  }
  else if(y_pos ==-10){
    for(int i=0;i<N;i++)x_pos[i] = -3+3*i;
  }
  else if(y_pos == -35){
    x_pos[0] = -6;
    x_pos[1] = -5;
    x_pos[2] = -4;
    x_pos[3] = -3;
    x_pos[4] = -2;
    x_pos[5] = 0;
    x_pos[6] = 2;
    x_pos[7] = 3;
    x_pos[8] = 4;
    x_pos[9] = 5;
    x_pos[10] = 6;
  }

  else if(y_pos==12){
    x_pos[0] = -4;
    x_pos[1] = -2;
    x_pos[2] = 0;
    x_pos[3] = 2;
    x_pos[4] = 4;
  }

  else if(y_pos==9999){
    x_pos[0] = 60;
    x_pos[1] = 70;
    x_pos[2] = 75;
    x_pos[3] = 80;
    x_pos[4] = 100;
  }

  else if(y_pos==13)x_pos[0]=0;
    

  Int_t start_pos = -6;
  TFile *file_po[N];
  TTree *data_po[N];

 
  Double_t ADCi[N][4];
  Double_t TDCi[N][4][16];
  Double_t ADCs[N][1];
  Double_t TDCs[N][1][16];
  Double_t Ta[N][4][1];
  Double_t Tt[N][4][1][16];

  TFile *file_simul[N][att];
  TTree *data_simul[N][att];

  Int_t nhmppc[N][4];

  Double_t tot[N];

  TH1D *hist_Ta[N][4];
  TH1D *hist_Tt[N][4];

  TH1D *hist_Ta_cut[N][4];
  TH1D *hist_Tt_cut[N][4];


  TH1D *hist_inda[N];
  TH1D *hist_suma[N];

  TH1D *hist_simul[N][att];

  TH1D *hist_inda_cut[N];
  TH1D *hist_suma_cut[N];

  TH1D *hist_indt[N][4];
  TH1D *hist_sumt[N];

  TH1D *hist_indt_cut[N][4];
  TH1D *hist_sumt_cut[N];

  TH2D *ind_sum[N];

  TH2D *sumt_2d[N];

  TH1D *hist_Ta_out[N][4];

  TF1 *f_a[N][4];
  TF1 *f_t[N][4];
  
  TF1 *f_inda[N];
  TF1 *f_suma[N];

  TF1 *f_simul[N][att];

  TF1 *f_indt[N][4];
  TF1 *f_sumt[N];

  Int_t att_start = 400;

  if(y_pos==0){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00313.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00314.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00315.root","read");
    //file_po[3] = new TFile("../../ELPH_data/exp_data/run00297.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00333.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00316.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00317.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00319.root","read");

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_0mm.root",x_pos[i]*10),"read");
      }
    }
  }


  else if(y_pos==-10){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00300.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00299.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00301.root","read");

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_-10mm.root",x_pos[i]*10),"read");
      }
    }
  }
  

  else if(y_pos==-11){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00065.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00064.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00062.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00063.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00061.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00060.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00059.root","read");

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_-11mm.root",x_pos[i]*10),"read");
      }
    }
  }
  

  else if(y_pos==-23){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00048.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00052.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00053.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00045.root","read");
    //file_po[3] = new TFile("../../ELPH_data/exp_data/run00079.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00054.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00055.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00056.root","read");

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_-23mm.root",x_pos[i]*10),"read");
      }
    }
  }


  else if(y_pos==-35){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00068.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00078.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00069.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00077.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00071.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00067.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00072.root","read");
    file_po[7] = new TFile("../../ELPH_data/exp_data/run00076.root","read");
    file_po[8] = new TFile("../../ELPH_data/exp_data/run00073.root","read");
    file_po[9] = new TFile("../../ELPH_data/exp_data/run00075.root","read");
    file_po[10] = new TFile("../../ELPH_data/exp_data/run00074.root","read");

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_-35mm.root",x_pos[i]*10),"read");
      }
    }
  }

  else if(y_pos==12){
    //file_po[0] = new TFile("../../ELPH_data/exp_data/run00328.root","read"); //-6
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00329.root","read"); //-4
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00330.root","read"); //-2
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00326.root","read"); //0
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00331.root","read"); //2
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00332.root","read"); //4

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_cal_att/elph_221020_%dmm_12mm.root",x_pos[i]*10),"read");
      }
    }
  }

  else if(y_pos==13){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00326.root","read"); //0

    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile(Form("../../ELPH_data/simul_acr/elph_221020_%dmm_12mm.root",x_pos[i]*10),"read");
      }
    }
  }

  else if(y_pos==9999){ //Threshold change
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00294.root","read"); //0
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00295.root","read"); //0
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00296.root","read"); //0
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00293.root","read"); //0
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00292.root","read"); //0

    
    for(int i=0;i<N;i++){
      for(int j=0;j<att;j++){
	file_simul[i][j] = new TFile("../../ELPH_data/simul_acr/elph_221020_-20mm_-10mm.root","read");
      }
    }
  }

  
    



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
      hist_Ta[i][j] = new TH1D(Form("hist_Ta%d_%dcm",j+4,x_pos[i]),Form("hist_Ta%d_%dcm",j+4,x_pos[i]),450,0,2000);
      hist_Ta_out[i][j] = new TH1D(Form("hist_Ta_out_%d_%dcm",j+4,x_pos[i]),Form("hist_Ta_out_%d_%dcm",j+4,x_pos[i]),450,0,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,x_pos[i]),Form("hist_T%dt_%dcm",j+4,x_pos[i]),250,600,850);

      hist_Ta_cut[i][j] = new TH1D(Form("hist_T%da_%dcm_cut",j+4,x_pos[i]),Form("hist_T%da_%dcm_cut",j+4,x_pos[i]),450,0,2000);
      hist_Tt_cut[i][j] = new TH1D(Form("hist_T%dt_%dcm_cut",j+4,x_pos[i]),Form("hist_T%dt_%dcm_cut",j+4,x_pos[i]),250,600,850);

      hist_indt[i][j] = new TH1D(Form("hist_indt_%dcm_BAC%d",x_pos[i],j+1),Form("hist_indt_%dcm_BAC%d",x_pos[i],j+1),500,0,1500);

      hist_indt_cut[i][j] = new TH1D(Form("hist_indt_cut_%dcm_BAC%d",x_pos[i],j+1),Form("hist_indt_cut_%dcm_BAC%d",x_pos[i],j+1),500,0,1500);
      
      f_a[i][j] = new TF1(Form("f_a%da_%dcm",j+4,x_pos[i]),"landau(0)",200,2000);
      f_t[i][j] = new TF1(Form("f_t%dt_%dcm",j+4,x_pos[i]),"gaus(0)",650,730);

      f_indt[i][j] = new TF1(Form("f_indt_%dcm_BAC%d",x_pos[i],j+1),"gaus(0)",0,1500);
    }
    
    hist_inda[i] = new TH1D(Form("hist_inda_%dcm",x_pos[i]),Form("hist_inda_%dcm",x_pos[i]),100,-10,100);
    hist_suma[i] = new TH1D(Form("hist_Suma_%dcm",x_pos[i]),Form("hist_Suma_%dcm",x_pos[i]),100,-20,100);
    hist_sumt[i] = new TH1D(Form("hist_Sumt_%dcm",x_pos[i]),Form("hist_Sumt_%dcm",x_pos[i]),500,0,1500);
    for(int j=0;j<att;j++){
      hist_simul[i][j] = new TH1D(Form("hist_simul_%dcm_%d",x_pos[i],att_start+50*j),Form("hist_simul_%dcm",x_pos[i]),100,-10,80);
    }

    hist_inda_cut[i] = new TH1D(Form("hist_inda_cut_%dcm",x_pos[i]),Form("hist_inda_cut_%dcm",x_pos[i]),100,-10,100);
    hist_suma_cut[i] = new TH1D(Form("hist_suma_cut_%dcm",x_pos[i]),Form("hist_suma_cut_%dcm",x_pos[i]),100,-10,100);
    hist_sumt_cut[i] = new TH1D(Form("hist_sumt_cut_%dcm",x_pos[i]),Form("hist_sumt_cut_%dcm",x_pos[i]),500,0,1500);

    ind_sum[i] = new TH2D(Form("ind_sum%d",i),Form("ind_sum%d",i),260,-100,2500,160,-100,1500);

    sumt_2d[i] = new TH2D(Form("sumt_2d%d",i),Form("sumt_2d%d",i),500,0,1500,500,0,1500);
    

    if(TMath::Abs(x_pos[i])<=5)f_inda[i] = new TF1(Form("f_inda_%dcm",x_pos[i]),"gaus(0)",-10,100);
    else if(TMath::Abs(x_pos[i])>5)f_inda[i] = new TF1(Form("f_inda_%dcm",x_pos[i]),"[0]*TMath::Poisson(x,[1])",0,100);
    f_sumt[i] = new TF1(Form("f_sumt_%dcm",x_pos[i]),"gaus(0)",0,1500);
    if(TMath::Abs(x_pos[i])<=5)f_suma[i] = new TF1(Form("f_suma_%dcm",x_pos[i]),"gaus(0)",-10,100);
    else if(TMath::Abs(x_pos[i])>5)f_suma[i] = new TF1(Form("f_suma_%dcm",x_pos[i]),"[0]*TMath::Poisson(x,[1])*TMath::Gaus(x,[2],[3])",-10,100);

    for(int j=0;j<att;j++){
      if(TMath::Abs(x_pos[i])<=5)f_simul[i][j] = new TF1(Form("f_simul_%dcm_%d",x_pos[i],att_start+j*50),"gaus(0)",-10,100);
      else if(TMath::Abs(x_pos[i])>5)f_simul[i][j] = new TF1(Form("f_simul_%dcm_%d",x_pos[i],att_start+j*50),"[0]*TMath::Poisson(x,[1])",0,10);
    }
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

  Int_t factor_a = 5;

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
	if(Ta[i][j][0]>pa_a[i][j][1]-factor_a*pa_a[i][j][2]&&Ta[i][j][0]<3840){
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

  Double_t one_photon= (15.4939+15.753+16.1096+16.0168)*0.945*0.5/4;
  Double_t one_photon_ind= (15.4939+15.753+16.1096+16.0168)*0.945/4;
  Double_t ind_gain[4];
  ind_gain[0] = 15.4939*0.959;
  ind_gain[1] = 15.753*0.995;
  ind_gain[2] = 16.1096*0.858;
  ind_gain[3] = 16.0168*0.968;

  Double_t numpho;
  Double_t rawadc;

  Int_t adc_factor = 5;

  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      numpho = 0;
      rawadc = 0;
      pass = 0;

      //Trigger counter cut condition
      if(y_pos!=12){
	for(int j=0;j<4;j++){
	  if(Ta[i][j][0]>pa_a[i][j][1]-factor_a*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	    if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
	      pass+=1;
	    
	    }
	  }
	}
	//if(pass==4){

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
      //}
      else if(y_pos==12){
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
  //Double_t pa_suma[N][3];
  Double_t pa_suma[N][4];
  Double_t pa_inda[N][3];
  for(int i=0;i<N;i++){
    c3->cd(i+1);
    hist_suma[i]->SetLineColor(kBlack);
    if(y_pos!=12){
      if(x_pos[i]==-6||x_pos[i]==6){
	f_suma[i]->SetParLimits(0,5000,50000);
	f_suma[i]->SetParLimits(1,0,5);
	f_suma[i]->SetParLimits(2,0,5);
	f_suma[i]->SetParLimits(3,0,100);
	hist_suma[i]->Fit(f_suma[i],"","",0,25);
      }
      else if(TMath::Abs(x_pos[i])<=5)hist_suma[i]->Fit(f_suma[i],"","",-10,60);
      
    }
    else if(y_pos==12){
      if(x_pos[i]!=-4)hist_suma[i]->Fit(f_suma[i],"Q","",0,60);
      else if(x_pos[i]==-4)hist_suma[i]->Fit(f_suma[i],"Q","",10,60);
    }
    f_suma[i]->GetParameters(pa_suma[i]);
    
    hist_inda[i]->SetLineColor(kBlue);
    hist_inda[i]->Draw("sames");
    if(y_pos!=12)hist_inda[i]->Fit(f_inda[i],"Q","",-10,60);
    if(y_pos==12)hist_inda[i]->Fit(f_inda[i],"Q","",0,60);
    f_inda[i]->GetParameters(pa_inda[i]);
  }
  c3->Close();

  //Fitting BAC TDC
  
  TCanvas *c3_1 = new TCanvas("c3_1","BAC TDC of Ind and Sum");
  c3_1->Divide(N);
  Double_t pa_sumt[N][3];
  Double_t pa_indt[N][4][3];
  for(int i=0;i<N;i++){
    c3_1->cd(i+1);
    hist_sumt[i]->SetLineColor(kBlack);
    hist_sumt[i]->Fit(f_sumt[i],"Q","",600,750);
    f_sumt[i]->GetParameters(pa_sumt[i]);
    
    for(int j=0;j<4;j++){
      hist_indt[i][j]->Fit(f_indt[i][j],"Q","",0,1500);
      f_indt[i][j]->GetParameters(pa_indt[i][j]);
    }
  }

  c3_1->Close();


  
  TCanvas *c4 = new TCanvas("c4","ADC of Ind channel and Sum channel",800,650);
  c4->Divide(4,2);
  for(int i=0;i<N;i++){
    ind_sum[i]->SetTitle("Comparision of Ind and Sum ADC;Ind [ADC Ch.];Sum [ADC Ch.]");
    c4->cd(i+1);
    ind_sum[i]->Draw("colz");
  }

  //auto eff = new TEfficiency("eff","Efficiency;X [mm];Efficiency",100,-60,70);
  auto eff = new TEfficiency("eff","Efficiency;Threshold [mV];Efficiency",100,50,110);
  
  Double_t numpho_sum;
  Double_t pass_inda;
  Int_t tot_evt[N];
  Int_t evt_suma[N];
  Int_t evt_sum_cuta[N];
  Int_t evt_inda[N];
  Int_t pass_t;
  Bool_t bPassed = 0;
  
  //Determine ADC cut
  for(int i=0;i<N;i++){
    tot_evt[i] = 0;
    evt_suma[i] = 0;
    evt_inda[i] = 0;
    evt_sum_cuta[i] = 0;
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      numpho = 0;
      pass = 0;
      pass_inda = 0;
      pass_t = 0;

      //Trigger counter cut condition
      if(y_pos!=12){
	for(int j=0;j<4;j++){
	  if(Ta[i][j][0]>pa_a[i][j][1]-factor_a*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	    if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
	      pass+=1;
	    
	    }
	  }
	}
	//if(pass==4){
	  tot_evt[i]+=1;
	  if(ADCs[i][0]<3840){
	    numpho_sum = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	    //if(numpho_sum>2){
	    for(int k=0;k<5;k++){
	      if(TDCs[i][0][k]>pa_sumt[i][1]-5*pa_sumt[i][2]&&TDCs[i][0][k]<pa_sumt[i][1]+5*pa_sumt[i][2])pass_t+=1;
	      
	    }
	    if(pass_t==1){
	      evt_suma[i]+=1;
	      hist_suma_cut[i]->Fill(numpho_sum);
	      for(int j=0;j<4;j++)hist_Ta_out[i][j]->Fill(Ta[i][j][0]);
	      hist_sumt_cut[i]->Fill(TDCs[i][0][0]);
	      bPassed=1;
	    }
	    else if(pass_t!=1)bPassed=0;
	    //}
	    if(numpho_sum>5)evt_sum_cuta[i]+=1;
	  }
	  else{bPassed=0;}
	  //eff->Fill(bPassed,x_pos[i]*10.0);
	  eff->Fill(bPassed,x_pos[i]*1.0);
	
	  for(int j=0;j<4;j++){
	    if(ADCi[i][j]<3840)pass_inda+=1;
	    numpho += (ADCi[i][j]-parameter_pe[j][1])/ind_gain[j];
	  }
	  if(pass_inda==4){
	    //if(numpho>2){
	    evt_inda[i]+=1;
	    hist_inda_cut[i]->Fill(numpho);
	    for(int j=0;j<4;j++){
	      hist_indt_cut[i][j]->Fill(TDCi[i][j][0]);
	    }
	    //}
	  }
	  //}
      }
      
      
      else if(y_pos==12){
	tot_evt[i]+=1;
	if(ADCs[i][0]<3840){
	  numpho_sum = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	  for(int k=0;k<5;k++){
	    if(TDCs[i][0][k]>pa_sumt[i][1]-5*pa_sumt[i][2]&&TDCs[i][0][k]<pa_sumt[i][1]+5*pa_sumt[i][2])pass_t+=1;
	  }
	  if(pass_t==1){
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
	  evt_inda[i]+=1;
	  hist_inda_cut[i]->Fill(numpho);
	  for(int j=0;j<4;j++){
	    hist_indt_cut[i][j]->Fill(TDCi[i][j][0]);
	  }
	}
      }
    }
  
  }


  TCanvas *c_eff = new TCanvas("c_eff","Efficiency using TEff",800,650);
  c_eff->cd();
  eff->Draw("");
  

  TCanvas *c1_1 = new TCanvas("c1_1","Trigger counters ADC histogram",800,650);
  c1_1->Divide(4,N);
  
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c1_1->cd(4*i+j+1);
      gPad->SetLogy();
      hist_Ta[i][j]->SetLineColor(kBlack);
      hist_Ta[i][j]->Draw();

      hist_Ta_out[i][j]->SetLineColor(kRed);
      hist_Ta_out[i][j]->SetFillColor(kRed);
      hist_Ta_out[i][j]->SetFillStyle(3001);
      hist_Ta_out[i][j]->Draw("sames");
      
    }
  }


  TCanvas *c5 = new TCanvas("c5","BAC Sum ADC histogram",800,650);
  c5->Divide(4,3);
  for(int i=0;i<N;i++){
    c5->cd(i+1);
    hist_suma[i]->SetLineColor(kBlack);
    hist_suma[i]->GetListOfFunctions()->Remove(f_suma[i]);
    hist_suma[i]->SetTitle(Form("Sum ADC %d cm;N_{p.e.};n",x_pos[i]));
    hist_suma[i]->Draw();
    
    hist_suma_cut[i]->SetLineColor(kRed);
    hist_suma_cut[i]->SetFillColor(kRed);
    hist_suma_cut[i]->SetFillStyle(3001);
    hist_suma_cut[i]->Draw("sames");
  }

  TCanvas *c6 = new TCanvas("c6","BAC Sum TDC histogram",800,650);
  c6->Divide(4,3);
  for(int i=0;i<N;i++){
    c6->cd(i+1);
    gPad->SetLogy();
    hist_sumt[i]->SetLineColor(kBlack);
    hist_sumt[i]->SetTitle(Form("Sum TDC %d cm;TDC [Ch.];n",x_pos[i]));
    hist_sumt[i]->Draw();
    
    hist_sumt_cut[i]->SetLineColor(kRed);
    hist_sumt_cut[i]->SetFillColor(kRed);
    hist_sumt_cut[i]->SetFillStyle(3001);
    hist_sumt_cut[i]->Draw("sames");
  }

  TCanvas *c7 = new TCanvas("c7","BAC Ind ADC histogram",800,650);
  c7->Divide(N);
  for(int i=0;i<N;i++){
    c7->cd(i+1);
    hist_inda[i]->SetLineColor(kBlack);
    hist_inda[i]->GetListOfFunctions()->Remove(f_inda[i]);
    hist_inda[i]->SetTitle(Form("Ind ADC %d cm;ADC [Ch.];n",x_pos[i]));
    hist_inda[i]->Draw();
    
    hist_inda_cut[i]->SetLineColor(kRed);
    hist_inda_cut[i]->SetFillColor(kRed);
    hist_inda_cut[i]->SetFillStyle(3001);
    hist_inda_cut[i]->Draw("sames");
  }

  TCanvas *c8 = new TCanvas("c8","BAC Ind TDC histogram",800,650);
  c8->Divide(N,4);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c8->cd(4*i+j+1);
      gPad->SetLogy();
      hist_indt[i][j]->SetLineColor(kBlack);
      hist_indt[i][j]->SetTitle(Form("Ind TDC %d cm;TDC [Ch.];n",x_pos[i]));
      hist_indt[i][j]->Draw();
      
      hist_indt_cut[i][j]->SetLineColor(kRed);
      hist_indt_cut[i][j]->SetFillColor(kRed);
      hist_indt_cut[i][j]->SetFillStyle(3001);
      hist_indt_cut[i][j]->Draw("sames");
    }
  }

  
  
  


  //Efficiency
  Double_t x_pos_e[N];
  Double_t x_error[N];
  Double_t effi_suma[N];
  Double_t effi_suma_err[N];
  Double_t effi_inda[N];
  Double_t npe_pos[N];
  Double_t npe_error[N];
  for(int i=0;i<N;i++){
    x_pos_e[i] = x_pos[i]*1.0;
    effi_suma[i] = 1.0*evt_suma[i]/tot_evt[i];
    effi_suma_err[i] = effi_suma[i]*TMath::Sqrt(1.0/evt_suma[i]+1.0/tot_evt[i]);
    effi_inda[i] = 1.0*evt_inda[i]/tot_evt[i];
    if(pa_suma[i][1]>=0)npe_pos[i] = pa_suma[i][1];
    else if(pa_suma[i][1]<0)npe_pos[i] = 0;
    npe_error[i] = f_suma[i]->GetParError(1);
    x_error[i] = 0.05;
  }

  for(int i=0;i<N;i++)cout<<"Efficiency by cutting ADC histogram : "<<1.0*evt_sum_cuta[i]/tot_evt[i]<<endl;

  

  //Npe simulation & experiment
  for(int i=0;i<N;i++){
    for(int j=0;j<att;j++){
      data_simul[i][j] = (TTree*)file_simul[i][j]->Get("tree");
      data_simul[i][j] ->SetBranchAddress("nhMppc",&nhmppc[i][j]);
    }
  }

  for(int i=0;i<N;i++){
    for(int j=0;j<att;j++){
      for(int n=0;n<5000;n++){
	data_simul[i][j]->GetEntry(n);
	hist_simul[i][j]->Fill(nhmppc[i][j]);
      }
    }
  }

  Double_t pa_simul[N][att][3];
  Double_t npe_pos_simul[att][N];
  Double_t npe_error_simul[att][N];
  TCanvas *c9 = new TCanvas("c9","BAC npe simulation",800,650);
  c9->Divide(N,att);
  for(int i=0;i<N;i++){
    for(int j=0;j<att;j++){
      c9->cd(att*i+j+1);
      hist_simul[i][j]->SetTitle(Form("Sum ADC %d cm %d;ADC [Ch.];n",x_pos[i],att_start+50*j));
      if(TMath::Abs(x_pos[i])>5){
	f_simul[i][j]->SetParLimits(0,5000,20000);
	f_simul[i][j]->SetParLimits(1,0,3);
	hist_simul[i][j]->Fit(f_simul[i][j],"","",0,10);
      }
      if(TMath::Abs(x_pos[i])<=5)hist_simul[i][j]->Fit(f_simul[i][j],"","",-10,80);
      

      f_simul[i][j]->GetParameters(pa_simul[i][j]);
      npe_pos_simul[j][i] = pa_simul[i][j][1];
      //if(pa_simul[i][j][1]>=0)npe_pos_simul[j][i] = TMath::Abs(pa_simul[i][j][1]-pa_suma[i][1])/pa_suma[i][1];
      npe_error_simul[j][i] = f_simul[i][j]->GetParError(1);
    }
  }
  
  
  TGraphErrors *npe0 = new TGraphErrors(N,x_pos_e,npe_pos,x_error,npe_error);
  TGraphErrors *npe[att];
  for(int i=0;i<att;i++){
    npe[i] = new TGraphErrors(N,x_pos_e,npe_pos_simul[i],x_error,npe_error_simul[i]);
    npe[i]->SetMarkerStyle(22);
    npe[i]->SetMarkerColor(2+i);
    npe[i]->SetLineColor(2+i);
    npe[i]->SetMarkerSize(2);
  }
  
  npe0->SetMarkerStyle(24);
  npe0->SetMarkerColor(1);
  npe0->SetLineColor(1);
  npe0->SetMarkerSize(2);

  auto mg_npe = new TMultiGraph();
  mg_npe->Add(npe0);
  for(int i=0;i<att;i++){
    mg_npe->Add(npe[i]);
  }

  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  le_s->AddEntry(npe0,"Experiment");
  for(int i=0;i<att;i++){
    le_s->AddEntry(npe[i],"Simulation");
  }
  
  TCanvas *c10 = new TCanvas("c10","Npe");
  c10->cd();
  mg_npe->SetTitle("Npe per each position;X [cm];N_{p.e.}");
  mg_npe->Draw("AP");
  le_s->Draw();
  

  TGraphErrors *eff_sum = new TGraphErrors(N,x_pos_e,effi_suma,x_error,effi_suma_err);
  eff_sum ->SetMarkerStyle(22);
  eff_sum ->SetMarkerSize(1.5);
  eff_sum ->SetMarkerColor(kBlack);
  TGraph *eff_ind = new TGraph(N,x_pos_e,effi_inda);
  eff_ind ->SetMarkerStyle(24);
  eff_ind ->SetMarkerColor(kBlue);
  auto mg = new TMultiGraph();
  mg->Add(eff_sum);
  mg->Add(eff_ind);
  TCanvas *c11 = new TCanvas("c11","Efficiency");
  c11->cd();
  eff_sum->SetTitle("Efficiency;X [cm];Efficiency");
  eff_sum->Draw("AP");


  

  Double_t mean[N];
  Double_t sigma[N];
  Double_t error[N];
  
  Double_t simul_mean[N];
  Double_t simul_sigma[N];
  Double_t simul_error[N];
  Int_t x_position[N];




		      
  TFile *file_save = new TFile(Form("parameter_%d.root",y_pos),"recreate");
  TTree* tree = new TTree("tree","tree");
  tree->Branch("N",&N,"N/I");
  tree->Branch("x_position",x_position,"x_position[N]/I");
  tree->Branch("mean",mean,"mean[N]/D");
  tree->Branch("sigma",sigma,"sigma[N]/D");
  tree->Branch("error",error,"error[N]/D");
  tree->Branch("simul_mean",simul_mean,"simul_mean[N]/D");
  tree->Branch("simul_sigma",simul_sigma,"simul_sigma[N]/D");
  tree->Branch("simul_error",simul_error,"simul_error[N]/D");
  if(y_pos!=12){
    tree->Branch("total_event",tot_evt,"total_event[N]/I");
    tree->Branch("selected_event",evt_suma,"selected_event[N]/I");
  }

  for(int i=0;i<N;i++){
    x_position[i] = x_pos[i];
    mean[i] = pa_suma[i][1];
    sigma[i] = pa_suma[i][2];
    error[i] = f_suma[i]->GetParError(1);
    
    simul_mean[i] = pa_simul[i][0][1];
    simul_sigma[i] = pa_simul[i][0][2];
    simul_error[i] = f_simul[i][0]->GetParError(1);

    
  }
  tree->Fill();
  tree->Write();
  if(y_pos!=12)eff->Write();
  file_save->Close();
  
 
    
      
	
    
    
      

  
} 
