Int_t y_pos = 0;

void trigger(){

  gStyle->SetOptStat(0);
  TFile *file_pe;
  if(y_pos == 12 ||y_pos == 0||y_pos ==-10)file_pe = new TFile("../../ELPH_data/exp_data/run00333.root","read");
  else if(y_pos ==-35 || y_pos ==-23 ||y_pos ==-11)file_pe = new TFile("../../ELPH_data/exp_data/run00079.root","read");
  
  Int_t N;  //x position
  Int_t A = 10;
  Int_t factor_a[A]; //adc cut
  for(int i=0;i<A;i++)factor_a[i] = 3+i;

  if(y_pos == 0)N=7;
  else if(y_pos ==-10)N=3;
  else if(y_pos ==-11)N=7;
  else if(y_pos ==-23)N=7;
  else if(y_pos ==-35)N=11;
  else if(y_pos == 12)N=6;
  
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
    x_pos[0] = -6;
    x_pos[1] = -4;
    x_pos[2] = -2;
    x_pos[3] = 0;
    x_pos[4] = 2;
    x_pos[5] = 4;
  }

  TFile *file_po[N];
  TTree *data_po[N];

  Double_t Ta[N][4][1];
  Double_t Tt[N][4][1][16];

  Double_t tot[N];

  TH1D *hist_Ta[N][4];
  TH1D *hist_Tt[N][4];

  TH1D *hist_Ta_cut[A][N][4];
  TH1D *hist_Tt_cut[A][N][4];

  TF1 *f_a[N][4];
  TF1 *f_t[N][4];

  if(y_pos==0){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00313.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00314.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00315.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00297.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00316.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00317.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00319.root","read");

  }


  else if(y_pos==-10){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00300.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00299.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00301.root","read");

    
  }
  

  else if(y_pos==-11){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00065.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00064.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00062.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00063.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00061.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00060.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00059.root","read");

    
  }
  

  else if(y_pos==-23){
    file_po[0] = new TFile("../../ELPH_data/exp_data/run00048.root","read");
    file_po[1] = new TFile("../../ELPH_data/exp_data/run00052.root","read");
    file_po[2] = new TFile("../../ELPH_data/exp_data/run00053.root","read");
    file_po[3] = new TFile("../../ELPH_data/exp_data/run00045.root","read");
    file_po[4] = new TFile("../../ELPH_data/exp_data/run00054.root","read");
    file_po[5] = new TFile("../../ELPH_data/exp_data/run00055.root","read");
    file_po[6] = new TFile("../../ELPH_data/exp_data/run00056.root","read");

    
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


  }

  for(int i=0;i<N;i++){
    data_po[i] = (TTree*)file_po[i]->Get("tree");
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
      hist_Ta[i][j] = new TH1D(Form("hist_T%da_%dcm",j+4,x_pos[i]),Form("hist_T%da_%dcm",j+4,x_pos[i]),450,0,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,x_pos[i]),Form("hist_T%dt_%dcm",j+4,x_pos[i]),250,600,850);
      for(int k=0;k<A;k++){
	hist_Ta_cut[k][i][j] = new TH1D(Form("hist_T%da_%dcm_factor%d",j+4,x_pos[i],factor_a[k]),Form("hist_T%da_%dcm_factor%d",j+4,x_pos[i],factor_a[k]),450,0,2000);
	hist_Tt_cut[k][i][j] = new TH1D(Form("hist_T%dt_%dcm_factor%d",j+4,x_pos[i],factor_a[k]),Form("hist_T%dt_%dcm_factor%d",j+4,x_pos[i],factor_a[k]),250,600,850);
      }
      
      
      
      f_a[i][j] = new TF1(Form("f_a%da_%dcm",j+4,x_pos[i]),"landau(0)",200,2000);
      f_t[i][j] = new TF1(Form("f_t%dt_%dcm",j+4,x_pos[i]),"gaus(0)",650,730);
      
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

  Double_t pa_a[N][4][3];
  Double_t pa_t[N][4][3];
  
  Int_t pass[A];
  

  Int_t evt[A][N];
  for(int i=0;i<A;i++){
    for(int j=0;j<N;j++){
      evt[i][j] = 0;
    }
  }
  
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      hist_Ta[i][j]->Fit(f_a[i][j],"Q0","",200,2000);
      f_a[i][j]->GetParameters(pa_a[i][j]);
      cout<<pa_a[i][j][2]<<endl;
      hist_Tt[i][j]->Fit(f_t[i][j],"Q0","",650,730);
      f_t[i][j]->GetParameters(pa_t[i][j]);
    }
  }

  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      for(int q=0;q<A;q++)pass[q] = 0;
      for(int j=0;j<4;j++){
	for(int k=0;k<A;k++){
	  if(Tt[i][j][0][0]>pa_t[i][j][1]-5*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+5*pa_t[i][j][2]){
	    if(Ta[i][j][0]>pa_a[i][j][1]-factor_a[k]*pa_a[i][j][2]&&Ta[i][j][0]<3840){
	    
	      pass[k]+=1;
	    }
	  }
	    

	}
      
	
      }
      for(int k=0;k<A;k++){
	if(pass[k]==4){
	  evt[k][i]+=1;
	  for(int j=0;j<4;j++){
	    
	    hist_Ta_cut[k][i][j]->Fill(Ta[i][j][0]);
	    hist_Tt_cut[k][i][j]->Fill(Tt[i][j][0][0]);
	    
	  }
	  
	}
      }
    
    }
  }


  TCanvas *c1[A];
  for(int k=0;k<A;k++){
    c1[k] = new TCanvas(Form("c1%d",k),Form("Trigger counters ADC histogram factor %d",factor_a[k]),800,650);
    //c1[k]->Divide(4,N);
  
    for(int i=0;i<1;i++){
      for(int j=0;j<1;j++){
	c1[k]->cd(4*i+j+1);
	gPad->SetLogy();
	hist_Ta[i][j]->SetTitle("ADC;ADC [Ch.];n");
	hist_Ta[i][j]->SetLineColor(kBlack);
	hist_Ta[i][j]->Draw();

	hist_Ta_cut[k][i][j]->SetLineColor(kRed);
	hist_Ta_cut[k][i][j]->SetFillColor(kRed);
	hist_Ta_cut[k][i][j]->SetFillStyle(3001);
	hist_Ta_cut[k][i][j]->Draw("sames");
      
      }
    }
  }

  TCanvas *c2[A];
  for(int k=0;k<A;k++){
    c2[k] = new TCanvas(Form("c2%d",k),Form("Trigger counters TDC histogram factor %d",factor_a[k]),800,650);
    //c2[k]->Divide(4,N);
    for(int i=0;i<1;i++){
      for(int j=0;j<1;j++){
	c2[k]->cd(4*i+j+1);
	gPad->SetLogy();
	hist_Tt[i][j]->SetTitle("TDC;TDC [Ch.];n");
	hist_Tt[i][j]->SetLineColor(kBlack);
	hist_Tt[i][j]->GetListOfFunctions()->Remove(f_t[i][j]);
	hist_Tt[i][j]->Draw();
	hist_Tt_cut[k][i][j]->SetLineColor(kRed);
	hist_Tt_cut[k][i][j]->SetFillColor(kRed);
	hist_Tt_cut[k][i][j]->SetFillStyle(3001);
	hist_Tt_cut[k][i][j]->Draw("sames");
	
      }
    }
  }


  TGraph *eff[A];
  TGraph *ef;
  Double_t efficiency[A][N];
  Double_t x_position[N];
  Double_t fact[A];
  Double_t ef_f[A];
  TMultiGraph *mg = new TMultiGraph();
  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  
  for(int i=0;i<A;i++){
    for(int j=0;j<N;j++){
      x_position[j] = x_pos[j]*10.0;
      efficiency[i][j] = evt[i][j]/tot[j]*100;
    }
    eff[i] = new TGraph(N,x_position,efficiency[i]);
    eff[i]->SetMarkerStyle(24);
    eff[i]->SetMarkerSize(1);
    eff[i]->SetMarkerColor(i+1);
    mg->Add(eff[i]);
    le->AddEntry(eff[i],Form("factor %d",factor_a[i]));
  }

  TCanvas *c3 = new TCanvas("c3","c3",800,650);
  c3->cd();
  mg->Draw("AP");
  le->Draw();
  
  for(int i=0;i<A;i++){
    fact[i] = factor_a[i]*pa_a[i][1][2];
    ef_f[i] = evt[i][3]/tot[3]*100;
  }
  ef = new TGraph(A,fact,ef_f);
    
  TCanvas *c4 = new TCanvas("c4","c4",800,650);
  c4->cd();
  ef->SetMarkerStyle(24);
  ef->SetMarkerSize(1);
  ef->Draw("AP");
  
  

  
      
      
  


}
