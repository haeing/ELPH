void withoutbeam(){
  Int_t N=2;
  TFile *file_pe[N];
  TTree *data_pe[N];
  
  file_pe[0] = new TFile("../../ELPH_data/exp_data/run00079.root","read"); //Before Beam (morning)
  file_pe[1] = new TFile("../../ELPH_data/exp_data/run00333.root","read"); //After Beam (Night)

  
  Double_t ADC_pe[N][4];
  Double_t TDC_pe[N][4][16];
  Double_t ADCs_pe[N][1];
  Double_t TDCs_pe[N][1][16];
  for(int i=0;i<N;i++){
    data_pe[i] = (TTree*)file_pe[i]->Get("tree");
    data_pe[i] ->SetBranchAddress("E72BACa",ADC_pe[i]);
    data_pe[i] ->SetBranchAddress("E72BACt",TDC_pe[i]);
    data_pe[i] ->SetBranchAddress("E72BACSUMa",ADCs_pe[i]);
    data_pe[i] ->SetBranchAddress("E72BACSUMt",TDCs_pe[i]);
  }
  
  Double_t total_pe = data_pe[0]->GetEntries();

  TH1D* hist_ADC_pe[N][4];
  TH1D* hist_TDC_pe[N][4];

  //Range
  Double_t ADC_start = 50;
  Double_t ADC_end = 500;
  Double_t TDC_start = 0;
  Double_t TDC_end = 1500;
  TH1D* hist_ADCs_pe[N];
  TH1D* hist_TDCs_pe[N];
  TF1* fit_ADCs_pe[N];
  TF1* fit_TDCs_pe[N];
  TF1* fit_ADC_pe[N][4];
  TF1* fit_TDC_pe[N][4];
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      hist_ADC_pe[i][j] = new TH1D(Form("hist_ADC_pe_BAC%d_%d",j+1,i),Form("hist_ADC_pe_BAC%d_%d",j+1,i),ADC_end-ADC_start,ADC_start,ADC_end);
      hist_TDC_pe[i][j] = new TH1D(Form("hist_TDC_pe_BAC%d_%d",j+1,i),Form("hist_TDC_pe_BAC%d_%d",j+1,i),(TDC_end-TDC_start)/10,TDC_start,TDC_end);
      fit_ADC_pe[i][j] = new TF1(Form("fit_ADC_pe_BAC%d_%d",j+1,i),"gaus(0)",ADC_start,ADC_end);
      fit_TDC_pe[i][j] = new TF1(Form("fit_TDC_pe_BAC%d_%d",j+1,i),"gaus(0)",TDC_start,TDC_end);
    }
    hist_ADCs_pe[i] = new TH1D(Form("hist_ADCs_pe_%d",i),Form("hist_ADCs_pe_%d",i),ADC_end-ADC_start,ADC_start,ADC_end);
    hist_TDCs_pe[i] = new TH1D(Form("hist_TDCs_pe%d",i),Form("hist_TDCs_pe_%d",i),(TDC_end-TDC_start)/10,TDC_start,TDC_end);
    fit_ADCs_pe[i] = new TF1(Form("fit_ADCs_pe_%d",i),"gaus(0)",ADC_start,ADC_end);
    fit_TDCs_pe[i] = new TF1(Form("fit_TDCs_pe_%d",i),"gaus(0)",TDC_start,TDC_end);
  }
  for(int j=0;j<N;j++){
    for(int n=0;n<total_pe;n++){
      data_pe[j]->GetEntry(n);
      for(int i=0;i<4;i++){
	hist_ADC_pe[j][i]->Fill(ADC_pe[j][i]);
	hist_TDC_pe[j][i]->Fill(TDC_pe[j][i][0]);
      }
      hist_ADCs_pe[j]->Fill(ADCs_pe[j][0]);
      hist_TDCs_pe[j]->Fill(TDCs_pe[j][0][0]);
    }
  }
    
  Double_t parameter_ADC_pe[N][5][3];
  Double_t parameter_TDC_pe[N][5][3];

  TCanvas *c1 = new TCanvas("c1","ADC TDC w/o beam",800,650);
  c1->Divide(5,2);
  for(int i=0;i<5;i++){
    if(i<4){
      c1->cd(i+1);
      hist_ADC_pe[0][i]->Fit(fit_ADC_pe[0][i],"","",ADC_start,ADC_end);
      fit_ADC_pe[0][i]->GetParameters(parameter_ADC_pe[0][i]);
      hist_ADC_pe[1][i]->SetLineColor(kRed);
      hist_ADC_pe[1][i]->SetFillColor(kRed);
      hist_ADC_pe[1][i]->SetFillStyle(3001);
      hist_ADC_pe[1][i]->Draw("sames");

      c1->cd(5+i+1);
      hist_TDC_pe[0][i]->Fit(fit_TDC_pe[0][i],"","",TDC_start,TDC_end);
      fit_TDC_pe[0][i]->GetParameters(parameter_TDC_pe[0][i]);
      hist_TDC_pe[1][i]->SetLineColor(kRed);
      hist_TDC_pe[1][i]->SetFillColor(kRed);
      hist_TDC_pe[1][i]->SetFillStyle(3001);
      hist_TDC_pe[1][i]->Draw("sames");
    }
    if(i==4){
      c1->cd(i+1);
      hist_ADCs_pe[0]->Fit(fit_ADCs_pe[0],"","",ADC_start,ADC_end);
      fit_ADCs_pe[0]->GetParameters(parameter_ADC_pe[0][i]);
      hist_ADCs_pe[1]->SetLineColor(kRed);
      hist_ADCs_pe[1]->SetFillColor(kRed);
      hist_ADCs_pe[1]->SetFillStyle(3001);
      hist_ADCs_pe[1]->Draw("sames");

      c1->cd(5+i+1);
      hist_TDCs_pe[0]->Fit(fit_TDCs_pe[0],"","",TDC_start,TDC_end);
      fit_TDCs_pe[0]->GetParameters(parameter_TDC_pe[0][i]);
      hist_TDCs_pe[1]->SetLineColor(kRed);
      hist_TDCs_pe[1]->SetFillColor(kRed);
      hist_TDCs_pe[1]->SetFillStyle(3001);
      hist_TDCs_pe[1]->Draw("sames");
    }
  }
}
