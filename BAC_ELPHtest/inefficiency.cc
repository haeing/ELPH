void inefficiency(){

  

  Int_t N = 3;
  
  TFile *file[N];
  TTree *tree[N];

  Double_t ADCs[N][1];
  Double_t TDCs[N][1][16];

  Double_t ADCi[N][4];
  Double_t TDCi[N][4][16];

  Double_t thre[N];
  
  thre[0] = 20.0;
  thre[1] = 80.0;
  thre[2] = 100.0;
  
  file[0] = new TFile("../../ELPH_data/exp_data/run00333.root","read");  //sum 60 mV
  file[1] = new TFile("../../ELPH_data/exp_data/run00079.root","read");  //sum 80 mV
  file[2] = new TFile("../../ELPH_data/exp_data/run00004.root","read");  //sum 100 mV
  

  
  Bool_t bPassed = 0;


  


  TEfficiency *ineff[4];
  for(int i=0;i<4;i++){
    ineff[i]  = new TEfficiency(Form("ineff%d",i),"Inefficiency;Threshold [mV];Inefficiency",100,10,110);
  }
  
  double total[N];
  for(int i=0;i<N;i++){
    tree[i] = (TTree*)file[i]->Get("tree");
    total[i] = tree[i]->GetEntries();
    tree[i]->SetBranchAddress("E72BACSUMa",ADCs[i]);
    tree[i]->SetBranchAddress("E72BACSUMt",TDCs[i]);

    tree[i]->SetBranchAddress("E72BACa",ADCi[i]);
    tree[i]->SetBranchAddress("E72BACt",TDCi[i]);
  }

  for(int i=0;i<N;i++){
    for(int n=0;n<total[i];n++){
      tree[i]->GetEntry(n);
      for(int k=0;k<4;k++){
	bPassed = 0;
	
	for(int j=0;j<16;j++)if(TDCi[i][k][j]>0&&TDCi[i][k][j]<1500)bPassed = 1;
      
	ineff[k]->Fill(bPassed,thre[i]);
      }
    }
  }

  //Efficiency
  TFile *file_eff = new TFile("parameter_9999.root","read");
  TEfficiency *eff = (TEfficiency*)file_eff->Get("eff");
  

  TLegend *le_e = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<4;i++){
    ineff[i]->SetMarkerStyle(3);
    ineff[i]->SetMarkerSize(1.5);
    ineff[i]->SetMarkerColor(i+1);
    ineff[i]->SetLineColor(i+1);
    le_e->AddEntry(ineff[i],Form("BAC Ind. %d",i));
  }
  //ineff->GetYaxis()->SetRange(0,1);

  eff->SetMarkerStyle(3);
  eff->SetMarkerSize(1.5);
  eff->SetMarkerColor(kRed);
  eff->SetLineColor(kRed);

  
  //le_e->AddEntry(eff,"Efficiency");

  
  
  TCanvas* c1 = new TCanvas("c1","c1",800,650);
  c1->cd();
  gPad->SetLogy();
  ineff[0]->Draw("AP");
  ineff[1]->Draw("same");
  ineff[2]->Draw("same");
  ineff[3]->Draw("same");
  //eff->SetTitle("Efficiency;Threshold [mV];efficiency");
  //eff->Draw("AP");
  le_e->Draw();
  


  
  }
  
    
    


