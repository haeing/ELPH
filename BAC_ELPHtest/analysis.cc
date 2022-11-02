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
  TCanvas *c_pe[5];

  for(int i=0;i<5;i++){
    c_pe[i] = new TCanvas(Form("c_pe_BAC%d",i+1),Form("c_pe_BAC%d",i+1),800,650);
    if(i<4){
    c_pe[i]->cd();
    hist_pe[i]->Fit(fit_pe[i],"Q","",50,500);
    fit_pe[i]->GetParameters(parameter_pe[i]);
    }
    if(i==4){
    c_pe[4]->cd();
    hist_sum->Fit(fit_pe_sum,"Q","",50,500);
    fit_pe_sum->GetParameters(parameter_pe[4]);
    }
  }

  Int_t N=5;
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
  TH1D *hist_Suma[N];
  TH1D *hist_Sumt[N];
  TF1 *f_a[N][4];
  TF1 *f_Tt[N][4];
  TF1 *f_sumt[N];
  TF1 *f_suma[N];
  
  file_po[0] = new TFile("data/run00314.root","read");
  file_po[1] = new TFile("data/run00315.root","read");
  file_po[2] = new TFile("data/run00296.root","read");
  file_po[3] = new TFile("data/run00316.root","read");
  file_po[4] = new TFile("data/run00317.root","read");
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
      hist_Ta[i][j] = new TH1D(Form("hist_T%da_%dcm",j+4,-4+2*i),Form("hist_T%da_%dcm",j+4,-4+2*i),450,200,2000);
      hist_Tt[i][j] = new TH1D(Form("hist_T%dt_%dcm",j+4,-4+2*i),Form("hist_T%dt_%dcm",j+4,-4+2*i),80,650,730);
      f_a[i][j] = new TF1(Form("f_Tt%da_%dcm",j+4,-4+2*i),"landau(0)",200,2000);
      f_Tt[i][j] = new TF1(Form("f_Tt%dt_%dcm",j+4,-4+2*i),"gaus(0)",650,730);
			  
			  
    }
    hist_Suma[i] = new TH1D(Form("hist_Suma_%dcm",-4+2*i),Form("hist_Suma_%dcm",-4+2*i),100,0,60);
    hist_Sumt[i] = new TH1D(Form("hist_Sumt_%dcm",-4+2*i),Form("hist_Sumt_%dcm",-4+2*i),200,650,850);
    f_sumt[i] = new TF1(Form("f_sumt_%dcm",-4+2*i),"gaus(0)",650,850);
    f_suma[i] = new TF1(Form("f_suma_%dcm",-4+2*i),"landau(0)",200,2000);
    

    
  }


  
  for(int i=0;i<N;i++){
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      for(int j=0;j<4;j++){
	hist_Ta[i][j]->Fill(Ta[i][j][0]);
	hist_Tt[i][j]->Fill(Tt[i][j][0][0]);

      }
      hist_Sumt[i]->Fill(TDCs[i][0][0]);

    }
  }

  Double_t pa_a[N][4][3];
  Double_t pa_t[N][4][3];
  Double_t pa_sumt[N][3];
  Double_t pa_suma[N][3];
  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  c1->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c1->cd(4*i+j+1);
      hist_Ta[i][j]->Fit(f_a[i][j],"Q","",200,2000);
      f_a[i][j]->GetParameters(pa_a[i][j]);
    }
  }

  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  c2->Divide(4,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<4;j++){
      c2->cd(4*i+j+1);
      hist_Tt[i][j]->Fit(f_Tt[i][j],"Q","",650,730);
      f_Tt[i][j]->GetParameters(pa_t[i][j]);
    }
  }

  TCanvas *c3 = new TCanvas("c3","c3",800,650);
  c3->Divide(N);
  for(int i=0;i<N;i++){
    c3->cd(i+1);
    hist_Sumt[i]->Fit(f_sumt[i],"Q","",650,860);
    f_sumt[i]->GetParameters(pa_sumt[i]);
  }

  TFile *file3[N];
  TTree *data3[N];
  Double_t ADC_sum[N];
  for(int i=0;i<N;i++){
    file3[i] = new TFile(Form("data/sum_%dcm.root",-4+2*i),"recreate");
    data3[i] = new TTree("tree","tree");
    data3[i]->Branch("ADC",&ADC_sum[i],"ADC/D");  //Branch(branchname, &value);
  }

  
  

  Int_t pass;
  Double_t one_photon= (15.4939+15.753+16.1096+16.0168)*0.9/4;
  for(int i=0;i<N;i++){
    file3[i]->cd();
    for(int n=0;n<tot[i];n++){
      data_po[i]->GetEntry(n);
      pass = 0;
	for(int j=0;j<4;j++){
	  if(Ta[i][j][0]>pa_a[i][j][1]-2*pa_a[i][j][2]&&Ta[i][j][0]<pa_a[i][j][1]+20*pa_a[i][j][2]){
	    if(Tt[i][j][0][0]>pa_t[i][j][1]-3*pa_t[i][j][2]&&Tt[i][j][0][0]<pa_t[i][j][1]+3*pa_t[i][j][2]){
	      pass+=1;
	    }
	  }
	}
	if(pass==4){
	  if(TDCs[i][0][0]>pa_sumt[i][1]-3*pa_sumt[i][2]&&TDCs[i][0][0]<pa_sumt[i][1]+3*pa_sumt[i][2]){
	    //hist_Suma[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/15.8);
	    hist_Suma[i]->Fill((ADCs[i][0]-parameter_pe[4][1])/one_photon);

	    //ADC_sum = (ADCs[i][0]-parameter_pe[4][1])/15.8;
	    ADC_sum[i] = (ADCs[i][0]-parameter_pe[4][1])/one_photon;
	    
	    data3[i]->Fill();
	  }
	
	  
	}
	//pass=0;
	  
	
    }
    data3[i]->Write();
    file3[i]->Close();
  }
  

  

  /*
  TLegend *lea = new TLegend(0.8,0.5,0.48,0.6);
  lea->AddEntry(hist_Suma[0],"(X,Y) = (-40,0)");
  lea->AddEntry(hist_Suma[1],"(X,Y) = (-20,0)");
  lea->AddEntry(hist_Suma[2],"(X,Y) = (0,0)");
  lea->AddEntry(hist_Suma[3],"(X,Y) = (20,0)");
  lea->AddEntry(hist_Suma[4],"(X,Y) = (40,0)");
  */

  TCanvas *c4 = new TCanvas("c4","c4",800,650);
  c4->Divide(N);
  Double_t xpos[N];
  Double_t x_error[N];
  Double_t npe_pos[N];
  Double_t npe_error[N];
  Double_t factor = 1.;
  for(int i=0;i<N;i++){
    c4->cd(i+1);
    //hist_Suma[i]->Scale(factor,"width");
    hist_Suma[i]->Fit(f_suma[i],"Q","",0,100);
    f_suma[i]->GetParameters(pa_suma[i]);
    
    xpos[i] = -40+i*20;
    x_error[i] = 0.5;
    npe_pos[i] = pa_suma[i][1];
    npe_error[i] = f_suma[i]->GetParError(1);
  }

  //lea->Draw();

  //Simulation - position dependence
  TFile *file_simul[N];
  TTree *data_simul[N];
  Int_t simul_npe[N];
  TH1D *hist_simul[N];
  TF1 *fit_simul[N];
  
  file_simul[0] = new TFile("../../simul_data/221024_ELPH_m40.root","read");
  file_simul[1] = new TFile("../../simul_data/221024_ELPH_m20.root","read");
  file_simul[2] = new TFile("../../simul_data/221024_ELPH_50_re.root","read");
  file_simul[3] = new TFile("../../simul_data/221024_ELPH_20.root","read");
  file_simul[4] = new TFile("../../simul_data/221024_ELPH_40.root","read");

  for(int i=0;i<N;i++){
    data_simul[i] = (TTree*)file_simul[i]->Get("tree");
    data_simul[i] ->SetBranchAddress("nhMppc",&simul_npe[i]);
    hist_simul[i] = new TH1D(Form("hist_simul%d",i),Form("hist_simul%d",i),60,0,60);
    fit_simul[i] = new TF1(Form("fit_simul%d",i),"landau(0)",0,60);
  }
    
  for(int n=0;n<999;n++){
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
    npe_simul_error[i] = fit_simul[i]->GetParError(1);
  }

  

  TGraphErrors *npe0 = new TGraphErrors(N,xpos,npe_pos,x_error,npe_error);
  
  TGraphErrors *npe1 = new TGraphErrors(N,xpos,npe_simul_pos,x_error,npe_simul_error);
  npe0->SetMarkerStyle(24);
  npe0->SetMarkerColor(1);
  npe0->SetMarkerSize(1);
  npe1->SetMarkerStyle(24);
  npe1->SetMarkerColor(2);
  npe1->SetMarkerSize(1);

  TMultiGraph *mg = new TMultiGraph();
  TLegend *le_s = new TLegend(0.8,0.5,0.48,0.6);
  mg->Add(npe0);
  mg->Add(npe1);
  le_s->AddEntry(npe0,"ELPH test");
  le_s->AddEntry(npe1,"Simulation");
  

  TCanvas *c10 = new TCanvas("c10","c10",800,650);
  c10->cd();
  mg->SetTitle("N_{p.e.} for various X position;X [mm];N_{p.e.}");
  mg->Draw("AP");
  le_s->Draw();


  

  TFile *file_com[3];
  TTree *data_com[3];
  Double_t npe[1];
  Int_t nhMppc[2];
  TH1D *hist_com[3];
  file_com[0] = new TFile("data/sum_0cm.root","read");
  file_com[1] = new TFile("../../simul_data/221024_ELPH_50_re.root","read");
  file_com[2] = new TFile("../../BACSimul/build/221024_E72_50.root","read");
  for(int i=0;i<3;i++){
    data_com[i] = (TTree*)file_com[i]->Get("tree");
    hist_com[i] = new TH1D(Form("hist_com%d",i),Form("hist_com%d",i),100,0,100);
  }
  
  data_com[0] ->SetBranchAddress("ADC",&npe[0]);
  data_com[1] ->SetBranchAddress("nhMppc",&nhMppc[0]);
  data_com[2] ->SetBranchAddress("nhMppc",&nhMppc[1]);


  
  for(int n=0;n<999;n++){
    for(int i=0;i<3;i++){
      data_com[i]->GetEntry(n);
      
    }
    hist_com[0]->Fill(npe[0]);
    hist_com[1]->Fill(nhMppc[0]);
    hist_com[2]->Fill(nhMppc[1]);
    hist_com[0]->SetLineColor(kBlack);
    hist_com[1]->SetLineColor(kRed);
    hist_com[2]->SetLineColor(kBlue);
    
  }

  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  le->AddEntry(hist_com[0],"ELPH Experiment");
  le->AddEntry(hist_com[1],"Simulation for ELPH");
  //le->AddEntry(hist_com[2],"Simulation for E72");
  

  TCanvas *c5 = new TCanvas("c5","c5",800,650);
  c5->cd();
  hist_com[0]->SetTitle(";N_{p.e.};n");
  hist_com[0]->Draw();
  hist_com[1]->Draw("sames");
  //hist_com[2]->Draw("sames");
  le->Draw();

  

  

  
      
      
      
  
  
  

  
     
}
