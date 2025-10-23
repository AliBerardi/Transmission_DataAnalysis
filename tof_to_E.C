
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The first macro converts the transmission histogram from tof domain to energy domain.
//The second macro plots the cross-section as a function of energy. It also draws ORELA data.
//
//When calling the functions, you can choose the flight path (L) and the number of bins/decade of the histogram (bins).
//For example, run with:
// root -l -b -q 'tof_to_E.C(182.1, 200)'
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.h"
#ifndef __CINT__
#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.cxx"
#endif

using namespace std;


TH1D* convert_to_E(double L, int bins){

  double t0=0.0;

  TFile *file_TOTtransm = TFile::Open(Form("./OUTPUT/Transmission/Total/Transmission_total_%dbin.root", bins));



  TCanvas *c = (TCanvas*)file_TOTtransm->Get("C_Transm_final");
  //c->GetListOfPrimitives()->ls(); //check what is inside the canvas
  TH1D *h_TT = (TH1D*)c->GetListOfPrimitives()->FindObject("Total transmission");
  h_TT->GetXaxis()->SetRangeUser(0, 1e8);
  h_TT->SetTitle("Transmission");



  TH1D *h_E_TT = cnv_hst_t2e(h_TT, L, t0);

  TCanvas *c_TT_E= new TCanvas("C_Transm_final", "Canvas Transm final", 2000, 1500);
  h_E_TT->Draw("HIST");
  c_TT_E->SetLogx();

  string output_path = Form("./OUTPUT/Transmission/Total/Transmission_energy__L%.2f__%dbin.root",L, bins);
  c_TT_E->SaveAs(output_path.c_str());

  return h_E_TT;

}

void compute_CS(double L, int bins){

  TH1D *h_E_TT = convert_to_E(L, bins);


  TH1D* h_CS = (TH1D*)h_E_TT->Clone("h_CS");
  
  h_CS->GetXaxis()->SetRangeUser(0,1e7);
  
  double areal_density=0.05174;
  for (Int_t i = 1; i <= h_CS->GetNbinsX(); i++)  {
      Double_t tra = h_CS->GetBinContent(i);
      Double_t tra_err = h_CS->GetBinError(i)/h_CS->GetBinContent(i);
      Double_t lntra =  -TMath::Log(tra)/areal_density;
      if(tra>0.0){ 
      h_CS->SetBinContent(i, lntra);
      h_CS->SetBinError(i, lntra*tra_err);
      }
      if(tra<=0.0){ 
      h_CS->SetBinContent(i, 0.0);
      h_CS->SetBinError(i, 0.1);
      }
	}

  

TGraphErrors * gJH = new TGraphErrors("./input_files/JHarvey.dat","%lg %lg %lg");


//convert the histogram into a graph to have a line that connects the dots
 int n = h_CS->GetNbinsX();
 TGraphErrors *ge = new TGraphErrors(n);
for (int i = 1; i <= n; ++i) {
    double x  = h_CS->GetBinCenter(i);
    double y  = h_CS->GetBinContent(i);
    double ex = 0.0;
    double ey = h_CS->GetBinError(i);
    ge->SetPoint(i-1, x, y);
    ge->SetPointError(i-1, ex, ey);
}



  
  TCanvas* Ccs=new TCanvas("Ccs","Cross section",2000,1500);
  Ccs->cd();
  h_CS->SetTitle("(n, natCu) total cross section");
  h_CS->GetYaxis()->SetTitle("Cross section (barn)");
  h_CS->GetXaxis()->SetTitle("Neutron energy (eV)");
  h_CS->SetLineColor(kBlue);

  ge->SetTitle("(n, natCu) total cross section");
  ge->GetYaxis()->SetTitle("Cross section (barn)");
  ge->GetXaxis()->SetTitle("Neutron energy (eV)");
  ge->SetLineColor(kBlue);
  //ge->Draw("PL"); // punti + linea + assi
  
  
  
  h_CS->Draw("e1");

  gPad->SetLogx(1);gPad->SetLogy(1);

  gJH->SetLineColor(kGreen+3);
  gJH->Draw("PLsame");
   Ccs->Update();
   

 

   TLegend* leg = new TLegend(.6,.7,.9,.9);  leg->SetTextSize(0.04);   leg->SetBorderSize(0);  
  leg->SetFillColor(0); 
  leg->AddEntry(h_CS,"n_TOF","l");
  leg->AddEntry(gJH,"ORELA","l");
  leg->Draw();
  Ccs->Update();
  
  Ccs->SaveAs(Form("./OUTPUT/Transmission/Total/CrossSection__L%.2f__%dbin_hist.root",L, bins));

 
  
}


void tof_to_E(double L, int bins){

  //double L = 181.9;

  //convert_to_E(L, bins);
  compute_CS(L, bins);
  

}



