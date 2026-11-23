
#include "TBox.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <TArrayD.h>
#include <THStack.h>
#include <TStyle.h>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.h"
#ifndef __CINT__
#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.cxx"
#endif

#include "./config/ConfigReader.h"

using namespace std;

void DIC6_tof() {

  int BinPerDecade = 1500; //SETTALO

  double L = 177; //flight path for tof-E conversion

  string cmnd_filename =
      "input_files/Transmission_ratio_final.cmnd"; // for all the data (with
                                                   // selection from efficiency
                                                   // study)

  // import variables from txt file
  ConfigReader cfg(cmnd_filename);

  // Data path
  string prefix = cfg.getString("prefix");
  string suffix = cfg.getString("suffix");

  // Vectors with run numbers
  vector<int> SIN = cfg.getIntVector("SIN");
  vector<int> SOUT = cfg.getIntVector("SOUT");


  vector<int> VEC; // QUI SELEZIONI IL VETTORE DI RUN CHE VUOI ANALIZZARE
  VEC.insert(VEC.end(), SIN.begin(), SIN.end());
  VEC.insert(VEC.end(), SOUT.begin(), SOUT.end());
  std::cout << "Selected runs: ";
  for (const auto &run : VEC) {
    std::cout << run << " ";
  }
  std::cout << std::endl;




   // setting up through unordered set
  auto to_set = [](const vector<int> &v) {
    return std::unordered_set<int>(v.begin(), v.end());
  };

  const auto vec_SET = to_set(VEC);

    // define the variables corresponding to the tree leaves that will be
  // extracted
  //int detn;
  double tof;
  double tflash;
  double tflash_p;
  float amp;
  int BunchNumber;
  
  int BunchNumberPK;
  float PulseIntensity, PulseIntensityPK;
  int PSpulse, PSpulsePK;

  // define logaritmic binning for the ToF histogram
  int N_DEC_tof = 8;
  // take BinPerDecade from the input
  int n_vec_tof = (int)BinPerDecade * N_DEC_tof;
  double xbins_tof[100100];
  double step_tof = (double)N_DEC_tof / n_vec_tof;
  double value_tof = 1;
  for (int filler = 0; filler <= n_vec_tof; ++filler) {
    value_tof = step_tof * (double)filler;
    xbins_tof[filler] = (double)pow(10., value_tof);
  }


  // binning log in E da 1e-2 a 1e10
  double E_LOG_MIN = -2;  // 10^-2 eV
  double E_LOG_MAX = 10;  // 10^10 eV
  int N_DEC_E = (int)(E_LOG_MAX - E_LOG_MIN);   // 12 decadi
  int n_vec_E = BinPerDecade * N_DEC_E;
  double xbins_E[100100];
  double step_E = (double)N_DEC_E / n_vec_E;
  for (int fE = 0; fE <= n_vec_E; ++fE) {
    xbins_E[fE] = pow(10., E_LOG_MIN + step_E * fE);
  }
  

  TH1D *HSin = new TH1D("Sin", "", n_vec_tof, xbins_tof);
   HSin->Sumw2();

   TH1D* HPROVA = new TH1D("HPROVA", "HPROVA", n_vec_E, xbins_E);


  int BinPerDecade_G = 10000; //SETTALO
  int N_DEC_tof_G = 8;
  int n_vec_tof_G = BinPerDecade_G * N_DEC_tof_G;
   double xbins_tof_G[100100];
  double step_tof_G = (double)N_DEC_tof_G / n_vec_tof_G;
  double value_tof_G = 1;
  for (int filler_G = 0; filler_G <= n_vec_tof_G; ++filler_G) {
    value_tof_G = step_tof_G * (double)filler_G;
    xbins_tof_G[filler_G] = (double)pow(10., value_tof_G);
  }
  

   TH1D *HSin_flash = new TH1D("Sin_flash", "", n_vec_tof_G, xbins_tof_G);
   //HSin_flash->GetXaxis()->SetRangeUser(5.e+3, 5.e+4);

   TH2D *H2Din = new TH2D("H2Din", "", n_vec_tof, xbins_tof, 1500, 0, 5.e+3);

   TH2D *H2DE = new TH2D("H2DE", "", n_vec_E, xbins_E, 1500, 0, 5.e+3);
   //H2DE->GetXaxis()->SetRangeUser(1.e-2, 1.e+10);

   TH1D *H_amp = new TH1D("H_amp", "", 1500, 0, 5.e+3);


   //histo for tof spectra with amplitude cuts
   const double ampMin = 60;
   const double ampMax = 1000;
   TH1D *H_tof_cut = new TH1D("H_tof_cut", "", n_vec_tof, xbins_tof);
   H_tof_cut->GetXaxis()->SetTitle("ToF-Tg (ns)");
   H_tof_cut->GetYaxis()->SetTitle("Counts");
   H_tof_cut->SetTitle("DIC6 DETECTOR ToF distribution for total data set WITH AMP CUTS");

   const double ampMin_b = 100;
   const double ampMax_b = 1000;
      TH1D *H_tof_cut_b = new TH1D("H_tof_cut_b", "", n_vec_tof, xbins_tof);
   H_tof_cut_b->GetXaxis()->SetTitle("ToF-Tg (ns)");
   H_tof_cut_b->GetYaxis()->SetTitle("Counts");
   H_tof_cut_b->SetTitle("DIC6 DETECTOR ToF distribution for total data set WITH AMP CUTS");
   



   
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ANALYSIS ON VEC SET

  int n_runs_SiTOT = VEC.size(); // length of the total Sin array
  int Run_SiTOT =
      0; // variable to store the number of the current run in the loop
  Char_t R_iTOT[100]; // variable to store the whole path of the current run

  TTree *t_iTOT, *t_pkup_iTOT; // ROOT trees for the 2 detectors of interest

  // loop over all the runs of Sin
  for (int k = 0; k < n_runs_SiTOT; ++k) {
    Run_SiTOT = VEC[k];
    sprintf(R_iTOT, "%s%i%s", prefix.c_str(), Run_SiTOT, suffix.c_str());
    cout << "Now Sorting " << R_iTOT << endl;

    /*
    // check to see if the current run is included in the set of runs of the
    // various detectors
    bool check_ = vec_SET.count(Run_SiTOT) > 0;

    cout << "CHECK:  " << check_ << endl;*/

    // open the ROOT file corresponding to the current run
    TFile *f_i = TFile::Open(R_iTOT, "READ");

    if (!f_i || f_i->IsZombie() || !f_i->IsOpen()) {
      std::cerr << "ERROR: cannot open file: " << R_iTOT << "\n";
      if (f_i) {
        delete f_i;
        f_i = nullptr;
      }
      continue; // skip this run
    }


     // read PKUP tree and extract the necessary leaves
    t_pkup_iTOT = (TTree *)f_i->Get("PKUP");

    // set read only branches
    t_pkup_iTOT->SetBranchStatus("*", 0);
    t_pkup_iTOT->SetBranchStatus("tflash", 1);
    t_pkup_iTOT->SetBranchStatus("BunchNumber", 1);
    t_pkup_iTOT->SetBranchStatus("PulseIntensity", 1);
    t_pkup_iTOT->SetBranchStatus("PSpulse", 1);

    t_pkup_iTOT->SetBranchAddress("tflash", &tflash_p);
    t_pkup_iTOT->SetBranchAddress("BunchNumber", &BunchNumberPK);
    t_pkup_iTOT->SetBranchAddress("PulseIntensity", &PulseIntensityPK);
    t_pkup_iTOT->SetBranchAddress("PSpulse", &PSpulsePK);

    // extract total number of protons (pulse intensity) from PKUP

    Long64_t nentryPK =
        t_pkup_iTOT->GetEntriesFast(); // get the total entries of the tree
    if (nentryPK < 0)
      nentryPK = t_pkup_iTOT->GetEntries();


    // read DIC6 tree and extract relevant leaves
    t_iTOT = (TTree *)f_i->Get("DIC6");

    t_iTOT->SetBranchStatus("tof", 1);
     t_iTOT->SetBranchStatus("BunchNumber", 1);
     t_iTOT->SetBranchStatus("amp", 1);
     
     t_iTOT->SetBranchAddress("tof", &tof);
     t_iTOT->SetBranchAddress("tflash", &tflash);
     t_iTOT->SetBranchAddress("BunchNumber", &BunchNumber);
     t_iTOT->SetBranchAddress("amp", &amp);

      // this method allows to read the 2 trees in parallel
    t_iTOT->AddFriend(t_pkup_iTOT);
    t_pkup_iTOT->BuildIndex("BunchNumber");

     // extract the total number of entries of the tree
    Long64_t nentry = t_iTOT->GetEntriesFast();
    if (nentry < 0)
      nentry = t_iTOT->GetEntries();
    cout << "Number of Entries: " << nentry << endl;

    int counter = 0;

    // loop over the entries of the tree
    for (Long64_t i = 0; i < nentry; ++i) {

      t_iTOT->GetEntry(i);
      

       // loop over the bunches
      //
      if (BunchNumberPK == BunchNumber) {

        counter += 1;
      

        //if(check_ == true) {

	  if(tflash>=11070 && tflash<=11130 && tof-tflash > 0){
	  
	  HSin_flash->Fill(tflash);
	  
	  
    HSin->Fill(tof-tflash);
	  
	  H2Din->Fill(tof-tflash, amp);

	  H2DE->Fill(t2e_relat(tof-tflash, L, 0.0), amp);

	  HPROVA->Fill(t2e_relat(tof-tflash, L, 0.0));

    H_amp->Fill(amp);


    if (amp >= ampMin && amp <= ampMax) {
      H_tof_cut->Fill(tof - tflash);
    }

    if (amp >= ampMin_b && amp <= ampMax_b) {
      H_tof_cut_b->Fill(tof - tflash);
    }

	  }
	  
        //}
      }
        
    } // close loop on tree entries

    //cout << "Processed entries: " << counter << '\n';

  } // close loop on the runs


  HSin->GetXaxis()->SetTitle("ToF-Tg (ns)");
  HSin->GetYaxis()->SetTitle("Counts");
  HSin->SetTitle("DIC6 DETECTOR ToF distribution for total data set");
  //HSin->GetXaxis()->SetRangeUser(xmin_plot, 1e8);

  HSin_flash->GetXaxis()->SetTitle("Tgamma (ns)");
  HSin_flash->GetYaxis()->SetTitle("Counts");
  HSin_flash->SetTitle("DIC6 DETECTOR Tgamma distribution for total data set");

  H2Din->GetXaxis()->SetTitle("ToF-Tg (ns)");
  H2Din->GetYaxis()->SetTitle("Amp (channels)");
  H2Din->SetTitle("DIC6 DETECTOR 2D ToF-amp distribution for total data set");

  H2DE->GetXaxis()->SetTitle("E (eV)");
  H2DE->GetYaxis()->SetTitle("Amp (channels)");
  H2DE->SetTitle("DIC6 DETECTOR 2D En-amp distribution for total data set");

  HPROVA->GetXaxis()->SetTitle("E (eV)");
  HPROVA->GetYaxis()->SetTitle("Counts");
  HPROVA->SetTitle("DIC6 DETECTOR energy distribution for total data set");

  H_amp->GetXaxis()->SetTitle("Amp (channels)");
  H_amp->GetYaxis()->SetTitle("Counts");  
  H_amp->SetTitle("DIC6 DETECTOR amplitude distribution for total data set");

 

   TCanvas *c_Sin =
      new TCanvas("C 1D tof", "C 1D tof", 2400, 1700);
  HSin->Draw("HIST");
  c_Sin->SetLogx();
  c_Sin->SetLogy();

   TCanvas *c_Sin_gamma =
      new TCanvas("C gflash t", "C gflash t", 2400, 1700);
  HSin_flash->Draw("HIST");
  c_Sin_gamma->SetLogx();
  c_Sin_gamma->SetLogy();

  TCanvas *c_2Din =
      new TCanvas("C 2D tof", "C 2D tof", 2400, 1700);
  H2Din->Draw("COLZ");
  c_2Din->SetLogx();
  c_2Din->SetLogz();


  TCanvas *c_2D_E = new TCanvas("c 2D E","c_2D_E", 2400, 1700);
  H2DE->Draw("COLZ");
  c_2D_E->SetLogx();
  c_2D_E->SetLogz();

  TCanvas *c_E_1D = new TCanvas("c 1D E","c 1D E", 2400, 1700);
  //HPROVA->GetXaxis()->SetRangeUser(1e-2, 1e10);
  HPROVA->Draw("HIST");
  c_E_1D->SetLogx();
  c_E_1D->SetLogy();

  TCanvas *c_amp_1D = new TCanvas("c 1D amp","c 1D amp", 2400, 1700);
  H_amp->Draw("HIST");
  c_amp_1D->SetLogy();


  TCanvas *c_tof_cut = new TCanvas("c 1D tof cut a","c 1D tof cut a", 2400, 1700);
  H_tof_cut->Draw("HIST");
  TLegend *cutsLegend = new TLegend(0.62, 0.76, 0.89, 0.88);
    cutsLegend->SetBorderSize(1);
    cutsLegend->SetFillStyle(0);
    cutsLegend->SetTextSize(0.045);
    cutsLegend->AddEntry(static_cast<const TObject *>(nullptr),
                         Form("Amplitude cut: %.0f - %.0f channels", ampMin, ampMax),
                         "");
    cutsLegend->Draw();
  c_tof_cut->SetLogx();
  c_tof_cut->SetLogy();

  TCanvas *c_tof_cut_b = new TCanvas("c 1D tof cut b","c 1D tof cut b", 2400, 1700);  
  H_tof_cut_b->Draw("HIST");
  TLegend *cutsLegend_b = new TLegend(0.62, 0.76, 0.89, 0.88);
    cutsLegend_b->SetBorderSize(1);
    cutsLegend_b->SetFillStyle(0);
    cutsLegend_b->SetTextSize(0.045);
    cutsLegend_b->AddEntry(static_cast<const TObject *>(nullptr),
                         Form("Amplitude cut: %.0f - %.0f channels", ampMin_b, ampMax_b),
                         "");
    cutsLegend_b->Draw();
  c_tof_cut_b->SetLogx();
  c_tof_cut_b->SetLogy(); 

  c_2D_E->SaveAs("OUTPUT/DIC6/DIC6_2D_E_morebins.root");
  /*
  TFile *f = TFile::Open("OUTPUT/DIC6/DIC6_AllData_morebins.root", "RECREATE");
  c_Sin_gamma->Write();
  c_Sin->Write();
  c_E_1D->Write();
  c_2Din->Write();
  c_2D_E->Write();
  c_amp_1D->Write();
  
  TFile *f = TFile::Open("OUTPUT/DIC6/DIC6_TofSpectra_AmpCuts.root", "RECREATE");
  c_Sin->Write();
  c_tof_cut->Write();
  c_tof_cut_b->Write();
  */

}//close macro  
