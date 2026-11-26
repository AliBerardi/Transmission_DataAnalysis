#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.h"
#ifndef __CINT__
#include "/eos/user/a/alberard/TTOFSort/TTOFSort_nTOF.cxx"
#endif

#include "./config/ConfigReader.h"
#include <algorithm>
#include <cmath>
#include <limits>

// TAGLI DI AMP E TOF
const double ampMin = 200;
const double ampMax = 260;
const double tofMin = 1.e5;
const double tofMax = 1.e8;

void correlation_dic6()
{

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

  vector<int> ALL_DATA;
  ALL_DATA.insert(ALL_DATA.end(), SIN.begin(), SIN.end());
  ALL_DATA.insert(ALL_DATA.end(), SOUT.begin(), SOUT.end());

  // setting up through unordered set
  auto to_set = [](const vector<int> &v)
  {
    return std::unordered_set<int>(v.begin(), v.end());
  };

  const auto vec_ALL_DATA = to_set(ALL_DATA);

  // define the variables corresponding to the tree leaves that will be extracted
  double tof;
  double tflash;
  float amp;
  int BunchNumber;

  int BunchNumberPK;
  float PulseIntensity, PulseIntensityPK;
  double tflash_p;
  int PSpulse, PSpulsePK;

  TTree *t_DIC6, *t_PKUP;

  int n_runs = ALL_DATA.size();
  int Run_i = 0;       // variable to store the number of the current run in the loop
  Char_t path_ri[100]; // variable to store the whole path of the current run

  // define logaritmic binning for the ToF histogram
  int BinPerDecade = 500;
  int N_DEC_tof = 8;
  // take BinPerDecade from the input
  int n_vec_tof = (int)BinPerDecade * N_DEC_tof;
  double xbins_tof[100100];
  double step_tof = (double)N_DEC_tof / n_vec_tof;
  double value_tof = 1;
  for (int filler = 0; filler <= n_vec_tof; ++filler)
  {
    value_tof = step_tof * (double)filler;
    xbins_tof[filler] = (double)pow(10., value_tof);
  }

  TH2D *hAmpTof = new TH2D("hAmpTof", "hAmpTof", n_vec_tof, xbins_tof, 500, 0, 5.e+3);
  hAmpTof->Sumw2();

  TGraphErrors *graph_A1 = new TGraphErrors();
  graph_A1->SetName("graph_A1");
  graph_A1->SetTitle(";Proton Intensity; DIC6 counts");
  graph_A1->SetMarkerStyle(20);

  TGraphErrors *graph_2 = new TGraphErrors();
  graph_2->SetTitle(";DIC6 counts; DIC6 counts/Proton Intensity");
  graph_2->SetMarkerStyle(20);

  TGraphErrors *graph_3 = new TGraphErrors();
  graph_3->SetTitle(";Proton Intensity; DIC6 counts/Proton Intensity");
  graph_3->SetMarkerStyle(20);

  
  // loop over all the runs of DIC6
  for (int k = 0; k < n_runs; ++k)
  {

    float PI_pkup = 0;          // variable to store the total number of protons in the current run
    Long64_t counts_DIC_A1 = 0; // variabile per storare il numero di conteggi del diamond nel range richiesto

    hAmpTof->Reset(); // reset histogram for each run

    Run_i = ALL_DATA[k];
    sprintf(path_ri, "%s%i%s", prefix.c_str(), Run_i, suffix.c_str());
    cout << "Now Sorting " << path_ri << endl;

    // open the ROOT file corresponding to the current run
    TFile *f_i = TFile::Open(path_ri, "READ");
    if (!f_i || f_i->IsZombie() || !f_i->IsOpen())
    {
      std::cerr << "ERROR: cannot open file: " << path_ri << "\n";
      if (f_i)
      {
        delete f_i;
        f_i = nullptr;
      }
      continue; // skip this run
    }

    // read PKUP tree and extract the necessary leaves
    t_PKUP = (TTree *)f_i->Get("PKUP");

    // set read only branches
    t_PKUP->SetBranchStatus("*", 0);
    t_PKUP->SetBranchStatus("tflash", 1);
    t_PKUP->SetBranchStatus("BunchNumber", 1);
    t_PKUP->SetBranchStatus("PulseIntensity", 1);
    t_PKUP->SetBranchStatus("PSpulse", 1);

    t_PKUP->SetBranchAddress("tflash", &tflash_p);
    t_PKUP->SetBranchAddress("BunchNumber", &BunchNumberPK);
    t_PKUP->SetBranchAddress("PulseIntensity", &PulseIntensityPK);
    t_PKUP->SetBranchAddress("PSpulse", &PSpulsePK);

    // extract entries of PKUP
    Long64_t nentryPK =
        t_PKUP->GetEntriesFast(); // get the total entries of the tree
    if (nentryPK < 0)
      nentryPK = t_PKUP->GetEntries();

    // read DIC6 tree and extract relevant leaves
    t_DIC6 = (TTree *)f_i->Get("DIC6");

    t_DIC6->SetBranchStatus("*", 0);
    t_DIC6->SetBranchStatus("tof", 1);
    t_DIC6->SetBranchStatus("tflash", 1);
    t_DIC6->SetBranchStatus("BunchNumber", 1);
    t_DIC6->SetBranchStatus("amp", 1);
    t_DIC6->SetBranchStatus("BunchNumber", 1);
    t_DIC6->SetBranchStatus("PulseIntensity", 1);
    t_DIC6->SetBranchStatus("PSpulse", 1);

    t_DIC6->SetBranchAddress("tof", &tof);
    t_DIC6->SetBranchAddress("tflash", &tflash);
    t_DIC6->SetBranchAddress("amp", &amp);
    t_DIC6->SetBranchAddress("BunchNumber", &BunchNumber);
    t_DIC6->SetBranchAddress("PulseIntensity", &PulseIntensity);
    t_DIC6->SetBranchAddress("PSpulse", &PSpulse);

    // extract entries of DIC6
    Long64_t nentry_DIC6 =
        t_DIC6->GetEntriesFast(); // get the total entries of the tree
    if (nentry_DIC6 < 0)
      nentry_DIC6 = t_DIC6->GetEntries();

    std::cout << "Number of Entries in PKUP: " << nentryPK << std::endl;
    std::cout << "Number of Entries in DIC6: " << nentry_DIC6 << std::endl;

    // this method allows to read the 2 trees in parallel
    t_DIC6->AddFriend(t_PKUP);
    t_PKUP->BuildIndex("BunchNumber");

    // loop over all the entries of pkup to get total number of protons in the current run
    for (Long64_t iP = 0; iP < nentryPK; ++iP)
    {
      t_PKUP->GetEntry(iP);
      PI_pkup += PulseIntensityPK;
    }

    std::cout << "Total number of protons in current run " << PI_pkup << std::endl;

    // filtra per il numero di bunch e fai t-tgamma per avere un confronto migliore con le regioni energetiche
    //->GetEntries

    // loop over the entries of the DIC6 tree
    int counter = 0; // counter per controllare di processare tutte le entries
    for (Long64_t i = 0; i < nentry_DIC6; ++i)
    {

      t_DIC6->GetEntry(i);

      // loop over the bunches
      if (BunchNumberPK == BunchNumber)
      { // loop over bunches

        counter += 1;

        if (tflash >= 11070 && tflash <= 11130 && tof - tflash >= tofMin && tof - tflash <= tofMax && amp >= ampMin && amp <= ampMax)
        {

          hAmpTof->Fill(tof - tflash, amp);
        }

      } // close loop over bunches

    } // close loop on tree entries

    cout << "Processed entries: " << counter << '\n';

    double dic_counts = hAmpTof->GetEntries();

    std::cout << "Counts of DIC6 detetcor in range tof:[" << tofMin << ", " << tofMax << "] and amp:[" << ampMin << ", " << ampMax << "] for current run: " << dic_counts << std::endl;

    if (PI_pkup > 0)
    {
      double ratio = static_cast<double>(dic_counts) / (PI_pkup);
      std::cout << "\n";
      std::cout << "ratio: " << ratio;
      std::cout << "\n";
      graph_A1->SetPoint(graph_A1->GetN(), PI_pkup, dic_counts);
      graph_2->SetPoint(graph_2->GetN(), dic_counts, ratio);
      graph_3->SetPoint(graph_3->GetN(), PI_pkup, ratio);
      // std::cout << "DIC6 counts/N protons PKUP for A1: " << static_cast<double>(counts_DIC_A1) / PI_pkup << std::endl;
      // std::cout << "" << std::endl;

      graph_A1->SetPointError(graph_A1->GetN() - 1, sqrt(PI_pkup), sqrt(dic_counts));
      graph_2->SetPointError(graph_2->GetN() - 1, 0, ratio * sqrt(1. / dic_counts + 1. / PI_pkup));
      graph_3->SetPointError(graph_3->GetN() - 1, 0, ratio * sqrt(1. / dic_counts + 1. / PI_pkup));
    }
    else
    {
      std::cerr << "WARNING: PI_PKUP is zero for run " << Run_i
                << "; skipping graph point." << std::endl;
    }

    // TCanvas *c1 = new TCanvas("c1", "c1", 2400, 1700);
    // c1->SetLogx();
    // c1->SetLogz();
    // hAmpTof->Draw("colz");
    // c1->SaveAs(Form("OUTPUT/DIC6/counts_DIC6_run_121%d.root", Run_i));

    // PROVA A DISEGNARE DELLE RIGHE (O ZONE COLORATE) IN CORRISPONDENZA DELLE TUE REGIONI DI TAGLIO

  } // close loop on the runs

  //graph_A1->Fit();

  TCanvas *cRatio = new TCanvas("DIC_PI", "DIC6 counts per PI_PKUP", 1200, 800);
  graph_A1->Draw("AP");

  TLegend *cutsLegend = new TLegend(0.62, 0.76, 0.89, 0.88);
  cutsLegend->SetBorderSize(0);
  cutsLegend->SetFillStyle(0);
  cutsLegend->SetTextSize(0.045);
  cutsLegend->AddEntry(static_cast<const TObject *>(nullptr),
                       Form("Amplitude cut: %.0f - %.0f channels", ampMin, ampMax),
                       "");
  cutsLegend->AddEntry(static_cast<const TObject *>(nullptr),
                       Form("TOF cut: %.1e - %.1e ns", tofMin, tofMax),
                       "");
  cutsLegend->Draw();

  TCanvas *c2 = new TCanvas("DIC_PI/DIC", "DIC_PI/DIC", 1200, 800);
  graph_2->Draw("AP");
  cutsLegend->Draw();

  TCanvas *c3 = new TCanvas("DIC_PI/PI", "DIC_PI/PI", 1200, 800);
  graph_3->Draw("AP");
  cutsLegend->Draw();

  TFile *f = TFile::Open("OUTPUT/DIC6/correlations/correlation_DIC6_PI.root", "RECREATE");
  cRatio->Write();
  c2->Write();
  c3->Write();
}