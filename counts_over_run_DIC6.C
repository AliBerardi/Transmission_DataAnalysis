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
const double tofMin_A1 = 1.e5;
const double tofMax_A1 = 1.e8;
/*
void count_single_run()
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
  graph_A1->SetTitle(";Run number; DIC6 counts / N protons from PKUP");
  graph_A1->SetMarkerStyle(20);

  // loop over all the runs of DIC6
  for (int k = 0; k < n_runs; ++k)
  { /// CAMBIA QUI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        if (tflash >= 11070 && tflash <= 11130 && tof - tflash > 0)
        {

          hAmpTof->Fill(tof - tflash, amp);
        }

      } // close loop over bunches

    } // close loop on tree entries

    cout << "Processed entries: " << counter << '\n';

    const int xFirst = hAmpTof->GetXaxis()->FindBin(tofMin_A1);
    const int xLast = hAmpTof->GetXaxis()->FindBin(tofMax_A1) - 1;
    const int yFirst = hAmpTof->GetYaxis()->FindBin(ampMin);
    const int yLast = hAmpTof->GetYaxis()->FindBin(ampMax) - 1;

    counts_DIC_A1 = static_cast<Long64_t>(
        hAmpTof->Integral(xFirst, xLast, yFirst, yLast));

    std::cout << "Counts of DIC6 detetcor in range tof:[" << tofMin_A1 << ", " << tofMax_A1 << "] and amp:[" << ampMin << ", " << ampMax << "] for current run: " << counts_DIC_A1 << std::endl;

    if (PI_pkup > 0)
    {
      double ratio = static_cast<double>(counts_DIC_A1) / PI_pkup;
      graph_A1->SetPoint(graph_A1->GetN(), k + 1,
                         ratio);
      std::cout << "DIC6 counts/N protons PKUP for A1: " << static_cast<double>(counts_DIC_A1) / PI_pkup << std::endl;
      std::cout << "" << std::endl;
      graph_A1->SetPointError(graph_A1->GetN() - 1, 0, ratio * sqrt(1. / counts_DIC_A1 + 1. / PI_pkup));
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

  TCanvas *cRatio = new TCanvas("cRatio", "DIC6 counts per PI_PKUP", 1200, 800);
  graph_A1->Draw("AP");
  TLegend *cutsLegend = new TLegend(0.62, 0.76, 0.89, 0.88);
  cutsLegend->SetBorderSize(0);
  cutsLegend->SetFillStyle(0);
  cutsLegend->SetTextSize(0.045);
  cutsLegend->AddEntry(static_cast<const TObject *>(nullptr),
                       Form("Amplitude cut: %.0f - %.0f channels", ampMin, ampMax),
                       "");
  cutsLegend->AddEntry(static_cast<const TObject *>(nullptr),
                       Form("TOF cut: %.1e - %.1e ns", tofMin_A1, tofMax_A1),
                       "");
  cutsLegend->Draw();
  // cRatio->SaveAs("OUTPUT/DIC6/counts_over_run_DIC6_e5_e8_ampcut_a.root");
}
*/

void count_pairs(int runsPerGroup)
{

  std::cout << "\n Processing runs in groups of " << runsPerGroup << "\n" << std::endl;


  if (runsPerGroup <= 0)
  {
    std::cerr << "ERROR: runsPerGroup must be greater than zero."
              << std::endl;
    return;
  }

  string cmnd_filename =
      "input_files/Transmission_ratio_final.cmnd";
  ConfigReader cfg(cmnd_filename);

  string prefix = cfg.getString("prefix");
  string suffix = cfg.getString("suffix");

  vector<int> SIN = cfg.getIntVector("SIN");
  vector<int> SOUT = cfg.getIntVector("SOUT");

  vector<int> ALL_DATA;
  ALL_DATA.insert(ALL_DATA.end(), SIN.begin(), SIN.end());
  ALL_DATA.insert(ALL_DATA.end(), SOUT.begin(), SOUT.end());
  std::sort(ALL_DATA.begin(), ALL_DATA.end());

  double tof;
  double tflash;
  float amp;
  int BunchNumber;

  int BunchNumberPK;
  float PulseIntensity, PulseIntensityPK;
  double tflash_p;
  int PSpulse, PSpulsePK;

  TTree *t_DIC6, *t_PKUP;

  const int n_runs = ALL_DATA.size();
  int Run_i = 0;
  Char_t path_ri[100];

  int BinPerDecade = 500;
  int N_DEC_tof = 8;
  int n_vec_tof = BinPerDecade * N_DEC_tof;
  double xbins_tof[100100];
  double step_tof = static_cast<double>(N_DEC_tof) / n_vec_tof;
  for (int filler = 0; filler <= n_vec_tof; ++filler)
  {
    xbins_tof[filler] =
        pow(10., step_tof * static_cast<double>(filler));
  }

  TH2D *hAmpTof =
      new TH2D("hAmpTof_groups", "hAmpTof_groups", n_vec_tof, xbins_tof,
               500, 0, 5.e+3);
  hAmpTof->Sumw2();

  TGraphErrors *graph_A1 = new TGraphErrors();
  graph_A1->SetName("graph_A1_groups");
  graph_A1->SetTitle(Form("Group of %d runs;DIC6 counts / N protons from PKUP", runsPerGroup));
  graph_A1->GetYaxis()->SetTitle("DIC6 counts / N protons from PKUP");
  graph_A1->GetXaxis()->SetTitle(Form("Group of %d runs", runsPerGroup));
  graph_A1->SetMarkerStyle(20);

  for (int k = 0; k < n_runs; k += runsPerGroup)
  {//ciclo sulle run
    const int groupEnd =
        (k + runsPerGroup < n_runs) ? k + runsPerGroup : n_runs;
    double PI_pkup = 0;
    Long64_t counts_DIC_A1 = 0;

    hAmpTof->Reset();

    for (int runIndex = k; runIndex < groupEnd; ++runIndex)
    {//ciclo sulle run del gruppo
      Run_i = ALL_DATA[runIndex];
      sprintf(path_ri, "%s%i%s", prefix.c_str(), Run_i, suffix.c_str());
      std::cout << "Now Sorting " << path_ri << std::endl;
      double PI_run = 0;

      TFile *f_i = TFile::Open(path_ri, "READ");
      if (!f_i || f_i->IsZombie() || !f_i->IsOpen())
      {
        std::cerr << "ERROR: cannot open file: " << path_ri << "\n";
        if (f_i)
        {
          delete f_i;
        }
        continue;
      }

      t_PKUP = (TTree *)f_i->Get("PKUP");
      t_PKUP->SetBranchStatus("*", 0);
      t_PKUP->SetBranchStatus("tflash", 1);
      t_PKUP->SetBranchStatus("BunchNumber", 1);
      t_PKUP->SetBranchStatus("PulseIntensity", 1);
      t_PKUP->SetBranchStatus("PSpulse", 1);

      t_PKUP->SetBranchAddress("tflash", &tflash_p);
      t_PKUP->SetBranchAddress("BunchNumber", &BunchNumberPK);
      t_PKUP->SetBranchAddress("PulseIntensity", &PulseIntensityPK);
      t_PKUP->SetBranchAddress("PSpulse", &PSpulsePK);

      Long64_t nentryPK = t_PKUP->GetEntriesFast();
      if (nentryPK < 0)
        nentryPK = t_PKUP->GetEntries();

      t_DIC6 = (TTree *)f_i->Get("DIC6");
      t_DIC6->SetBranchStatus("*", 0);
      t_DIC6->SetBranchStatus("tof", 1);
      t_DIC6->SetBranchStatus("tflash", 1);
      t_DIC6->SetBranchStatus("BunchNumber", 1);
      t_DIC6->SetBranchStatus("amp", 1);
      t_DIC6->SetBranchStatus("PulseIntensity", 1);
      t_DIC6->SetBranchStatus("PSpulse", 1);

      t_DIC6->SetBranchAddress("tof", &tof);
      t_DIC6->SetBranchAddress("tflash", &tflash);
      t_DIC6->SetBranchAddress("amp", &amp);
      t_DIC6->SetBranchAddress("BunchNumber", &BunchNumber);
      t_DIC6->SetBranchAddress("PulseIntensity", &PulseIntensity);
      t_DIC6->SetBranchAddress("PSpulse", &PSpulse);

      Long64_t nentry_DIC6 = t_DIC6->GetEntriesFast();
      if (nentry_DIC6 < 0)
        nentry_DIC6 = t_DIC6->GetEntries();

      std::cout << "Number of Entries in PKUP: " << nentryPK << std::endl;
      std::cout << "Number of Entries in DIC6: " << nentry_DIC6 << std::endl;

      t_DIC6->AddFriend(t_PKUP);
      t_PKUP->BuildIndex("BunchNumber");

      for (Long64_t iP = 0; iP < nentryPK; ++iP)
      {
        t_PKUP->GetEntry(iP);
        PI_run += PulseIntensityPK;
      }
      PI_pkup += PI_run;

      std::cout << "Total number of protons in run 121" << Run_i << ": "
                << PI_run << std::endl;

      int counter = 0;
      for (Long64_t i = 0; i < nentry_DIC6; ++i)
      {
        t_DIC6->GetEntry(i);

        if (BunchNumberPK == BunchNumber)
        {
          ++counter;
          if (tflash >= 11070 && tflash <= 11130 && tof - tflash > 0)
          {
            hAmpTof->Fill(tof - tflash, amp);
          }
        }
      }

      std::cout << "Processed entries for run 121" << Run_i << ": "
                << counter << '\n';

      f_i->Close();
      delete f_i;
    }//chiude ciclo sulle run del gruppo

    const int xFirst = hAmpTof->GetXaxis()->FindBin(tofMin_A1);
    const int xLast = hAmpTof->GetXaxis()->FindBin(tofMax_A1) - 1;
    const int yFirst = hAmpTof->GetYaxis()->FindBin(ampMin);
    const int yLast = hAmpTof->GetYaxis()->FindBin(ampMax) - 1;

    counts_DIC_A1 = static_cast<Long64_t>(
        hAmpTof->Integral(xFirst, xLast, yFirst, yLast));

    std::cout << "Counts of DIC6 detector in range tof:[" << tofMin_A1
              << ", " << tofMax_A1 << "] and amp:[" << ampMin << ", "
              << ampMax << "] for run group [" << ALL_DATA[k] << ", "
              << ALL_DATA[groupEnd - 1] << "]: " << counts_DIC_A1
              << std::endl;

    if (PI_pkup > 0)
    {
      const double ratio = static_cast<double>(counts_DIC_A1) / PI_pkup;
      graph_A1->SetPoint(graph_A1->GetN(), k / runsPerGroup + 1, ratio);
      graph_A1->SetPointError(
          graph_A1->GetN() - 1, 0,
          ratio * sqrt(1. / counts_DIC_A1 + 1. / PI_pkup));
      std::cout << "DIC6 counts/N protons PKUP for run group ["
                << ALL_DATA[k] << ", " << ALL_DATA[groupEnd - 1]
                << "]: " << ratio << std::endl;
                std::cout << "" << std::endl;
    }
    else
    {
      std::cerr << "WARNING: PI_PKUP is zero for run group ["
                << ALL_DATA[k] << ", " << ALL_DATA[groupEnd - 1]
                << "]; skipping graph point." << std::endl;
    }
  }
  std::cout << "" << std::endl;

  
  /*
  // derive histogram range and binning from the ratio values
  double ratioMin =  std::numeric_limits<double>::infinity();
  double ratioMax = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < graph_A1->GetN(); ++i)
  {
    double x, ratio;
    graph_A1->GetPoint(i, x, ratio);
    ratioMin = std::min(ratioMin, ratio);
    ratioMax = std::max(ratioMax, ratio);
  }

  if (graph_A1->GetN() == 0)
  {
    std::cerr << "ERROR: no valid ratio values available; "
                 "cannot create ratio histogram."
              << std::endl;
    return;
  }

  double ratioSpan = ratioMax - ratioMin;
  if (ratioSpan == 0.)
  {
    ratioSpan = std::max(std::abs(ratioMin) * 0.10,
                         1.e-15);
  }
  const double ratioMargin = 0.10 * ratioSpan;
  ratioMin -= ratioMargin;
  ratioMax += ratioMargin;

  const int nRatioBins =
      std::max(10, std::min(100,
                            static_cast<int>(std::ceil(
                                std::sqrt(graph_A1->GetN())))));
                              */
  const int nRatioBins = 20;
  const double ratioMin = 100.e-15;
  const double ratioMax = 110.e-15;
  TH1D *hRatio =
      new TH1D("hRatio",
               Form("Ratio values for groups of %d runs", runsPerGroup),
               nRatioBins, ratioMin, ratioMax);
  hRatio->GetXaxis()->SetTitle("DIC6 counts / N protons from PKUP");
  hRatio->GetYaxis()->SetTitle("Counts");

  
  for (int i = 0; i < graph_A1->GetN(); ++i)
  {
    double x, ratio;
    graph_A1->GetPoint(i, x, ratio);
    hRatio->Fill(ratio);
  }

  TFitResultPtr Results = hRatio->Fit("gaus", "S", "");//, 190.e-15, 400.e-15);
  double mean = Results->Parameter(1);
  double sigma = Results->Parameter(2);


  //canvas per disegnare il grafico ratio over run
  TCanvas *cRatio =
      new TCanvas("cRatio_groups", "DIC6 counts per PI_PKUP", 1200, 800);
  graph_A1->Draw("AP");
  graph_A1->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  TLegend *cutsLegend = new TLegend(0.62, 0.76, 0.89, 0.88);
  cutsLegend->SetBorderSize(0);
  cutsLegend->SetFillStyle(0);
  cutsLegend->SetTextSize(0.045);
  cutsLegend->AddEntry(
      static_cast<const TObject *>(nullptr),
      Form("Amplitude cut: %.0f - %.0f channels", ampMin, ampMax), "");
  cutsLegend->AddEntry(
      static_cast<const TObject *>(nullptr),
      Form("TOF cut: %.1e - %.1e ns", tofMin_A1, tofMax_A1), "");
  cutsLegend->Draw();
  //disegno linea in corrispondenza del valore medio estratto dal fit gaussiano
  TLine *meanLine = new TLine(0, mean, graph_A1->GetN() + 1, mean);
  meanLine->SetLineColor(kRed);
  meanLine->SetLineWidth(2);
  meanLine->Draw();
  //disegno linee in corrispondenza del valore medio +- sigma
  TLine *sigma_up = new TLine(0, mean + (2*sigma), graph_A1->GetN() + 1, mean + (2*sigma));
  sigma_up->SetLineColor(kRed);
  sigma_up->SetLineWidth(2);
  sigma_up->SetLineStyle(2); //dashed line
  sigma_up->Draw();
  TLine *sigma_down = new TLine(0, mean - (2*sigma), graph_A1->GetN() + 1, mean - (2*sigma));
  sigma_down->SetLineColor(kRed);
  sigma_down->SetLineWidth(2);
  sigma_down->SetLineStyle(2); //dashed line
  sigma_down->Draw();

  

  TFile *fr = TFile::Open(Form("OUTPUT/DIC6/counts_over_run_DIC6_group%druns_e5e8_solotrizio.root", runsPerGroup), "RECREATE");
  cRatio->Write();


  //disegno l'istogramma che raggruppa i valori di ratio e estrae media e sigma
  TCanvas *cFit = new TCanvas("cFit", "Gaussian fit of Ratio values", 1200, 800);
  hRatio->Draw();
  TLegend *fitLegend = new TLegend(0.62, 0.76, 0.89, 0.88);
  fitLegend->SetBorderSize(0);
  fitLegend->SetFillStyle(0);
  fitLegend->SetTextSize(0.045);
fitLegend->AddEntry(static_cast<const TObject *>(nullptr),
                    "Gaussian fit:", "");
fitLegend->AddEntry(static_cast<const TObject *>(nullptr),
                    Form("mean = %.2e", mean), "");
fitLegend->AddEntry(static_cast<const TObject *>(nullptr),
                    Form("sigma = %.2e", sigma), "");
fitLegend->AddEntry(static_cast<const TObject *>(nullptr),
                    Form("sigma/mean = %.2e", sigma / mean), "");
fitLegend->Draw();

  cFit->Write();
}

void counts_over_run_DIC6(int RunsPerGroup)
{

  // count_single_run();
  count_pairs(RunsPerGroup);
}