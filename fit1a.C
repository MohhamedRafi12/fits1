#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TCanvas.h"

#include <iostream>
#include <vector>

struct results {
  float params[4];   // [0]=A, [1]=mean, [2]=sigma, [3]=reduced chi2
  float errs[3];     // errors on params 0..2
  double ndof;       // NDF
  double prob;       // chi2 probability (p-value)
};

void plot_results_basic(const std::vector<results>& res, const char* outpdf="result1.pdf") {
  gStyle->SetOptStat(1110);

  const int N = (int)res.size();
  if (N == 0) return;

  // Fixed/basic ranges
  TH1F* h_chi2   = new TH1F("h_chi2",   "Reduced #chi^{2} distribution;#chi^{2}/ndf;Counts", 50, 0.0, 3.0);
  TH1F* h_mean   = new TH1F("h_mean",   "Distribution of mean from fits;mean;Counts",         50, 40.0, 60.0);
  TH1F* h_prob   = new TH1F("h_prob",   "#chi^{2} Probability;P(#chi^{2});Counts",            50, 0.0, 1.0);
  TH1F* h_err_mu = new TH1F("h_err_mu", "Error on mean from fits;#sigma_{mean};Counts",       50, 0.0, 2.0);

  for (const auto& r : res) {
    h_chi2  ->Fill(r.params[3]); // reduced chi2
    h_mean  ->Fill(r.params[1]); // mean
    h_prob  ->Fill(r.prob);      // p-value
    h_err_mu->Fill(r.errs[1]);   // error on mean
  }

  TCanvas* c = new TCanvas("c", "Fit result distributions", 900, 800);
  c->Divide(2,2);
  c->cd(1); h_chi2  ->Draw();
  c->cd(2); h_mean  ->Draw();
  c->cd(3); h_prob  ->Draw();
  c->cd(4); h_err_mu->Draw();
  c->SaveAs(outpdf);
}

results fit1(int entries=1000, bool save=false) {
  gROOT->Reset(); // usually avoid in compiled code

  TFile *tf = nullptr;
  if (save) tf = new TFile("histo.root","recreate");

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 rng(0);

  for (int i=0 ; i<entries ; ++i)
    randomHist1->Fill(rng.Gaus(50,10)); // mean=50, sigma=10

  gStyle->SetOptFit(1111);
  randomHist1->Fit("gaus");
  TF1 *fitfunc = randomHist1->GetFunction("gaus");

  results res{};
  for (int i = 0; i < 3; ++i) {
    res.params[i] = fitfunc->GetParameter(i);
    res.errs[i]   = fitfunc->GetParError(i);
  }

  const double chi2 = fitfunc->GetChisquare();
  const double ndf  = fitfunc->GetNDF();
  res.params[3] = chi2/ndf; // reduced chi2
  res.ndof      = ndf;
  res.prob      = fitfunc->GetProb();

  if (save && tf) { tf->Write(); tf->Close(); }

  return res;
}

void fit1a() {
  std::vector<results> res;
  res.reserve(1e6);

  for (int i = 0; i < 1e6; ++i)
    res.push_back(fit1(1000));

  plot_results_basic(res, "result1.pdf");
}
