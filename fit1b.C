#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TCanvas.h"

#include <iostream>
#include <vector>
#include <string>

struct results {
  float params[4];   // [0]=A, [1]=mean, [2]=sigma, [3]=reduced chi2
  float errs[3];     // errors on params 0..2
  double ndof;       // NDF
  double prob;       // chi2 probability (p-value)
};

void plot_means_2x1(const std::vector<results>& res_chi2, const std::vector<results>& res_nll, const char* outpdf = "result2.pdf")
{
  if (res_chi2.empty() && res_nll.empty()) return;

  gStyle->SetOptStat(1110);

  TH1F h_mu_chi2("h_mu_chi2",
                 "Distribution of mean from #chi^{2} fits;#mu;Counts",
                 100, 0, 100);
  TH1F h_mu_nll ("h_mu_nll",
                 "Distribution of mean from NLL fits;#mu;Counts",
                 100, 0, 100);

  for (const auto& r : res_chi2) h_mu_chi2.Fill(r.params[1]); // mean
  for (const auto& r : res_nll)  h_mu_nll .Fill(r.params[1]); // mean

  h_mu_chi2.SetLineWidth(2);
  h_mu_nll .SetLineWidth(2);

  TCanvas c("c_result2", "result2", 900, 700);
  c.Divide(2,1);
  c.cd(1); h_mu_chi2.Draw("HIST");
  c.cd(2); h_mu_nll .Draw("HIST");
  c.SaveAs(outpdf);
}

results fit1(int entries=1000, bool nll=false, bool save=false) {
  gROOT->Reset(); // usually avoid in compiled code

  TFile *tf = nullptr;
  if (save) tf = new TFile("histo.root","recreate");

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 rng(0);

  for (int i=0 ; i<entries ; ++i)
    randomHist1->Fill(rng.Gaus(50,10)); // mean=50, sigma=10

  gStyle->SetOptFit(1111);
  
  if (nll) randomHist1->Fit("gaus","LQ");
  else randomHist1->Fit("gaus", "Q");

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

void fit1b() {
  std::vector<results> res_normal;
  res_normal.reserve(1000);
  std::vector<results> res_nll;
  res_nll.reserve(1000);

  for (int i = 0; i < 1000; ++i) {
    res_normal.push_back(fit1(10));
    res_nll.push_back(fit1(10, true));
  }
  plot_means_2x1(res_normal, res_nll);

}
