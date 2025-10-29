#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMath.h"
#include "TFile.h"
#include "TKey.h"
#include "TIterator.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <algorithm>

// =====================================================
// --- Poisson negative log-likelihood for fixed mean ---
// =====================================================
static double NLL_for_mean(const TH1* h, TF1& f, double mean_val) {
  f.SetParameter(1, mean_val);
  const int nb = h->GetNbinsX();
  double nll = 0.0;

  for (int i = 1; i <= nb; ++i) {
    const double n   = h->GetBinContent(i);
    const double xl  = h->GetBinLowEdge(i);
    const double xr  = xl + h->GetBinWidth(i);
    double mu = f.Integral(xl, xr);
    if (mu <= 0) mu = 1e-12;
    nll += (mu - n * std::log(mu));
  }
  return nll;
}

// =====================================================
// --- Chi2 for fixed mean ---
// =====================================================
static double Chi2_for_mean(const TH1* h, TF1& f, double mean_val) {
  f.SetParameter(1, mean_val);
  const int nb = h->GetNbinsX();
  double chi2 = 0.0;

  for (int i = 1; i <= nb; ++i) {
    const double n   = h->GetBinContent(i);
    const double xl  = h->GetBinLowEdge(i);
    const double xr  = xl + h->GetBinWidth(i);
    double mu = f.Integral(xl, xr);
    if (mu <= 0) mu = 1e-12;
    const double diff = n - mu;
    chi2 += diff * diff / mu;
  }
  return chi2;
}

// =====================================================
// --- Get first histogram from file ---
// =====================================================
static TH1* GetFirstTH1(TFile* f) {
  if (!f) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj->InheritsFrom(TH1::Class())) return (TH1*)obj;
  }
  return nullptr;
}

// =====================================================
// --- Main plotting function ---
// =====================================================
void histfit2(const char* infile="histo1k.root") {
  // Load histogram
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open " << infile << "\n";
    return;
  }
  TH1* h = GetFirstTH1(fin);
  if (!h) { std::cerr << "No histogram found\n"; return; }

  // Gaussian model
  TF1 f("gaus","gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  f.SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  const double mu0 = f.GetParameter(1);
  const double s0  = f.GetParameter(2);

  // =====================================================
  // --- 1. NLL scan around minimum pm4 ---
  // =====================================================
  const int NPTS1 = 200;
  const double muMin1 = mu0 - 2.0 * s0; // wide enough around min
  const double muMax1 = mu0 + 2.0 * s0;
  std::vector<double> xs1(NPTS1), ys1(NPTS1);

  for (int i = 0; i < NPTS1; ++i) {
    const double mu = muMin1 + (muMax1 - muMin1) * (i / double(NPTS1 - 1));
    xs1[i] = mu;
    ys1[i] = NLL_for_mean(h, f, mu);
  }

  const double nllMin = *std::min_element(ys1.begin(), ys1.end());
  for (auto& y : ys1) y = 2.0 * (y - nllMin); // -2lnL form

  TGraph gNLL(NPTS1, xs1.data(), ys1.data());
  gNLL.SetTitle("-2#Delta ln L vs Mean;Mean; -2#Delta ln L");
  gNLL.SetLineWidth(2);
  gNLL.SetLineColor(kBlue+1);

  // Draw NLL plot
  TCanvas c1("c1","NLL vs Mean",900,700);
  gNLL.Draw("AL");

  // Limit y-range to pm4 region
  gNLL.GetYaxis()->SetRangeUser(0, 4.5);
  TLine l1(mu0, 1, mu0 + s0, 1); // horizontal guides
  TLine l4(mu0, 4, mu0 + s0, 4);
  l1.SetLineStyle(2); l4.SetLineStyle(2);
  l1.Draw("same"); l4.Draw("same");

  c1.SaveAs("results4_1k_nll.pdf");

  // =====================================================
  // --- 2. Chi2 scan centered at 53, range pm1 ---
  // =====================================================
  const double center = 50;
  const int NPTS2 = 200;
  const double muMin2 = center - 3.0;
  const double muMax2 = center + 3.0;

  std::vector<double> xs2(NPTS2), ys2(NPTS2);
  for (int i = 0; i < NPTS2; ++i) {
    const double mu = muMin2 + (muMax2 - muMin2) * (i / double(NPTS2 - 1));
    xs2[i] = mu;
    ys2[i] = Chi2_for_mean(h, f, mu);
  }

  TGraph gChi2(NPTS2, xs2.data(), ys2.data());
  gChi2.SetTitle("#chi^{2} vs Mean;Mean;#chi^{2}");
  gChi2.SetLineWidth(2);
  gChi2.SetLineColor(kRed+1);

  // Draw Chi2 plot
  TCanvas c2("c2","Chi2 vs Mean",900,700);
  gChi2.Draw("AL");

  // Vertical guide at 50.8
  TLine lcenter(center, gChi2.GetYaxis()->GetXmin(),
                center, gChi2.GetYaxis()->GetXmax());
  lcenter.SetLineColor(kBlack);
  lcenter.SetLineStyle(2);
  lcenter.Draw();

  c2.SaveAs("results4_1k_chi2.pdf");

  fin->Close();
}
