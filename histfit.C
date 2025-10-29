#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include <algorithm>
#include "TKey.h"
#include "TIterator.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include <iostream>
#include <vector>

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

// Poisson negative log-likelihood (up to additive constants):
// NLL = sum_i [ mu_i - n_i * ln(mu_i) ]
// Use bin-integrals for mu_i for better asymptotic behavior.
static double NLL_from_hist_and_model(const TH1* h, TF1* f) {
  const int nb = h->GetNbinsX();
  double nll = 0.0;
  for (int i = 1; i <= nb; ++i) {
    const double n   = h->GetBinContent(i);
    const double bw  = h->GetBinWidth(i);
    const double xl  = h->GetBinLowEdge(i);
    const double xr  = xl + bw;
    double mu = f->Integral(xl, xr);
    if (mu <= 0) mu = 1e-12; // guard
    nll += (mu - n * std::log(mu));
  }
  return nll;
}

void histFit(const char* infile="histo25.root",
             const char* outpdf="result3.pdf",
             int ntoy=10000)
{
  gStyle->SetOptStat(1110);

  // Load histogram
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open file " << infile << "\n";
    return;
  }
  TH1* hdata_in = GetFirstTH1(fin);
  if (!hdata_in) {
    std::cerr << "ERROR: no TH1 found in " << infile << "\n";
    return;
  }

  TH1* hdata = (TH1*)hdata_in->Clone("hdata_detached");
  hdata->SetDirectory(nullptr);
  fin->Close();

 
  if (hdata->GetListOfFunctions()) hdata->GetListOfFunctions()->Clear();

  TF1 gaus("gaus","gaus", hdata->GetXaxis()->GetXmin(), hdata->GetXaxis()->GetXmax());
  TFitResultPtr r = hdata->Fit(&gaus, "LQ0NS"); // store covariance
  const double mu_err = r ? r->ParError(1) : gaus.GetParError(1);

  const double nll_data = NLL_from_hist_and_model(hdata, &gaus);

  TRandom3 rng(0);
  const int nb = hdata->GetNbinsX();

  std::vector<double> mu(nb+1, 0.0); // 1..nb
  for (int i = 1; i <= nb; ++i) {
    const double xl = hdata->GetXaxis()->GetBinLowEdge(i);
    const double xr = hdata->GetXaxis()->GetBinUpEdge(i);
    mu[i] = gaus.Integral(xl, xr);
    if (mu[i] <= 0) mu[i] = 1e-12;
  }

  const double span = std::max(10.0, 3.0 * std::sqrt(2.0 * nb));
  TH1F hNLL("hNLL", "NLL distribution from pseudo-experiments;NLL;Counts", 80,
            nll_data - span, nll_data + span);

  int geq = 0;
  TH1D htoy("htoy","toy", nb, hdata->GetXaxis()->GetXmin(), hdata->GetXaxis()->GetXmax());

  for (int t = 0; t < ntoy; ++t) {
    for (int i = 1; i <= nb; ++i) {
      const double n = rng.Poisson(mu[i]);
      htoy.SetBinContent(i, n);
    }
    const double nll_toy = NLL_from_hist_and_model(&htoy, &gaus);
    hNLL.Fill(nll_toy);
    if (nll_toy >= nll_data) ++geq;
  }

  const double pval = (double)geq / (double)ntoy;

  // Draw and annotate
  TCanvas c("c","result3", 900, 700);
  hNLL.SetLineWidth(2);
  hNLL.Draw("HIST");

  TLine line(nll_data, 0, nll_data, hNLL.GetMaximum()*1.05);
  line.SetLineColor(kRed+1);
  line.SetLineWidth(3);
  line.Draw();

  TLegend leg(0.55, 0.78, 0.88, 0.90);
  leg.AddEntry(&hNLL, "Toys NLL", "l");
  leg.AddEntry(&line, Form("Data NLL = %.2f", nll_data), "l");
  leg.Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.55, 0.72, Form("p-value (NLL_{toy} #geq NLL_{data}) = %.3f", pval));

  c.SaveAs(outpdf);

  std::cout << "Data NLL = " << nll_data << "\n";
  std::cout << "p-value  = " << pval << "  (fraction of toys with NLL >= data)\n";
}
