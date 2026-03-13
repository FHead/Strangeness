#include <algorithm>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {
constexpr double kFitMin = 1.000;
constexpr double kFitMax = 1.050;

TF1 buildTripleGauss(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
        "[5]*exp(-0.5*((x-[1])/[6])^2)",
        kFitMin, kFitMax);
  f.SetParNames("A1", "mean", "sigma1", "A2", "sigma2", "A3", "sigma3");
  f.SetParameters(std::max(50.0, 0.5 * h->GetMaximum()), 1.0195, 0.0018,
                  std::max(25.0, 0.3 * h->GetMaximum()), 0.0038,
                  std::max(10.0, 0.15 * h->GetMaximum()), 0.0065);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(4, 0.0003, 0.03);
  f.SetParLimits(6, 0.0003, 0.05);
  return f;
}

void styleHistogram(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.70);
  h->SetLineWidth(2);
}

void stylePull(TH1D& hPull) {
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.60);
  hPull.SetLineColor(kBlack);
  hPull.SetMarkerColor(kBlack);
  hPull.GetYaxis()->SetTitle("Pull");
  hPull.GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hPull.GetYaxis()->SetNdivisions(505);
  hPull.GetYaxis()->SetTitleSize(0.11);
  hPull.GetYaxis()->SetLabelSize(0.085);
  hPull.GetYaxis()->SetTitleOffset(0.55);
  hPull.GetXaxis()->SetTitleSize(0.11);
  hPull.GetXaxis()->SetLabelSize(0.085);
  hPull.GetXaxis()->SetTitleOffset(1.05);
  hPull.SetMinimum(-5.0);
  hPull.SetMaximum(5.0);
}

void drawCategoryPanel(TPad* padTop, TPad* padBottom, TH1D* h, const std::string& title, int color) {
  padTop->cd();
  styleHistogram(h, color);
  h->SetTitle("");
  h->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  h->GetYaxis()->SetTitle("Candidates / bin");

  TF1 fit = buildTripleGauss(std::string("f_") + h->GetName(), h);
  h->Fit(&fit, "RQ0");
  const double originalWidth = fit.GetXmax() > fit.GetXmin() ? (kFitMax - kFitMin) / 200.0 : 1.0;
  const double displayWidth = h->GetXaxis()->GetBinWidth(1);
  const double scale = displayWidth / ((1.06 - 0.99) / 280.0);
  fit.SetParameter(0, fit.GetParameter(0) * scale);
  fit.SetParameter(3, fit.GetParameter(3) * scale);
  fit.SetParameter(5, fit.GetParameter(5) * scale);
  fit.SetLineColor(kOrange + 7);
  fit.SetLineWidth(4);

  h->Draw("E");
  fit.Draw("same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.16, 0.88, title.c_str());
  latex.SetTextSize(0.036);
  latex.DrawLatex(0.16, 0.81, Form("Triple-Gauss fit, #chi^{2}/ndf = %.2f", fit.GetChisquare() / fit.GetNDF()));
  latex.DrawLatex(0.16, 0.74, Form("mean = %.5f GeV", fit.GetParameter(1)));

  TLegend leg(0.64, 0.76, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h, "MC histogram", "lep");
  leg.AddEntry(&fit, "Best-fit curve", "l");
  leg.Draw();

  padBottom->cd();
  TH1D hPull((std::string("hPull_") + h->GetName()).c_str(), "",
             h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  stylePull(hPull);

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < kFitMin || x > kFitMax) continue;
    const double y = h->GetBinContent(i);
    const double ey = h->GetBinError(i);
    if (ey <= 0.0) continue;
    hPull.SetBinContent(i, (y - fit.Eval(x)) / ey);
    hPull.SetBinError(i, 1.0);
  }

  hPull.Draw("E");
  TLine line0(kFitMin, 0.0, kFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");
}
}

int PlotPhiSignalOnlyBestFits(std::string input = "PhiSignalOnlyHistograms.root",
                              std::string outputDir = "Plots") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open " << input << std::endl;
    return 1;
  }

  TH1D* h1 = nullptr;
  TH1D* h2 = nullptr;
  file.GetObject("hPhiMass1Tag", h1);
  file.GetObject("hPhiMass2Tag", h2);
  if (h1 == nullptr || h2 == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  TH1D* h1Draw = static_cast<TH1D*>(h1->Clone("h1DrawBest"));
  TH1D* h2Draw = static_cast<TH1D*>(h2->Clone("h2DrawBest"));
  h1Draw->Rebin(4);
  h2Draw->Rebin(4);

  TCanvas c("c", "c", 850, 980);
  TPad p1Top("p1Top", "p1Top", 0.0, 0.56, 1.0, 1.0);
  TPad p1Bottom("p1Bottom", "p1Bottom", 0.0, 0.39, 1.0, 0.56);
  TPad p2Top("p2Top", "p2Top", 0.0, 0.17, 1.0, 0.39);
  TPad p2Bottom("p2Bottom", "p2Bottom", 0.0, 0.00, 1.0, 0.17);

  for (TPad* pad : {&p1Top, &p1Bottom, &p2Top, &p2Bottom}) {
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.05);
  }
  p1Top.SetBottomMargin(0.02);
  p1Bottom.SetTopMargin(0.03);
  p1Bottom.SetBottomMargin(0.22);
  p2Top.SetTopMargin(0.05);
  p2Top.SetBottomMargin(0.02);
  p2Bottom.SetTopMargin(0.03);
  p2Bottom.SetBottomMargin(0.22);

  p1Top.Draw();
  p1Bottom.Draw();
  p2Top.Draw();
  p2Bottom.Draw();

  drawCategoryPanel(&p1Top, &p1Bottom, h1Draw, "1-tag category", kBlue + 1);
  drawCategoryPanel(&p2Top, &p2Bottom, h2Draw, "2-tag category", kRed + 1);

  c.SaveAs((outputDir + "/phi_signal_only_best_fits.pdf").c_str());

  delete h1Draw;
  delete h2Draw;
  file.Close();
  return 0;
}
