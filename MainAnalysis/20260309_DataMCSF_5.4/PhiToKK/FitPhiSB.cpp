#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kThreshold = 2.0 * kKaonMass;
constexpr double kFitMin = 0.99;
constexpr double kFitMax = 1.06;
constexpr int kDisplayRebin = 4;

struct SignalShape {
  double mean = 1.0195;
  double sigmaG = 0.0023;
  double ratio = 0.2;
  double sigmaCB = 0.0060;
  double alpha = 1.5;
  double n = 8.0;
};

struct FitSummary {
  std::string category;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 0.0;
  double signalAmp = 0.0;
  double bgGaussAmp1 = 0.0;
  double bgGaussMean1 = 0.0;
  double bgGaussSigma1 = 0.0;
  double bgGaussAmp2 = 0.0;
  double bgGaussMean2 = 0.0;
  double bgGaussSigma2 = 0.0;
  double thresholdNorm = 0.0;
  double thresholdPower = 0.0;
  double thresholdSlope = 0.0;
};

SignalShape gSignalShape;

double RightTailCBUnit(double x, double mean, double sigma, double alpha, double n) {
  if (sigma <= 0.0 || alpha <= 0.0 || n <= 1.0) return 0.0;
  const double t = (x - mean) / sigma;
  if (t <= alpha) return std::exp(-0.5 * t * t);
  const double A = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double B = n / alpha - alpha;
  return A / std::pow(B + t, n);
}

double SignalOnlyShape(double* x, double* p) {
  const double xx = x[0];
  const double amp = p[0];
  const double gaussian = std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.sigmaG, 2));
  const double cb = RightTailCBUnit(xx, gSignalShape.mean, gSignalShape.sigmaCB, gSignalShape.alpha, gSignalShape.n);
  return amp * (gaussian + gSignalShape.ratio * cb);
}

double BackgroundShape(double* x, double* p) {
  const double xx = x[0];
  const double gauss1 = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  const double gauss2 = p[3] * std::exp(-0.5 * std::pow((xx - p[4]) / p[5], 2));
  double threshold = 0.0;
  if (xx > kThreshold)
    threshold = p[6] * std::pow(xx - kThreshold, p[7]) * std::exp(p[8] * xx);
  return gauss1 + gauss2 + threshold;
}

double TotalShape(double* x, double* p) {
  double signalP[1] = {p[0]};
  double backgroundP[9] = {p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]};
  return SignalOnlyShape(x, signalP) + BackgroundShape(x, backgroundP);
}

SignalShape deriveSignalShape(TH1D* h) {
  TF1 f("fSignalShape", [](double* x, double* p) {
    const double xx = x[0];
    const double gNorm = p[0];
    const double mean = p[1];
    const double sigmaG = p[2];
    const double cbNorm = p[3];
    const double sigmaCB = p[4];
    const double alpha = p[5];
    const double n = p[6];
    const double tg = (xx - mean) / sigmaG;
    const double gauss = gNorm * std::exp(-0.5 * tg * tg);
    double cb = 0.0;
    if (sigmaCB > 0.0 && alpha > 0.0 && n > 1.0) {
      const double tc = (xx - mean) / sigmaCB;
      if (tc <= alpha)
        cb = cbNorm * std::exp(-0.5 * tc * tc);
      else {
        const double A = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
        const double B = n / alpha - alpha;
        cb = cbNorm * A / std::pow(B + tc, n);
      }
    }
    return gauss + cb;
  }, 1.000, 1.050, 7);

  f.SetParameters(std::max(60.0, 0.6 * h->GetMaximum()), 1.0195, 0.0016,
                  std::max(30.0, 0.3 * h->GetMaximum()), 0.0035, 1.6, 8.0);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(4, 0.0003, 0.03);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  h->Fit(&f, "RQ0");

  SignalShape shape;
  shape.mean = f.GetParameter(1);
  shape.sigmaG = f.GetParameter(2);
  shape.ratio = (f.GetParameter(0) != 0.0) ? f.GetParameter(3) / f.GetParameter(0) : 0.0;
  shape.sigmaCB = f.GetParameter(4);
  shape.alpha = f.GetParameter(5);
  shape.n = f.GetParameter(6);
  return shape;
}

FitSummary fitCategory(TH1D* hSignal, TH1D* hSB, const std::string& category,
                       const std::string& outputDir) {
  gSignalShape = deriveSignalShape(hSignal);

  TF1 total(("fTotal_" + category).c_str(), TotalShape, kFitMin, kFitMax, 10);
  total.SetParNames("S", "G1amp", "G1mean", "G1sigma", "G2amp", "G2mean", "G2sigma", "N", "p", "b1");
  total.SetParameters(std::max(100.0, 0.4 * hSB->GetMaximum()),
                      std::max(100.0, 0.08 * hSB->GetMaximum()), 1.000, 0.003,
                      std::max(60.0, 0.05 * hSB->GetMaximum()), 1.004, 0.006,
                      std::max(100.0, 0.2 * hSB->GetMaximum()), 0.8, -2.0);
  total.SetParLimits(0, 0.0, 1e7);
  total.SetParLimits(1, 0.0, 1e7);
  total.SetParLimits(2, 0.995, 1.010);
  total.SetParLimits(3, 0.0005, 0.02);
  total.SetParLimits(4, 0.0, 1e7);
  total.SetParLimits(5, 0.995, 1.015);
  total.SetParLimits(6, 0.0005, 0.03);
  total.SetParLimits(7, 0.0, 1e9);
  total.SetParLimits(8, 0.0, 10.0);
  total.SetParLimits(9, -100.0, 20.0);
  hSB->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hSB->Clone((std::string(hSB->GetName()) + "_disp").c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double displayScale = hDisp->GetXaxis()->GetBinWidth(1) / hSB->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw(("fSignal_" + category).c_str(), SignalOnlyShape, kFitMin, kFitMax, 1);
  signalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 backgroundDraw(("fBackground_" + category).c_str(), BackgroundShape, kFitMin, kFitMax, 9);
  backgroundDraw.SetParameters(total.GetParameter(1) * displayScale, total.GetParameter(2), total.GetParameter(3),
                               total.GetParameter(4) * displayScale, total.GetParameter(5), total.GetParameter(6),
                               total.GetParameter(7) * displayScale, total.GetParameter(8), total.GetParameter(9));
  backgroundDraw.SetLineColor(kGreen + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(2);

  TF1 totalDraw(total);
  totalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  totalDraw.SetParameter(1, total.GetParameter(1) * displayScale);
  totalDraw.SetParameter(4, total.GetParameter(4) * displayScale);
  totalDraw.SetParameter(7, total.GetParameter(7) * displayScale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TCanvas c(("c_" + category).c_str(), ("c_" + category).c_str(), 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle((category + " reco-only S+B fit").c_str());
  hDisp->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Pairs / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  signalDraw.Draw("same");
  backgroundDraw.Draw("same");

  TLegend leg(0.56, 0.64, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event pairs", "lep");
  leg.AddEntry(&totalDraw, "Signal + background", "l");
  leg.AddEntry(&signalDraw, "Fixed signal shape", "l");
  leg.AddEntry(&backgroundDraw, "Background component", "l");
  leg.Draw();

  TPaveText txt(0.12, 0.60, 0.52, 0.88, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("Signal mean fixed at %.5f GeV", gSignalShape.mean));
  txt.AddText(Form("Bkg G1 mean = %.5f GeV, G2 mean = %.5f GeV", total.GetParameter(2), total.GetParameter(5)));
  txt.AddText(Form("Threshold p = %.3f, b_{1} = %.3f", total.GetParameter(8), total.GetParameter(9)));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull(("hPull_" + category).c_str(), "", hDisp->GetNbinsX(),
             hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
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

  for (int i = 1; i <= hDisp->GetNbinsX(); ++i) {
    const double xMin = hDisp->GetBinLowEdge(i);
    const double xMax = xMin + hDisp->GetBinWidth(i);
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0.0) continue;
    const double yfit = total.Integral(xMin, xMax) / hSB->GetXaxis()->GetBinWidth(1);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(kFitMin, 0.0, kFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/phi_sb_" + category + "_fit.pdf").c_str());
  delete hDisp;

  FitSummary summary;
  summary.category = category;
  summary.chi2 = total.GetChisquare();
  summary.ndf = total.GetNDF();
  summary.chi2ndf = total.GetChisquare() / total.GetNDF();
  summary.signalAmp = total.GetParameter(0);
  summary.bgGaussAmp1 = total.GetParameter(1);
  summary.bgGaussMean1 = total.GetParameter(2);
  summary.bgGaussSigma1 = total.GetParameter(3);
  summary.bgGaussAmp2 = total.GetParameter(4);
  summary.bgGaussMean2 = total.GetParameter(5);
  summary.bgGaussSigma2 = total.GetParameter(6);
  summary.thresholdNorm = total.GetParameter(7);
  summary.thresholdPower = total.GetParameter(8);
  summary.thresholdSlope = total.GetParameter(9);
  return summary;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput = argc > 1 ? argv[1] : "PhiSignalOnlyHistograms.root";
  const std::string sbInput = argc > 2 ? argv[2] : "PhiSBHistograms.root";
  const std::string outputDir = argc > 3 ? argv[3] : "SBFitResults";

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile signalFile(signalInput.c_str(), "READ");
  TFile sbFile(sbInput.c_str(), "READ");
  if (signalFile.IsZombie() || sbFile.IsZombie()) {
    std::cerr << "Error opening input ROOT files" << std::endl;
    return 1;
  }

  TH1D* hSignal1 = nullptr;
  TH1D* hSignal2 = nullptr;
  TH1D* hSB1 = nullptr;
  TH1D* hSB2 = nullptr;
  signalFile.GetObject("hPhiMass1Tag", hSignal1);
  signalFile.GetObject("hPhiMass2Tag", hSignal2);
  sbFile.GetObject("hPhiSBMass1Tag", hSB1);
  sbFile.GetObject("hPhiSBMass2Tag", hSB2);
  if (hSignal1 == nullptr || hSignal2 == nullptr || hSB1 == nullptr || hSB2 == nullptr) {
    std::cerr << "Missing required histograms" << std::endl;
    return 1;
  }

  FitSummary s1 = fitCategory(hSignal1, hSB1, "1tag", outputDir);
  FitSummary s2 = fitCategory(hSignal2, hSB2, "2tag", outputDir);

  std::ofstream out(outputDir + "/phi_sb_fit_summary.csv");
  out << "category,chi2,ndf,chi2ndf,signalAmp,bgGaussAmp1,bgGaussMean1,bgGaussSigma1,bgGaussAmp2,bgGaussMean2,bgGaussSigma2,thresholdNorm,thresholdPower,thresholdSlope\n";
  for (const FitSummary& s : {s1, s2})
    out << s.category << "," << s.chi2 << "," << s.ndf << "," << s.chi2ndf << ","
        << s.signalAmp << "," << s.bgGaussAmp1 << "," << s.bgGaussMean1 << "," << s.bgGaussSigma1 << ","
        << s.bgGaussAmp2 << "," << s.bgGaussMean2 << "," << s.bgGaussSigma2 << ","
        << s.thresholdNorm << "," << s.thresholdPower << "," << s.thresholdSlope << "\n";
  out.close();

  std::cout << "1tag chi2/ndf = " << s1.chi2ndf << std::endl;
  std::cout << "2tag chi2/ndf = " << s2.chi2ndf << std::endl;
  return 0;
}
