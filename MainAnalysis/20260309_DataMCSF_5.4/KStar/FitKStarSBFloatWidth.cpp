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
constexpr double kPionMass = 0.13957039;
constexpr double kThreshold = kKaonMass + kPionMass;
constexpr double kFitMin = 0.70;
constexpr double kFitMax = 1.10;
constexpr int kDisplayRebin = 4;

struct SignalShape {
  double mean = 0.8955;
  double sigma = 0.020;
  double alphaL = 1.5;
  double nL = 5.0;
  double alphaR = 1.5;
  double nR = 5.0;
  double gaussFrac = 0.2;
  double gaussSigma = 0.05;
};

SignalShape gSignalShape;

double DoubleSidedCrystalBallUnit(double x, double mean, double sigma,
                                  double alphaL, double nL,
                                  double alphaR, double nR) {
  if (sigma <= 0.0 || alphaL <= 0.0 || nL <= 1.0 || alphaR <= 0.0 || nR <= 1.0) return 0.0;
  const double t = (x - mean) / sigma;
  if (t > -alphaL && t < alphaR) return std::exp(-0.5 * t * t);
  if (t <= -alphaL) {
    const double A = std::pow(nL / alphaL, nL) * std::exp(-0.5 * alphaL * alphaL);
    const double B = nL / alphaL - alphaL;
    return A / std::pow(B - t, nL);
  }
  const double A = std::pow(nR / alphaR, nR) * std::exp(-0.5 * alphaR * alphaR);
  const double B = nR / alphaR - alphaR;
  return A / std::pow(B + t, nR);
}

double SignalOnlyShape(double* x, double* p) {
  const double xx = x[0];
  const double amp = p[0];
  const double widthScale = p[1];
  const double dscb =
      DoubleSidedCrystalBallUnit(xx, gSignalShape.mean, gSignalShape.sigma * widthScale,
                                 gSignalShape.alphaL, gSignalShape.nL,
                                 gSignalShape.alphaR, gSignalShape.nR);
  const double gauss =
      std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / (gSignalShape.gaussSigma * widthScale), 2));
  return amp * (dscb + gSignalShape.gaussFrac * gauss);
}

double ThresholdExp1(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) * std::exp(p[2] * xx);
}

double TotalThresholdExp1FloatWidth(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[3] = {p[2], p[3], p[4]};
  return SignalOnlyShape(x, sp) + ThresholdExp1(x, bp);
}

SignalShape deriveSignalShape(TH1D* h) {
  TF1 f("fSignalShapeKStarFloat", [](double* x, double* p) {
    const double xx = x[0];
    const double mean = p[1];
    const double sigma = p[2];
    const double alphaL = p[3];
    const double nL = p[4];
    const double alphaR = p[5];
    const double nR = p[6];
    const double gaussNorm = p[7];
    const double gaussSigma = p[8];
    const double dscb = p[0] * DoubleSidedCrystalBallUnit(xx, mean, sigma, alphaL, nL, alphaR, nR);
    const double gauss = gaussNorm * std::exp(-0.5 * std::pow((xx - mean) / gaussSigma, 2));
    return dscb + gauss;
  }, kFitMin, kFitMax, 9);

  f.SetParameters(std::max(100.0, 0.8 * h->GetMaximum()), 0.8955, 0.020, 1.5, 5.0, 1.5, 5.0,
                  std::max(30.0, 0.2 * h->GetMaximum()), 0.05);
  f.SetParLimits(1, 0.86, 0.93);
  f.SetParLimits(2, 0.002, 0.10);
  f.SetParLimits(3, 0.2, 8.0);
  f.SetParLimits(4, 1.2, 80.0);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  f.SetParLimits(8, 0.002, 0.20);
  h->Fit(&f, "RQ0");

  SignalShape shape;
  shape.mean = f.GetParameter(1);
  shape.sigma = f.GetParameter(2);
  shape.alphaL = f.GetParameter(3);
  shape.nL = f.GetParameter(4);
  shape.alphaR = f.GetParameter(5);
  shape.nR = f.GetParameter(6);
  shape.gaussFrac = (f.GetParameter(0) != 0.0) ? f.GetParameter(7) / f.GetParameter(0) : 0.0;
  shape.gaussSigma = f.GetParameter(8);
  return shape;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput = argc > 1 ? argv[1] : "KStarCombinedAssignmentHistograms.root";
  const std::string sbInput = argc > 2 ? argv[2] : "KStarSBHistograms.root";
  const std::string outputDir = argc > 3 ? argv[3] : "SBFitResultsFloatWidth";

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile signalFile(signalInput.c_str(), "READ");
  TFile sbFile(sbInput.c_str(), "READ");
  if (signalFile.IsZombie() || sbFile.IsZombie()) {
    std::cerr << "Error opening input files" << std::endl;
    return 1;
  }

  TH1D* hSignal = nullptr;
  TH1D* hSB = nullptr;
  signalFile.GetObject("hKStarMassKaonTag", hSignal);
  sbFile.GetObject("hKStarSBMassKaonTag", hSB);
  if (hSignal == nullptr || hSB == nullptr) {
    std::cerr << "Error: required kaon-tag histograms not found" << std::endl;
    return 1;
  }

  gSignalShape = deriveSignalShape(hSignal);

  TF1 total("fTotalKaonTagFloatWidth", TotalThresholdExp1FloatWidth, kFitMin, kFitMax, 5);
  total.SetParNames("S", "WidthScale", "N", "p", "b1");
  total.SetParameters(0.08 * hSB->GetMaximum(), 1.0, 0.9 * hSB->GetMaximum(), 0.8, -4.0);
  total.SetParLimits(0, 0.0, 1e9);
  total.SetParLimits(1, 0.85, 1.15);
  total.SetParLimits(2, 0.0, 1e12);
  total.SetParLimits(3, 0.0, 10.0);
  total.SetParLimits(4, -30.0, 10.0);
  hSB->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hSB->Clone("hKStarSBKaonTagFloatWidthDisp"));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double displayScale = hDisp->GetXaxis()->GetBinWidth(1) / hSB->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw("fSignalKaonTagFloatWidth", SignalOnlyShape, kFitMin, kFitMax, 2);
  signalDraw.SetParameters(total.GetParameter(0) * displayScale, total.GetParameter(1));
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 totalDraw(total);
  totalDraw.SetParameter(0, totalDraw.GetParameter(0) * displayScale);
  totalDraw.SetParameter(2, totalDraw.GetParameter(2) * displayScale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TF1 backgroundDraw("fBackgroundKaonTagFloatWidth", ThresholdExp1, kFitMin, kFitMax, 3);
  backgroundDraw.SetParameters(total.GetParameter(2) * displayScale, total.GetParameter(3), total.GetParameter(4));
  backgroundDraw.SetLineColor(kGreen + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(2);

  TCanvas c("cKaonTagFloatWidth", "cKaonTagFloatWidth", 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle("kaon_tag reco-only S+B fit with floating width");
  hDisp->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hDisp->GetYaxis()->SetTitle("Assignments / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  signalDraw.Draw("same");
  backgroundDraw.Draw("same");

  TLegend leg(0.52, 0.63, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event assignments", "lep");
  leg.AddEntry(&totalDraw, "Signal + background", "l");
  leg.AddEntry(&signalDraw, "Signal with width scale", "l");
  leg.AddEntry(&backgroundDraw, "Threshold background", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.14, 0.44, 0.32, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("Width scale = %.3f", total.GetParameter(1)));
  txt.AddText(Form("Signal mean fixed at %.5f GeV", gSignalShape.mean));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull("hPullKaonTagFloatWidth", "", hDisp->GetNbinsX(),
             hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.GetYaxis()->SetTitle("Pull");
  hPull.GetXaxis()->SetTitle("m(K#pi) [GeV]");
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

  c.SaveAs((outputDir + "/kstar_sb_kaon_tag_thresholdexp1_floatwidth_fit.pdf").c_str());

  std::ofstream out(outputDir + "/kstar_sb_kaon_tag_thresholdexp1_floatwidth_summary.csv");
  out << "category,model,chi2,ndf,chi2ndf,signalAmp,widthScale,thresholdNorm,thresholdPower,thresholdSlope\n";
  out << "kaon_tag,ThresholdExp1FloatWidth," << total.GetChisquare() << "," << total.GetNDF() << ","
      << total.GetChisquare() / total.GetNDF() << "," << total.GetParameter(0) << ","
      << total.GetParameter(1) << "," << total.GetParameter(2) << ","
      << total.GetParameter(3) << "," << total.GetParameter(4) << "\n";
  out.close();

  std::cout << "kaon_tag ThresholdExp1 with floating width: chi2/ndf = "
            << total.GetChisquare() / total.GetNDF()
            << ", width scale = " << total.GetParameter(1) << std::endl;
  std::cout << "Wrote outputs to " << outputDir << std::endl;
  return 0;
}
