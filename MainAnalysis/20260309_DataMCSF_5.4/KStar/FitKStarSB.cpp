#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
constexpr int kDisplayRebin = 4;
double gFitMin = 0.70;
double gFitMax = 1.10;

struct SignalShape {
  double mean = 0.8955;
  double sigma = 0.020;
  double alphaL = 1.5;
  double nL = 5.0;
  double alphaR = 1.5;
  double nR = 5.0;
  double gaussFrac = 0.2;
  double gaussSigma = 0.05;
  double gaussFrac2 = 0.0;
  double gaussSigma2 = 0.0;
  double gaussFrac3 = 0.0;
  double gaussSigma3 = 0.0;
};

struct FitSummary {
  std::string category;
  std::string model;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 1e9;
  double signalAmp = 0.0;
};

SignalShape gSignalShape;
std::string gActiveSignalModel = "DoubleSidedCBPlusGauss";

std::string getArgument(int argc, char* argv[], const std::string& option,
                        const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i)
    if (argv[i] == option) return argv[i + 1];
  return defaultValue;
}

double getDoubleArgument(int argc, char* argv[], const std::string& option, double defaultValue) {
  const std::string value = getArgument(argc, argv, option, "");
  return value.empty() ? defaultValue : std::stod(value);
}

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
  if (gActiveSignalModel == "QuadGaussian") {
    const double g1 =
        std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.sigma, 2));
    const double g2 =
        std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.gaussSigma, 2));
    const double g3 =
        std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.gaussSigma2, 2));
    const double g4 =
        std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.gaussSigma3, 2));
    return amp * (g1 + gSignalShape.gaussFrac * g2 +
                  gSignalShape.gaussFrac2 * g3 + gSignalShape.gaussFrac3 * g4);
  }
  const double dscb =
      DoubleSidedCrystalBallUnit(xx, gSignalShape.mean, gSignalShape.sigma,
                                 gSignalShape.alphaL, gSignalShape.nL,
                                 gSignalShape.alphaR, gSignalShape.nR);
  const double gauss =
      std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / gSignalShape.gaussSigma, 2));
  return amp * (dscb + gSignalShape.gaussFrac * gauss);
}

double ThresholdExp1(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) * std::exp(p[2] * xx);
}

double ThresholdExp2(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) * std::exp(p[2] * xx + p[3] * xx * xx);
}

double GaussThresholdExp1(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double tp[3] = {p[3], p[4], p[5]};
  return g + ThresholdExp1(x, tp);
}

double Gauss080ThresholdExp1(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double tp[3] = {p[3], p[4], p[5]};
  return g + ThresholdExp1(x, tp);
}

double GaussThresholdExp2(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double tp[4] = {p[3], p[4], p[5], p[6]};
  return g + ThresholdExp2(x, tp);
}

double DoubleGaussThresholdExp1(double* x, double* p) {
  const double xx = x[0];
  const double g1 = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  const double g2 = p[3] * std::exp(-0.5 * std::pow((xx - p[4]) / p[5], 2));
  double tp[3] = {p[6], p[7], p[8]};
  return g1 + g2 + ThresholdExp1(x, tp);
}

double TotalThresholdExp1(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[3] = {p[1], p[2], p[3]};
  return SignalOnlyShape(x, sp) + ThresholdExp1(x, bp);
}

double TotalThresholdExp2(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[4] = {p[1], p[2], p[3], p[4]};
  return SignalOnlyShape(x, sp) + ThresholdExp2(x, bp);
}

double TotalGaussThresholdExp1(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[6] = {p[1], p[2], p[3], p[4], p[5], p[6]};
  return SignalOnlyShape(x, sp) + GaussThresholdExp1(x, bp);
}

double TotalGauss080ThresholdExp1(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[6] = {p[1], p[2], p[3], p[4], p[5], p[6]};
  return SignalOnlyShape(x, sp) + Gauss080ThresholdExp1(x, bp);
}

double TotalGaussThresholdExp2(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[7] = {p[1], p[2], p[3], p[4], p[5], p[6], p[7]};
  return SignalOnlyShape(x, sp) + GaussThresholdExp2(x, bp);
}

double TotalDoubleGaussThresholdExp1(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[9] = {p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]};
  return SignalOnlyShape(x, sp) + DoubleGaussThresholdExp1(x, bp);
}

SignalShape deriveSignalShape(TH1D* h, const std::string& signalModel) {
  SignalShape shape;
  gActiveSignalModel = signalModel;
  if (signalModel == "QuadGaussian") {
    TF1 f("fSignalShapeKStarQuadG",
          "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
          "[5]*exp(-0.5*((x-[1])/[6])^2)+[7]*exp(-0.5*((x-[1])/[8])^2)",
          gFitMin, gFitMax);
    f.SetParameters(std::max(100.0, 0.40 * h->GetMaximum()), 0.8955, 0.012,
                    std::max(50.0, 0.28 * h->GetMaximum()), 0.025,
                    std::max(30.0, 0.20 * h->GetMaximum()), 0.045,
                    std::max(10.0, 0.12 * h->GetMaximum()), 0.080);
    f.SetParLimits(1, 0.86, 0.93);
    f.SetParLimits(2, 0.002, 0.08);
    f.SetParLimits(4, 0.002, 0.12);
    f.SetParLimits(6, 0.002, 0.18);
    f.SetParLimits(8, 0.002, 0.24);
    h->Fit(&f, "RQ0");

    shape.mean = f.GetParameter(1);
    shape.sigma = f.GetParameter(2);
    shape.gaussSigma = f.GetParameter(4);
    shape.gaussSigma2 = f.GetParameter(6);
    shape.gaussSigma3 = f.GetParameter(8);
    shape.gaussFrac = (f.GetParameter(0) != 0.0) ? f.GetParameter(3) / f.GetParameter(0) : 0.0;
    shape.gaussFrac2 = (f.GetParameter(0) != 0.0) ? f.GetParameter(5) / f.GetParameter(0) : 0.0;
    shape.gaussFrac3 = (f.GetParameter(0) != 0.0) ? f.GetParameter(7) / f.GetParameter(0) : 0.0;
    shape.alphaL = shape.alphaR = 1.0;
    shape.nL = shape.nR = 2.0;
    return shape;
  }
  if (signalModel == "DoubleSidedCB") {
    TF1 f("fSignalShapeKStarDSCB", [](double* x, double* p) {
      return p[0] * DoubleSidedCrystalBallUnit(x[0], p[1], p[2], p[3], p[4], p[5], p[6]);
    }, gFitMin, gFitMax, 7);

    f.SetParameters(std::max(100.0, 0.8 * h->GetMaximum()), 0.8955, 0.020, 1.5, 5.0, 1.5, 5.0);
    f.SetParLimits(1, 0.86, 0.93);
    f.SetParLimits(2, 0.002, 0.10);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
    h->Fit(&f, "RQ0");

    shape.mean = f.GetParameter(1);
    shape.sigma = f.GetParameter(2);
    shape.alphaL = f.GetParameter(3);
    shape.nL = f.GetParameter(4);
    shape.alphaR = f.GetParameter(5);
    shape.nR = f.GetParameter(6);
    shape.gaussFrac = 0.0;
    shape.gaussSigma = shape.sigma;
    return shape;
  }

  TF1 f("fSignalShapeKStarDSCBG", [](double* x, double* p) {
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
  }, gFitMin, gFitMax, 9);

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

TF1 buildTotalModel(const std::string& model, const std::string& name, TH1D* hSB) {
  const double maxY = std::max(1000.0, hSB->GetMaximum());
  TF1 f;
  if (model == "ThresholdExp1") {
    f = TF1(name.c_str(), TotalThresholdExp1, gFitMin, gFitMax, 4);
    f.SetParNames("S", "N", "p", "b1");
    f.SetParameters(0.08 * maxY, 0.9 * maxY, 0.8, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e12);
    f.SetParLimits(2, 0.0, 10.0);
    f.SetParLimits(3, -30.0, 10.0);
  } else if (model == "ThresholdExp2") {
    f = TF1(name.c_str(), TotalThresholdExp2, gFitMin, gFitMax, 5);
    f.SetParNames("S", "N", "p", "b1", "b2");
    f.SetParameters(0.08 * maxY, 0.8 * maxY, 0.8, -4.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e12);
    f.SetParLimits(2, 0.0, 10.0);
    f.SetParLimits(3, -30.0, 10.0);
    f.SetParLimits(4, -20.0, 20.0);
  } else if (model == "GaussPlusThresholdExp1") {
    f = TF1(name.c_str(), TotalGaussThresholdExp1, gFitMin, gFitMax, 7);
    f.SetParNames("S", "Gamp", "Gmean", "Gsigma", "N", "p", "b1");
    f.SetParameters(0.08 * maxY, 0.15 * maxY, 1.02, 0.05, 0.5 * maxY, 0.8, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e9);
    f.SetParLimits(2, 0.80, 1.08);
    f.SetParLimits(3, 0.01, 0.20);
    f.SetParLimits(4, 0.0, 1e12);
    f.SetParLimits(5, 0.0, 10.0);
    f.SetParLimits(6, -30.0, 10.0);
  } else if (model == "Gauss080PlusThresholdExp1") {
    f = TF1(name.c_str(), TotalGauss080ThresholdExp1, gFitMin, gFitMax, 7);
    f.SetParNames("S", "Gamp", "Gmean", "Gsigma", "N", "p", "b1");
    f.SetParameters(0.08 * maxY, 0.12 * maxY, 0.80, 0.035, 0.5 * maxY, 0.8, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e9);
    f.SetParLimits(2, 0.76, 0.84);
    f.SetParLimits(3, 0.01, 0.10);
    f.SetParLimits(4, 0.0, 1e12);
    f.SetParLimits(5, 0.0, 10.0);
    f.SetParLimits(6, -30.0, 10.0);
  } else if (model == "GaussPlusThresholdExp2") {
    f = TF1(name.c_str(), TotalGaussThresholdExp2, gFitMin, gFitMax, 8);
    f.SetParNames("S", "Gamp", "Gmean", "Gsigma", "N", "p", "b1", "b2");
    f.SetParameters(0.08 * maxY, 0.15 * maxY, 1.02, 0.05, 0.5 * maxY, 0.8, -4.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e9);
    f.SetParLimits(2, 0.80, 1.08);
    f.SetParLimits(3, 0.01, 0.20);
    f.SetParLimits(4, 0.0, 1e12);
    f.SetParLimits(5, 0.0, 10.0);
    f.SetParLimits(6, -30.0, 10.0);
    f.SetParLimits(7, -20.0, 20.0);
  } else {
    f = TF1(name.c_str(), TotalDoubleGaussThresholdExp1, gFitMin, gFitMax, 10);
    f.SetParNames("S", "G1amp", "G1mean", "G1sigma", "G2amp", "G2mean", "G2sigma", "N", "p", "b1");
    f.SetParameters(0.08 * maxY, 0.12 * maxY, 0.78, 0.04, 0.10 * maxY, 1.03, 0.06, 0.4 * maxY, 0.8, -4.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e9);
    f.SetParLimits(2, 0.72, 0.90);
    f.SetParLimits(3, 0.01, 0.12);
    f.SetParLimits(4, 0.0, 1e9);
    f.SetParLimits(5, 0.92, 1.08);
    f.SetParLimits(6, 0.01, 0.18);
    f.SetParLimits(7, 0.0, 1e12);
    f.SetParLimits(8, 0.0, 10.0);
    f.SetParLimits(9, -30.0, 10.0);
  }
  return f;
}

void scaleDisplayedFunction(TF1& f, const std::string& model, double scale) {
  f.SetParameter(0, f.GetParameter(0) * scale);
  if (model == "ThresholdExp1") {
    f.SetParameter(1, f.GetParameter(1) * scale);
  } else if (model == "ThresholdExp2") {
    f.SetParameter(1, f.GetParameter(1) * scale);
  } else if (model == "GaussPlusThresholdExp1") {
    f.SetParameter(1, f.GetParameter(1) * scale);
    f.SetParameter(4, f.GetParameter(4) * scale);
  } else if (model == "Gauss080PlusThresholdExp1") {
    f.SetParameter(1, f.GetParameter(1) * scale);
    f.SetParameter(4, f.GetParameter(4) * scale);
  } else if (model == "GaussPlusThresholdExp2") {
    f.SetParameter(1, f.GetParameter(1) * scale);
    f.SetParameter(4, f.GetParameter(4) * scale);
  } else if (model == "DoubleGaussPlusThresholdExp1") {
    f.SetParameter(1, f.GetParameter(1) * scale);
    f.SetParameter(4, f.GetParameter(4) * scale);
    f.SetParameter(7, f.GetParameter(7) * scale);
  }
}

TF1 buildBackgroundDraw(const std::string& model, const std::string& name, TF1& total, double scale) {
  TF1 f;
  if (model == "ThresholdExp1") {
    f = TF1(name.c_str(), ThresholdExp1, gFitMin, gFitMax, 3);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3));
  } else if (model == "ThresholdExp2") {
    f = TF1(name.c_str(), ThresholdExp2, gFitMin, gFitMax, 4);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2),
                    total.GetParameter(3), total.GetParameter(4));
  } else if (model == "GaussPlusThresholdExp1") {
    f = TF1(name.c_str(), GaussThresholdExp1, gFitMin, gFitMax, 6);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3),
                    total.GetParameter(4) * scale, total.GetParameter(5), total.GetParameter(6));
  } else if (model == "Gauss080PlusThresholdExp1") {
    f = TF1(name.c_str(), Gauss080ThresholdExp1, gFitMin, gFitMax, 6);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3),
                    total.GetParameter(4) * scale, total.GetParameter(5), total.GetParameter(6));
  } else if (model == "GaussPlusThresholdExp2") {
    f = TF1(name.c_str(), GaussThresholdExp2, gFitMin, gFitMax, 7);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3),
                    total.GetParameter(4) * scale, total.GetParameter(5), total.GetParameter(6),
                    total.GetParameter(7));
  } else {
    f = TF1(name.c_str(), DoubleGaussThresholdExp1, gFitMin, gFitMax, 9);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3),
                    total.GetParameter(4) * scale, total.GetParameter(5), total.GetParameter(6),
                    total.GetParameter(7) * scale, total.GetParameter(8), total.GetParameter(9));
  }
  return f;
}

FitSummary runFit(TH1D* hSB, const std::string& category, const std::string& model,
                  const std::string& outputDir) {
  TF1 total = buildTotalModel(model, "fTotal_" + category + "_" + model, hSB);
  hSB->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hSB->Clone((std::string(hSB->GetName()) + "_disp_" + model).c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double displayScale = hDisp->GetXaxis()->GetBinWidth(1) / hSB->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw(("fSignal_" + category + "_" + model).c_str(), SignalOnlyShape, gFitMin, gFitMax, 1);
  signalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 totalDraw(total);
  scaleDisplayedFunction(totalDraw, model, displayScale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TF1 backgroundDraw = buildBackgroundDraw(model, "fBackground_" + category + "_" + model, total, displayScale);
  backgroundDraw.SetLineColor(kGreen + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(2);

  TCanvas c(("c_" + category + "_" + model).c_str(), ("c_" + category + "_" + model).c_str(), 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle((category + " reco-only S+B fit").c_str());
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
  leg.AddEntry(&signalDraw, "Fixed signal shape", "l");
  leg.AddEntry(&backgroundDraw, "Background component", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.14, 0.44, 0.30, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("Model: %s", model.c_str()));
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("Signal mean fixed at %.5f GeV", gSignalShape.mean));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull(("hPull_" + category + "_" + model).c_str(), "", hDisp->GetNbinsX(),
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
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/kstar_sb_" + category + "_" + model + "_fit.pdf").c_str());
  delete hDisp;

  FitSummary summary;
  summary.category = category;
  summary.model = model;
  summary.chi2 = total.GetChisquare();
  summary.ndf = total.GetNDF();
  summary.chi2ndf = total.GetChisquare() / total.GetNDF();
  summary.signalAmp = total.GetParameter(0);
  return summary;
}

FitSummary fitCategory(TH1D* hSignal, TH1D* hSB, const std::string& category,
                       const std::string& outputDir, std::ofstream& csv,
                       const std::string& signalModel) {
  gSignalShape = deriveSignalShape(hSignal, signalModel);

  const std::vector<std::string> models = {
      "ThresholdExp1", "ThresholdExp2", "GaussPlusThresholdExp1", "Gauss080PlusThresholdExp1",
      "GaussPlusThresholdExp2", "DoubleGaussPlusThresholdExp1"};
  std::vector<FitSummary> results;
  for (const std::string& model : models)
    results.push_back(runFit(hSB, category, model, outputDir));

  const FitSummary& best = *std::min_element(results.begin(), results.end(),
                                             [](const FitSummary& a, const FitSummary& b) {
                                               return a.chi2ndf < b.chi2ndf;
                                             });

  for (const FitSummary& r : results) {
    csv << r.category << "," << r.model << "," << r.chi2 << "," << r.ndf << ","
        << r.chi2ndf << "," << r.signalAmp << ","
        << (r.model == best.model ? "yes" : "no") << "\n";
  }
  return best;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput =
      getArgument(argc, argv, "--signal-input",
                  argc > 1 ? argv[1] : "KStarCombinedAssignmentHistograms.root");
  const std::string sbInput =
      getArgument(argc, argv, "--sb-input",
                  argc > 2 ? argv[2] : "KStarSBHistograms.root");
  const std::string outputDir =
      getArgument(argc, argv, "--output-dir",
                  argc > 3 ? argv[3] : "SBFitResults");
  const std::string signalModel =
      getArgument(argc, argv, "--signal-model", "DoubleSidedCBPlusGauss");
  gFitMin = getDoubleArgument(argc, argv, "--fit-min", gFitMin);
  gFitMax = getDoubleArgument(argc, argv, "--fit-max", gFitMax);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile signalFile(signalInput.c_str(), "READ");
  TFile sbFile(sbInput.c_str(), "READ");
  if (signalFile.IsZombie() || sbFile.IsZombie()) {
    std::cerr << "Error opening input files" << std::endl;
    return 1;
  }

  TH1D* hSignal1 = nullptr;
  TH1D* hSignal2 = nullptr;
  TH1D* hSB1 = nullptr;
  TH1D* hSB2 = nullptr;
  signalFile.GetObject("hKStarMassKaonTag", hSignal1);
  signalFile.GetObject("hKStarMassKaonPionTag", hSignal2);
  sbFile.GetObject("hKStarSBMassKaonTag", hSB1);
  sbFile.GetObject("hKStarSBMassKaonPionTag", hSB2);
  if (hSignal1 == nullptr || hSignal2 == nullptr || hSB1 == nullptr || hSB2 == nullptr) {
    std::cerr << "Error: required histograms not found" << std::endl;
    return 1;
  }

  std::ofstream out(outputDir + "/kstar_sb_fit_summary.csv");
  out << "category,model,chi2,ndf,chi2ndf,signalAmp,isBest\n";

  FitSummary s1 = fitCategory(hSignal1, hSB1, "kaon_tag", outputDir, out, signalModel);
  FitSummary s2 = fitCategory(hSignal2, hSB2, "kaon_pion_tag", outputDir, out, signalModel);
  out.close();

  std::cout << "kaon_tag best: " << s1.model << " chi2/ndf = " << s1.chi2ndf << std::endl;
  std::cout << "kaon_pion_tag best: " << s2.model << " chi2/ndf = " << s2.chi2ndf << std::endl;
  std::cout << "Wrote outputs to " << outputDir << std::endl;
  return 0;
}
