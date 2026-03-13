#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
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
constexpr double kPionMass = 0.13957039;
constexpr double kThreshold = 2.0 * kPionMass;
constexpr double kFitMin = 0.30;
constexpr double kFitMax = 0.60;
constexpr int kDisplayRebin = 4;

struct SignalShape {
  double mean = 0.4976;
  double sigma = 0.008;
  double alphaL = 1.5;
  double nL = 5.0;
  double alphaR = 1.5;
  double nR = 5.0;
  double gaussFrac = 0.2;
  double gaussSigma = 0.03;
};

struct FitSummary {
  std::string category;
  std::string model;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 1e9;
  double signalAmp = 0.0;
  double widthScale = 1.0;
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
  const double dscb = DoubleSidedCrystalBallUnit(xx, gSignalShape.mean, gSignalShape.sigma * widthScale,
                                                 gSignalShape.alphaL, gSignalShape.nL,
                                                 gSignalShape.alphaR, gSignalShape.nR);
  const double gauss =
      std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / (gSignalShape.gaussSigma * widthScale), 2));
  return amp * (dscb + gSignalShape.gaussFrac * gauss);
}

double ThresholdExp2(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) * std::exp(p[2] * xx + p[3] * xx * xx);
}

double ThresholdExp3(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) *
         std::exp(p[2] * xx + p[3] * xx * xx + p[4] * xx * xx * xx);
}

double GaussThresholdExp2(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double tp[4] = {p[3], p[4], p[5], p[6]};
  return g + ThresholdExp2(x, tp);
}

double DoubleGaussThresholdExp2(double* x, double* p) {
  const double xx = x[0];
  const double g1 = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  const double g2 = p[3] * std::exp(-0.5 * std::pow((xx - p[4]) / p[5], 2));
  double tp[4] = {p[6], p[7], p[8], p[9]};
  return g1 + g2 + ThresholdExp2(x, tp);
}

double GaussThresholdExp3(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double tp[5] = {p[3], p[4], p[5], p[6], p[7]};
  return g + ThresholdExp3(x, tp);
}

double DoubleGaussThresholdExp3(double* x, double* p) {
  const double xx = x[0];
  const double g1 = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  const double g2 = p[3] * std::exp(-0.5 * std::pow((xx - p[4]) / p[5], 2));
  double tp[5] = {p[6], p[7], p[8], p[9], p[10]};
  return g1 + g2 + ThresholdExp3(x, tp);
}

double TotalThresholdOnly(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[4] = {p[2], p[3], p[4], p[5]};
  return SignalOnlyShape(x, sp) + ThresholdExp2(x, bp);
}

double TotalGaussThreshold(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[7] = {p[2], p[3], p[4], p[5], p[6], p[7], p[8]};
  return SignalOnlyShape(x, sp) + GaussThresholdExp2(x, bp);
}

double TotalDoubleGaussThreshold(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[10] = {p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]};
  return SignalOnlyShape(x, sp) + DoubleGaussThresholdExp2(x, bp);
}

double TotalThresholdOnlyExp3(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[5] = {p[2], p[3], p[4], p[5], p[6]};
  return SignalOnlyShape(x, sp) + ThresholdExp3(x, bp);
}

double TotalGaussThresholdExp3(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[8] = {p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]};
  return SignalOnlyShape(x, sp) + GaussThresholdExp3(x, bp);
}

double TotalDoubleGaussThresholdExp3(double* x, double* p) {
  double sp[2] = {p[0], p[1]};
  double bp[11] = {p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12]};
  return SignalOnlyShape(x, sp) + DoubleGaussThresholdExp3(x, bp);
}

SignalShape deriveSignalShape(TH1D* h) {
  TF1 f("fSignalShapeKShort", [](double* x, double* p) {
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
  }, 0.35, 0.70, 9);

  f.SetParameters(std::max(60.0, 0.8 * h->GetMaximum()), 0.4976, 0.008, 1.5, 5.0, 1.5, 5.0,
                  std::max(20.0, 0.2 * h->GetMaximum()), 0.03);
  f.SetParLimits(1, 0.485, 0.510);
  f.SetParLimits(2, 0.001, 0.05);
  f.SetParLimits(3, 0.2, 8.0);
  f.SetParLimits(4, 1.2, 80.0);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  f.SetParLimits(8, 0.001, 0.15);
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

TF1 buildTotalModel(const std::string& model, const std::string& name, TH1D* hSB) {
  const double maxY = std::max(1000.0, hSB->GetMaximum());
  TF1 f;
  if (model == "ThresholdExp2") {
    f = TF1(name.c_str(), TotalThresholdOnly, kFitMin, kFitMax, 6);
    f.SetParNames("S", "WidthScale", "N", "p", "b1", "b2");
    f.SetParameters(0.10 * maxY, 1.0, 0.8 * maxY, 0.8, -4.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e12);
    f.SetParLimits(3, 0.0, 10.0);
    f.SetParLimits(4, -50.0, 20.0);
    f.SetParLimits(5, -50.0, 50.0);
  } else if (model == "GaussPlusThresholdExp2") {
    f = TF1(name.c_str(), TotalGaussThreshold, kFitMin, kFitMax, 9);
    f.SetParNames("S", "WidthScale", "Gamp", "Gmean", "Gsigma", "N", "p", "b1", "b2");
    f.SetParameters(0.10 * maxY, 1.0, 0.15 * maxY, 0.75, 0.10, 0.6 * maxY, 0.8, -4.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e9);
    f.SetParLimits(3, 0.55, 0.95);
    f.SetParLimits(4, 0.02, 0.30);
    f.SetParLimits(5, 0.0, 1e12);
    f.SetParLimits(6, 0.0, 10.0);
    f.SetParLimits(7, -50.0, 20.0);
    f.SetParLimits(8, -50.0, 50.0);
  } else if (model == "DoubleGaussPlusThresholdExp2") {
    f = TF1(name.c_str(), TotalDoubleGaussThreshold, kFitMin, kFitMax, 12);
    const double params[12] = {0.10 * maxY, 1.0,
                               0.10 * maxY, 0.35, 0.06,
                               0.12 * maxY, 0.80, 0.10,
                               0.5 * maxY, 0.8, -4.0, 0.0};
    f.SetParameters(params);
    const char* names[12] = {"S", "WidthScale", "G1amp", "G1mean", "G1sigma", "G2amp",
                             "G2mean", "G2sigma", "N", "p", "b1", "b2"};
    for (int i = 0; i < 12; ++i)
      f.SetParName(i, names[i]);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e9);
    f.SetParLimits(3, 0.28, 0.45);
    f.SetParLimits(4, 0.02, 0.12);
    f.SetParLimits(5, 0.0, 1e9);
    f.SetParLimits(6, 0.60, 0.95);
    f.SetParLimits(7, 0.03, 0.25);
    f.SetParLimits(8, 0.0, 1e12);
    f.SetParLimits(9, 0.0, 10.0);
    f.SetParLimits(10, -50.0, 20.0);
    f.SetParLimits(11, -50.0, 50.0);
  } else if (model == "ThresholdExp3") {
    f = TF1(name.c_str(), TotalThresholdOnlyExp3, kFitMin, kFitMax, 7);
    f.SetParNames("S", "WidthScale", "N", "p", "b1", "b2", "b3");
    f.SetParameters(0.10 * maxY, 1.0, 0.8 * maxY, 0.8, -4.0, 0.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e12);
    f.SetParLimits(3, 0.0, 10.0);
    f.SetParLimits(4, -50.0, 20.0);
    f.SetParLimits(5, -50.0, 50.0);
    f.SetParLimits(6, -100.0, 100.0);
  } else if (model == "GaussPlusThresholdExp3") {
    f = TF1(name.c_str(), TotalGaussThresholdExp3, kFitMin, kFitMax, 10);
    f.SetParNames("S", "WidthScale", "Gamp", "Gmean", "Gsigma", "N", "p", "b1", "b2", "b3");
    f.SetParameters(0.10 * maxY, 1.0, 0.15 * maxY, 0.75, 0.10, 0.6 * maxY, 0.8, -4.0, 0.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e9);
    f.SetParLimits(3, 0.55, 0.95);
    f.SetParLimits(4, 0.02, 0.30);
    f.SetParLimits(5, 0.0, 1e12);
    f.SetParLimits(6, 0.0, 10.0);
    f.SetParLimits(7, -50.0, 20.0);
    f.SetParLimits(8, -50.0, 50.0);
    f.SetParLimits(9, -100.0, 100.0);
  } else {
    f = TF1(name.c_str(), TotalDoubleGaussThresholdExp3, kFitMin, kFitMax, 13);
    const double params[13] = {0.10 * maxY, 1.0,
                               0.10 * maxY, 0.35, 0.06,
                               0.12 * maxY, 0.80, 0.10,
                               0.5 * maxY, 0.8, -4.0, 0.0, 0.0};
    f.SetParameters(params);
    const char* names[13] = {"S", "WidthScale", "G1amp", "G1mean", "G1sigma", "G2amp",
                             "G2mean", "G2sigma", "N", "p", "b1", "b2", "b3"};
    for (int i = 0; i < 13; ++i)
      f.SetParName(i, names[i]);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.7, 1.3);
    f.SetParLimits(2, 0.0, 1e9);
    f.SetParLimits(3, 0.28, 0.45);
    f.SetParLimits(4, 0.02, 0.12);
    f.SetParLimits(5, 0.0, 1e9);
    f.SetParLimits(6, 0.60, 0.95);
    f.SetParLimits(7, 0.03, 0.25);
    f.SetParLimits(8, 0.0, 1e12);
    f.SetParLimits(9, 0.0, 10.0);
    f.SetParLimits(10, -50.0, 20.0);
    f.SetParLimits(11, -50.0, 50.0);
    f.SetParLimits(12, -100.0, 100.0);
  }
  return f;
}

void scaleDisplayedFunction(TF1& f, const std::string& model, double scale) {
  f.SetParameter(0, f.GetParameter(0) * scale);
  if (model == "ThresholdExp2") {
    f.SetParameter(2, f.GetParameter(2) * scale);
  } else if (model == "GaussPlusThresholdExp2") {
    f.SetParameter(2, f.GetParameter(2) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
  } else if (model == "DoubleGaussPlusThresholdExp2") {
    f.SetParameter(2, f.GetParameter(2) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
    f.SetParameter(8, f.GetParameter(8) * scale);
  } else if (model == "ThresholdExp3") {
    f.SetParameter(2, f.GetParameter(2) * scale);
  } else if (model == "GaussPlusThresholdExp3") {
    f.SetParameter(2, f.GetParameter(2) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
  } else if (model == "DoubleGaussPlusThresholdExp3") {
    f.SetParameter(2, f.GetParameter(2) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
    f.SetParameter(8, f.GetParameter(8) * scale);
  }
}

TF1 buildBackgroundDraw(const std::string& model, const std::string& name, TF1& total, double scale) {
  TF1 f;
  if (model == "ThresholdExp2") {
    f = TF1(name.c_str(), ThresholdExp2, kFitMin, kFitMax, 4);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4), total.GetParameter(5));
  } else if (model == "GaussPlusThresholdExp2") {
    f = TF1(name.c_str(), GaussThresholdExp2, kFitMin, kFitMax, 7);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5) * scale, total.GetParameter(6), total.GetParameter(7), total.GetParameter(8));
  } else {
    f = TF1(name.c_str(), DoubleGaussThresholdExp2, kFitMin, kFitMax, 10);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5) * scale, total.GetParameter(6), total.GetParameter(7),
                    total.GetParameter(8) * scale, total.GetParameter(9), total.GetParameter(10), total.GetParameter(11));
  }
  if (model == "ThresholdExp3") {
    f = TF1(name.c_str(), ThresholdExp3, kFitMin, kFitMax, 5);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5), total.GetParameter(6));
  } else if (model == "GaussPlusThresholdExp3") {
    f = TF1(name.c_str(), GaussThresholdExp3, kFitMin, kFitMax, 8);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5) * scale, total.GetParameter(6), total.GetParameter(7),
                    total.GetParameter(8), total.GetParameter(9));
  } else if (model == "DoubleGaussPlusThresholdExp3") {
    f = TF1(name.c_str(), DoubleGaussThresholdExp3, kFitMin, kFitMax, 11);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5) * scale, total.GetParameter(6), total.GetParameter(7),
                    total.GetParameter(8) * scale, total.GetParameter(9), total.GetParameter(10),
                    total.GetParameter(11), total.GetParameter(12));
  } else if (model == "DoubleGaussPlusThresholdExp2") {
    f = TF1(name.c_str(), DoubleGaussThresholdExp2, kFitMin, kFitMax, 10);
    f.SetParameters(total.GetParameter(2) * scale, total.GetParameter(3), total.GetParameter(4),
                    total.GetParameter(5) * scale, total.GetParameter(6), total.GetParameter(7),
                    total.GetParameter(8) * scale, total.GetParameter(9), total.GetParameter(10), total.GetParameter(11));
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

  TF1 signalDraw(("fSignal_" + category + "_" + model).c_str(), SignalOnlyShape, kFitMin, kFitMax, 2);
  signalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  signalDraw.SetParameter(1, total.GetParameter(1));
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
  hDisp->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Pairs / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  signalDraw.Draw("same");
  backgroundDraw.Draw("same");

  TLegend leg(0.54, 0.64, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event pairs", "lep");
  leg.AddEntry(&totalDraw, "Signal + background", "l");
  leg.AddEntry(&signalDraw, "Fixed signal shape", "l");
  leg.AddEntry(&backgroundDraw, "Background component", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.14, 0.44, 0.34, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("Model: %s", model.c_str()));
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
  TH1D hPull(("hPull_" + category + "_" + model).c_str(), "", hDisp->GetNbinsX(),
             hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
  hPull.SetStats(0);
  hPull.SetMarkerStyle(20);
  hPull.SetMarkerSize(0.65);
  hPull.GetYaxis()->SetTitle("Pull");
  hPull.GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
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

  c.SaveAs((outputDir + "/kshort_sb_" + category + "_" + model + "_fit.pdf").c_str());
  delete hDisp;

  FitSummary summary;
  summary.category = category;
  summary.model = model;
  summary.chi2 = total.GetChisquare();
  summary.ndf = total.GetNDF();
  summary.chi2ndf = total.GetChisquare() / total.GetNDF();
  summary.signalAmp = total.GetParameter(0);
  summary.widthScale = total.GetParameter(1);
  return summary;
}

FitSummary fitCategory(TH1D* hSignal, TH1D* hSB, const std::string& category,
                       const std::string& outputDir, std::ofstream& csv) {
  gSignalShape = deriveSignalShape(hSignal);

  const std::vector<std::string> models = {
      "ThresholdExp2", "GaussPlusThresholdExp2", "DoubleGaussPlusThresholdExp2",
      "ThresholdExp3", "GaussPlusThresholdExp3", "DoubleGaussPlusThresholdExp3"};
  std::vector<FitSummary> results;
  for (const std::string& model : models)
    results.push_back(runFit(hSB, category, model, outputDir));

  for (const FitSummary& r : results) {
    csv << r.category << "," << r.model << "," << r.chi2 << "," << r.ndf << ","
        << r.chi2ndf << "," << r.signalAmp << "," << r.widthScale << ","
        << (r.model == std::min_element(results.begin(), results.end(),
          [](const FitSummary& a, const FitSummary& b) { return a.chi2ndf < b.chi2ndf; })->model ? "yes" : "no")
        << "\n";
  }

  return *std::min_element(results.begin(), results.end(),
                           [](const FitSummary& a, const FitSummary& b) {
                             return a.chi2ndf < b.chi2ndf;
                           });
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput = argc > 1 ? argv[1] : "KShortSignalOnlyHistograms.root";
  const std::string sbInput = argc > 2 ? argv[2] : "KShortSBHistograms.root";
  const std::string outputDir = argc > 3 ? argv[3] : "SBFitResults";

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
  signalFile.GetObject("hKShortMass1Tag", hSignal1);
  signalFile.GetObject("hKShortMass2Tag", hSignal2);
  sbFile.GetObject("hKShortSBMass1Tag", hSB1);
  sbFile.GetObject("hKShortSBMass2Tag", hSB2);
  if (hSignal1 == nullptr || hSignal2 == nullptr || hSB1 == nullptr || hSB2 == nullptr) {
    std::cerr << "Error: required histograms not found" << std::endl;
    return 1;
  }

  std::ofstream out(outputDir + "/kshort_sb_fit_summary.csv");
  out << "category,model,chi2,ndf,chi2ndf,signalAmp,widthScale,isBest\n";

  FitSummary s1 = fitCategory(hSignal1, hSB1, "1tag", outputDir, out);
  FitSummary s2 = fitCategory(hSignal2, hSB2, "2tag", outputDir, out);
  out.close();

  std::cout << "1tag best: " << s1.model << " chi2/ndf = " << s1.chi2ndf << std::endl;
  std::cout << "2tag best: " << s2.model << " chi2/ndf = " << s2.chi2ndf << std::endl;
  std::cout << "Wrote outputs to " << outputDir << std::endl;
  return 0;
}
