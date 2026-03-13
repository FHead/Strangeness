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
double gFitMin = 0.99;
double gFitMax = 1.06;
constexpr int kDisplayRebin = 4;

enum class SignalModel {
  GaussPlusRightTailCB,
  TripleGaussian
};

struct SignalShape {
  SignalModel model = SignalModel::GaussPlusRightTailCB;
  double mean = 1.0195;
  double sigma1 = 0.0023;
  double scale2 = 0.2;
  double sigma2 = 0.0060;
  double scale3 = 0.0;
  double sigma3 = 0.0100;
  double alpha = 1.5;
  double n = 8.0;
};

struct FitSummary {
  std::string category;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 0.0;
  double signalAmp = 0.0;
  double signalWidthScale = 1.0;
  double thresholdNorm = 0.0;
  double thresholdPower = 0.0;
  double thresholdSlope = 0.0;
};

SignalShape gSignalShape;

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

SignalModel parseSignalModel(const std::string& value) {
  if (value == "TripleGaussian" || value == "triplegaussian" || value == "triple")
    return SignalModel::TripleGaussian;
  return SignalModel::GaussPlusRightTailCB;
}

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
  const double widthScale = p[1];
  const double sigma1 = gSignalShape.sigma1 * widthScale;
  const double sigma2 = gSignalShape.sigma2 * widthScale;
  const double sigma3 = gSignalShape.sigma3 * widthScale;

  if (gSignalShape.model == SignalModel::TripleGaussian) {
    const double g1 = std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / sigma1, 2));
    const double g2 = std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / sigma2, 2));
    const double g3 = std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / sigma3, 2));
    return amp * (g1 + gSignalShape.scale2 * g2 + gSignalShape.scale3 * g3);
  }

  const double gaussian = std::exp(-0.5 * std::pow((xx - gSignalShape.mean) / sigma1, 2));
  const double cb = RightTailCBUnit(xx, gSignalShape.mean, sigma2, gSignalShape.alpha, gSignalShape.n);
  return amp * (gaussian + gSignalShape.scale2 * cb);
}

double BackgroundShape(double* x, double* p) {
  const double xx = x[0];
  if (xx <= kThreshold) return 0.0;
  return p[0] * std::pow(xx - kThreshold, p[1]) * std::exp(p[2] * xx);
}

double TotalShape(double* x, double* p) {
  double signalP[2] = {p[0], p[1]};
  double backgroundP[3] = {p[2], p[3], p[4]};
  return SignalOnlyShape(x, signalP) + BackgroundShape(x, backgroundP);
}

SignalShape deriveSignalShape(TH1D* h, SignalModel model) {
  if (model == SignalModel::TripleGaussian) {
    TF1 f("fSignalShapeDataTriple",
          "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+[5]*exp(-0.5*((x-[1])/[6])^2)",
          1.000, 1.050);
    f.SetParameters(std::max(50.0, 0.5 * h->GetMaximum()), 1.0195, 0.0018,
                    std::max(25.0, 0.3 * h->GetMaximum()), 0.0038,
                    std::max(10.0, 0.15 * h->GetMaximum()), 0.0065);
    f.SetParLimits(1, 1.015, 1.024);
    f.SetParLimits(2, 0.0003, 0.02);
    f.SetParLimits(4, 0.0003, 0.03);
    f.SetParLimits(6, 0.0003, 0.05);
    h->Fit(&f, "RQ0");

    SignalShape shape;
    shape.model = model;
    shape.mean = f.GetParameter(1);
    shape.sigma1 = f.GetParameter(2);
    shape.scale2 = (f.GetParameter(0) != 0.0) ? f.GetParameter(3) / f.GetParameter(0) : 0.0;
    shape.sigma2 = f.GetParameter(4);
    shape.scale3 = (f.GetParameter(0) != 0.0) ? f.GetParameter(5) / f.GetParameter(0) : 0.0;
    shape.sigma3 = f.GetParameter(6);
    return shape;
  }

  TF1 f("fSignalShapeData", [](double* x, double* p) {
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
  shape.model = model;
  shape.mean = f.GetParameter(1);
  shape.sigma1 = f.GetParameter(2);
  shape.scale2 = (f.GetParameter(0) != 0.0) ? f.GetParameter(3) / f.GetParameter(0) : 0.0;
  shape.sigma2 = f.GetParameter(4);
  shape.alpha = f.GetParameter(5);
  shape.n = f.GetParameter(6);
  return shape;
}

FitSummary fitCategory(TH1D* hSignal, TH1D* hData, const std::string& category,
                       const std::string& outputDir, SignalModel signalModel) {
  gSignalShape = deriveSignalShape(hSignal, signalModel);

  TF1 total(("fTotalData_" + category).c_str(), TotalShape, gFitMin, gFitMax, 5);
  total.SetParNames("S", "WidthScale", "N", "p", "b1");
  total.SetParameters(std::max(100.0, 0.3 * hData->GetMaximum()),
                      1.0, std::max(100.0, 0.2 * hData->GetMaximum()), 0.8, -2.0);
  total.SetParLimits(0, 0.0, 1e7);
  total.SetParLimits(1, 0.85, 1.15);
  total.SetParLimits(2, 0.0, 1e9);
  total.SetParLimits(3, 0.0, 10.0);
  total.SetParLimits(4, -100.0, 20.0);
  hData->Fit(&total, "RQ0");

  TH1D* hDisp = static_cast<TH1D*>(hData->Clone((std::string(hData->GetName()) + "_disp").c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double displayScale = hDisp->GetXaxis()->GetBinWidth(1) / hData->GetXaxis()->GetBinWidth(1);

  TF1 signalDraw(("fSignalData_" + category).c_str(), SignalOnlyShape, gFitMin, gFitMax, 2);
  signalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  signalDraw.SetParameter(1, total.GetParameter(1));
  signalDraw.SetLineColor(kBlue + 1);
  signalDraw.SetLineWidth(3);

  TF1 backgroundDraw(("fBackgroundData_" + category).c_str(), BackgroundShape, gFitMin, gFitMax, 3);
  backgroundDraw.SetParameters(total.GetParameter(2) * displayScale, total.GetParameter(3), total.GetParameter(4));
  backgroundDraw.SetLineColor(kGreen + 2);
  backgroundDraw.SetLineWidth(3);
  backgroundDraw.SetLineStyle(2);

  TF1 totalDraw(total);
  totalDraw.SetParameter(0, total.GetParameter(0) * displayScale);
  totalDraw.SetParameter(2, total.GetParameter(2) * displayScale);
  totalDraw.SetLineColor(kRed + 1);
  totalDraw.SetLineWidth(4);

  TCanvas c(("cData_" + category).c_str(), ("cData_" + category).c_str(), 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle((category + " data reco-only S+B fit").c_str());
  hDisp->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Pairs / bin");
  hDisp->Draw("E");
  totalDraw.Draw("same");
  signalDraw.Draw("same");
  backgroundDraw.Draw("same");

  TLegend leg(0.56, 0.64, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "Reco-only same-event data", "lep");
  leg.AddEntry(&totalDraw, "Signal + background", "l");
  leg.AddEntry(&signalDraw,
               signalModel == SignalModel::TripleGaussian ? "Triple-Gaussian signal" : "Nominal signal model", "l");
  leg.AddEntry(&backgroundDraw, "Threshold background", "l");
  leg.Draw();

  TPaveText txt(0.12, 0.60, 0.52, 0.88, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("Signal mean fixed at %.5f GeV", gSignalShape.mean));
  txt.AddText(Form("Signal width scale = %.3f", total.GetParameter(1)));
  txt.AddText(Form("Threshold p = %.3f, b_{1} = %.3f", total.GetParameter(3), total.GetParameter(4)));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull(("hPullData_" + category).c_str(), "", hDisp->GetNbinsX(),
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
    if (xMax < gFitMin || xMin > gFitMax) continue;
    const double yfit = total.Integral(xMin, xMax) / hData->GetXaxis()->GetBinWidth(1);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs((outputDir + "/phi_data_sb_" + category + "_fit.pdf").c_str());
  delete hDisp;

  FitSummary summary;
  summary.category = category;
  summary.chi2 = total.GetChisquare();
  summary.ndf = total.GetNDF();
  summary.chi2ndf = total.GetChisquare() / total.GetNDF();
  summary.signalAmp = total.GetParameter(0);
  summary.signalWidthScale = total.GetParameter(1);
  summary.thresholdNorm = total.GetParameter(2);
  summary.thresholdPower = total.GetParameter(3);
  summary.thresholdSlope = total.GetParameter(4);
  return summary;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput = getArgument(argc, argv, "--signal-input", "PhiSignalOnlyHistograms.root");
  const std::string dataInput = getArgument(argc, argv, "--data-input", "PhiSBHistogramsData.root");
  const std::string outputDir = getArgument(argc, argv, "--output-dir", "SBDataFitResults");
  const SignalModel signalModel = parseSignalModel(
      getArgument(argc, argv, "--signal-model", "GaussPlusRightTailCB"));
  gFitMin = getDoubleArgument(argc, argv, "--fit-min", gFitMin);
  gFitMax = getDoubleArgument(argc, argv, "--fit-max", gFitMax);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile signalFile(signalInput.c_str(), "READ");
  TFile dataFile(dataInput.c_str(), "READ");
  if (signalFile.IsZombie() || dataFile.IsZombie()) {
    std::cerr << "Error opening input ROOT files" << std::endl;
    return 1;
  }

  TH1D* hSignal1 = nullptr;
  TH1D* hSignal2 = nullptr;
  TH1D* hData1 = nullptr;
  TH1D* hData2 = nullptr;
  signalFile.GetObject("hPhiMass1Tag", hSignal1);
  signalFile.GetObject("hPhiMass2Tag", hSignal2);
  dataFile.GetObject("hPhiSBMass1Tag", hData1);
  dataFile.GetObject("hPhiSBMass2Tag", hData2);
  if (hSignal1 == nullptr || hSignal2 == nullptr || hData1 == nullptr || hData2 == nullptr) {
    std::cerr << "Missing required histograms" << std::endl;
    return 1;
  }

  FitSummary s1 = fitCategory(hSignal1, hData1, "1tag", outputDir, signalModel);
  FitSummary s2 = fitCategory(hSignal2, hData2, "2tag", outputDir, signalModel);

  std::ofstream out(outputDir + "/phi_data_sb_fit_summary.csv");
  out << "category,signalModel,fitMin,fitMax,chi2,ndf,chi2ndf,signalAmp,signalWidthScale,thresholdNorm,thresholdPower,thresholdSlope\n";
  for (const FitSummary& s : {s1, s2})
    out << s.category << ","
        << (signalModel == SignalModel::TripleGaussian ? "TripleGaussian" : "GaussPlusRightTailCB") << ","
        << gFitMin << "," << gFitMax << ","
        << s.chi2 << "," << s.ndf << "," << s.chi2ndf << ","
        << s.signalAmp << "," << s.signalWidthScale << "," << s.thresholdNorm << "," << s.thresholdPower << ","
        << s.thresholdSlope << "\n";
  out.close();

  std::cout << "1tag chi2/ndf = " << s1.chi2ndf << std::endl;
  std::cout << "2tag chi2/ndf = " << s2.chi2ndf << std::endl;
  return 0;
}
