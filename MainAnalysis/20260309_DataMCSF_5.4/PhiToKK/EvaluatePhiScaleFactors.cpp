#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kThreshold = 2.0 * kKaonMass;

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

struct FitYieldResult {
  double yield1Tag = 0.0;
  double yield2Tag = 0.0;
};

struct SummaryRow {
  std::string label;
  double mcYield1 = 0.0;
  double mcYield2 = 0.0;
  double dataYield1 = 0.0;
  double dataYield2 = 0.0;
  double mcEfficiency = 0.0;
  double dataEfficiency = 0.0;
  double scaleFactor = 0.0;
};

double RightTailCBUnit(double x, double mean, double sigma, double alpha, double n) {
  if (sigma <= 0.0 || alpha <= 0.0 || n <= 1.0) return 0.0;
  const double t = (x - mean) / sigma;
  if (t <= alpha) return std::exp(-0.5 * t * t);
  const double A = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double B = n / alpha - alpha;
  return A / std::pow(B + t, n);
}

SignalShape deriveSignalShape(TH1D* h, SignalModel model, double fitMin, double fitMax) {
  if (model == SignalModel::TripleGaussian) {
    TF1 f("fSignalTriple",
          "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+[5]*exp(-0.5*((x-[1])/[6])^2)",
          fitMin, fitMax);
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

  TF1 f("fSignalNominal", [](double* x, double* p) {
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
  }, fitMin, fitMax, 7);

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

double signalUnit(double x, const SignalShape& shape, double widthScale) {
  const double sigma1 = shape.sigma1 * widthScale;
  const double sigma2 = shape.sigma2 * widthScale;
  const double sigma3 = shape.sigma3 * widthScale;
  if (shape.model == SignalModel::TripleGaussian) {
    const double g1 = std::exp(-0.5 * std::pow((x - shape.mean) / sigma1, 2));
    const double g2 = std::exp(-0.5 * std::pow((x - shape.mean) / sigma2, 2));
    const double g3 = std::exp(-0.5 * std::pow((x - shape.mean) / sigma3, 2));
    return g1 + shape.scale2 * g2 + shape.scale3 * g3;
  }
  const double g = std::exp(-0.5 * std::pow((x - shape.mean) / sigma1, 2));
  const double cb = RightTailCBUnit(x, shape.mean, sigma2, shape.alpha, shape.n);
  return g + shape.scale2 * cb;
}

FitYieldResult fitMC(TH1D* hSignal1, TH1D* hSignal2, TH1D* hMC1, TH1D* hMC2,
                     SignalModel model, double fitMin, double fitMax) {
  FitYieldResult result;
  const double binWidth1 = hMC1->GetXaxis()->GetBinWidth(1);
  const double binWidth2 = hMC2->GetXaxis()->GetBinWidth(1);
  const SignalShape shape1 = deriveSignalShape(hSignal1, model, fitMin, fitMax);
  const SignalShape shape2 = deriveSignalShape(hSignal2, model, fitMin, fitMax);

  auto doFit = [&](TH1D* h, const SignalShape& shape, const std::string& name, double binWidth) {
    TF1 total(name.c_str(), [&](double* x, double* p) {
      const double xx = x[0];
      const double signal = p[0] * signalUnit(xx, shape, 1.0);
      const double gauss1 = p[1] * std::exp(-0.5 * std::pow((xx - p[2]) / p[3], 2));
      const double gauss2 = p[4] * std::exp(-0.5 * std::pow((xx - p[5]) / p[6], 2));
      double threshold = 0.0;
      if (xx > kThreshold)
        threshold = p[7] * std::pow(xx - kThreshold, p[8]) * std::exp(p[9] * xx);
      return signal + gauss1 + gauss2 + threshold;
    }, fitMin, fitMax, 10);
    total.SetParameters(std::max(100.0, 0.4 * h->GetMaximum()),
                        std::max(100.0, 0.08 * h->GetMaximum()), 1.000, 0.003,
                        std::max(60.0, 0.05 * h->GetMaximum()), 1.004, 0.006,
                        std::max(100.0, 0.2 * h->GetMaximum()), 0.8, -2.0);
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
    h->Fit(&total, "RQ0");

    TF1 signalOnly("signalOnly", [&](double* x, double* p) {
      return p[0] * signalUnit(x[0], shape, 1.0);
    }, fitMin, fitMax, 1);
    signalOnly.SetParameter(0, total.GetParameter(0));
    return signalOnly.Integral(fitMin, fitMax) / binWidth;
  };

  result.yield1Tag = doFit(hMC1, shape1, "fMCTotal1", binWidth1);
  result.yield2Tag = doFit(hMC2, shape2, "fMCTotal2", binWidth2);
  return result;
}

FitYieldResult fitData(TH1D* hSignal1, TH1D* hSignal2, TH1D* hData1, TH1D* hData2,
                       SignalModel model, double fitMin, double fitMax) {
  FitYieldResult result;
  const double binWidth1 = hData1->GetXaxis()->GetBinWidth(1);
  const double binWidth2 = hData2->GetXaxis()->GetBinWidth(1);
  const SignalShape shape1 = deriveSignalShape(hSignal1, model, fitMin, fitMax);
  const SignalShape shape2 = deriveSignalShape(hSignal2, model, fitMin, fitMax);

  auto doFit = [&](TH1D* h, const SignalShape& shape, const std::string& name, double binWidth) {
    TF1 total(name.c_str(), [&](double* x, double* p) {
      const double xx = x[0];
      const double signal = p[0] * signalUnit(xx, shape, p[1]);
      double threshold = 0.0;
      if (xx > kThreshold)
        threshold = p[2] * std::pow(xx - kThreshold, p[3]) * std::exp(p[4] * xx);
      return signal + threshold;
    }, fitMin, fitMax, 5);
    total.SetParameters(std::max(100.0, 0.3 * h->GetMaximum()), 1.0,
                        std::max(100.0, 0.2 * h->GetMaximum()), 0.8, -2.0);
    total.SetParLimits(0, 0.0, 1e7);
    total.SetParLimits(1, 0.85, 1.15);
    total.SetParLimits(2, 0.0, 1e9);
    total.SetParLimits(3, 0.0, 10.0);
    total.SetParLimits(4, -100.0, 20.0);
    h->Fit(&total, "RQ0");

    TF1 signalOnly("signalOnlyData", [&](double* x, double* p) {
      return p[0] * signalUnit(x[0], shape, p[1]);
    }, fitMin, fitMax, 2);
    signalOnly.SetParameters(total.GetParameter(0), total.GetParameter(1));
    return signalOnly.Integral(fitMin, fitMax) / binWidth;
  };

  result.yield1Tag = doFit(hData1, shape1, "fDataTotal1", binWidth1);
  result.yield2Tag = doFit(hData2, shape2, "fDataTotal2", binWidth2);
  return result;
}

double efficiency(double n1, double n2) {
  const double den = n1 + 2.0 * n2;
  return (den > 0.0) ? 2.0 * n2 / den : 0.0;
}

double efficiencyError(double n1, double n2) {
  const double den = n1 + 2.0 * n2;
  if (den <= 0.0) return 0.0;
  const double d1 = -2.0 * n2 / (den * den);
  const double d2 = 2.0 / den - 4.0 * n2 / (den * den);
  return std::sqrt(d1 * d1 * n1 + d2 * d2 * n2);
}

double scaleFactorError(const SummaryRow& row) {
  const double dataErr = efficiencyError(row.dataYield1, row.dataYield2);
  const double mcErr = efficiencyError(row.mcYield1, row.mcYield2);
  if (row.dataEfficiency <= 0.0 || row.mcEfficiency <= 0.0) return 0.0;
  return row.scaleFactor * std::sqrt(std::pow(dataErr / row.dataEfficiency, 2) +
                                     std::pow(mcErr / row.mcEfficiency, 2));
}

SummaryRow evaluate(const std::string& label,
                    const std::string& signalFileName,
                    const std::string& mcFileName,
                    const std::string& dataFileName,
                    SignalModel model,
                    double fitMin, double fitMax) {
  TFile signalFile(signalFileName.c_str(), "READ");
  TFile mcFile(mcFileName.c_str(), "READ");
  TFile dataFile(dataFileName.c_str(), "READ");

  TH1D *hSignal1 = nullptr, *hSignal2 = nullptr, *hMC1 = nullptr, *hMC2 = nullptr, *hData1 = nullptr, *hData2 = nullptr;
  signalFile.GetObject("hPhiMass1Tag", hSignal1);
  signalFile.GetObject("hPhiMass2Tag", hSignal2);
  mcFile.GetObject("hPhiSBMass1Tag", hMC1);
  mcFile.GetObject("hPhiSBMass2Tag", hMC2);
  dataFile.GetObject("hPhiSBMass1Tag", hData1);
  dataFile.GetObject("hPhiSBMass2Tag", hData2);

  SummaryRow row;
  row.label = label;
  const FitYieldResult mc = fitMC(hSignal1, hSignal2, hMC1, hMC2, model, fitMin, fitMax);
  const FitYieldResult data = fitData(hSignal1, hSignal2, hData1, hData2, model, fitMin, fitMax);
  row.mcYield1 = mc.yield1Tag;
  row.mcYield2 = mc.yield2Tag;
  row.dataYield1 = data.yield1Tag;
  row.dataYield2 = data.yield2Tag;
  row.mcEfficiency = efficiency(row.mcYield1, row.mcYield2);
  row.dataEfficiency = efficiency(row.dataYield1, row.dataYield2);
  row.scaleFactor = (row.mcEfficiency > 0.0) ? row.dataEfficiency / row.mcEfficiency : 0.0;
  return row;
}
}  // namespace

int main() {
  gROOT->SetBatch(kTRUE);

  const std::string nominalSignal = "PhiSignalOnlyHistograms.root";
  const std::string nominalMC = "PhiSBHistograms.root";
  const std::string nominalData = "PhiSBHistogramsData.root";
  const std::string matchSignal = "Systematics/MatchingAngle/Match0025/PhiSignalOnlyHistograms.root";

  std::vector<SummaryRow> rows;
  rows.push_back(evaluate("Nominal", nominalSignal, nominalMC, nominalData,
                          SignalModel::GaussPlusRightTailCB, 0.99, 1.06));
  rows.push_back(evaluate("SignalFunctionTripleGaussian", nominalSignal, nominalMC, nominalData,
                          SignalModel::TripleGaussian, 0.99, 1.06));
  rows.push_back(evaluate("FitRange100to105", nominalSignal, nominalMC, nominalData,
                          SignalModel::GaussPlusRightTailCB, 1.00, 1.05));
  rows.push_back(evaluate("MatchingAngle0025", matchSignal, nominalMC, nominalData,
                          SignalModel::GaussPlusRightTailCB, 0.99, 1.06));

  std::ofstream out("PhiScaleFactorSummary.csv");
  out << "label,mcYield1,mcYield2,dataYield1,dataYield2,mcEfficiency,dataEfficiency,scaleFactor,statError\n";
  for (const SummaryRow& row : rows) {
    out << row.label << ","
        << row.mcYield1 << "," << row.mcYield2 << ","
        << row.dataYield1 << "," << row.dataYield2 << ","
        << row.mcEfficiency << "," << row.dataEfficiency << ","
        << row.scaleFactor << "," << scaleFactorError(row) << "\n";
  }
  out.close();

  for (const SummaryRow& row : rows) {
    std::cout << row.label
              << " SF=" << row.scaleFactor
              << " stat=" << scaleFactorError(row)
              << " dataEff=" << row.dataEfficiency
              << " mcEff=" << row.mcEfficiency << std::endl;
  }
  return 0;
}
