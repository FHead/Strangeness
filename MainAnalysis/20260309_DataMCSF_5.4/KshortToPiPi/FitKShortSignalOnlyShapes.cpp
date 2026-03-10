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
double gFitMin = 0.30;
double gFitMax = 1.00;
constexpr int kDisplayRebin = 4;

struct ModelResult {
  std::string name;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 1e9;
  double mean = 0.0;
  double meanErr = 0.0;
  double width1 = 0.0;
  double width1Err = 0.0;
  double width2 = 0.0;
  double width2Err = 0.0;
  double width3 = 0.0;
  double width3Err = 0.0;
};

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

double CrystalBallLeftTail(double* x, double* p) {
  const double xx = x[0];
  const double norm = p[0];
  const double mean = p[1];
  const double sigma = p[2];
  const double alpha = p[3];
  const double n = p[4];
  if (sigma <= 0.0 || alpha <= 0.0 || n <= 1.0) return 0.0;

  const double t = (xx - mean) / sigma;
  if (t >= -alpha)
    return norm * std::exp(-0.5 * t * t);

  const double A = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double B = n / alpha - alpha;
  return norm * A / std::pow(B - t, n);
}

double CrystalBallRightTail(double* x, double* p) {
  const double xx = x[0];
  const double norm = p[0];
  const double mean = p[1];
  const double sigma = p[2];
  const double alpha = p[3];
  const double n = p[4];
  if (sigma <= 0.0 || alpha <= 0.0 || n <= 1.0) return 0.0;

  const double t = (xx - mean) / sigma;
  if (t <= alpha)
    return norm * std::exp(-0.5 * t * t);

  const double A = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
  const double B = n / alpha - alpha;
  return norm * A / std::pow(B + t, n);
}

double DoubleSidedCrystalBall(double* x, double* p) {
  const double xx = x[0];
  const double norm = p[0];
  const double mean = p[1];
  const double sigma = p[2];
  const double alphaL = p[3];
  const double nL = p[4];
  const double alphaR = p[5];
  const double nR = p[6];
  if (sigma <= 0.0 || alphaL <= 0.0 || nL <= 1.0 || alphaR <= 0.0 || nR <= 1.0) return 0.0;

  const double t = (xx - mean) / sigma;
  if (t > -alphaL && t < alphaR)
    return norm * std::exp(-0.5 * t * t);

  if (t <= -alphaL) {
    const double A = std::pow(nL / alphaL, nL) * std::exp(-0.5 * alphaL * alphaL);
    const double B = nL / alphaL - alphaL;
    return norm * A / std::pow(B - t, nL);
  }

  const double A = std::pow(nR / alphaR, nR) * std::exp(-0.5 * alphaR * alphaR);
  const double B = nR / alphaR - alphaR;
  return norm * A / std::pow(B + t, nR);
}

double GaussPlusLeftCB(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double cbp[5] = {p[3], p[1], p[4], p[5], p[6]};
  return g + CrystalBallLeftTail(x, cbp);
}

double GaussPlusRightCB(double* x, double* p) {
  const double xx = x[0];
  const double g = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  double cbp[5] = {p[3], p[1], p[4], p[5], p[6]};
  return g + CrystalBallRightTail(x, cbp);
}

double DoubleGaussPlusLeftCB(double* x, double* p) {
  const double xx = x[0];
  const double g1 = p[0] * std::exp(-0.5 * std::pow((xx - p[1]) / p[2], 2));
  const double g2 = p[3] * std::exp(-0.5 * std::pow((xx - p[1]) / p[4], 2));
  double cbp[5] = {p[5], p[1], p[6], p[7], p[8]};
  return g1 + g2 + CrystalBallLeftTail(x, cbp);
}

double DoubleSidedCBPlusGauss(double* x, double* p) {
  double dscbp[7] = {p[0], p[1], p[2], p[3], p[4], p[5], p[6]};
  const double cb = DoubleSidedCrystalBall(x, dscbp);
  const double xx = x[0];
  const double g = p[7] * std::exp(-0.5 * std::pow((xx - p[1]) / p[8], 2));
  return cb + g;
}

TF1 buildModel(const std::string& model, const std::string& name, TH1D* h) {
  const double maxY = std::max(100.0, h->GetMaximum());
  TF1 f;
  if (model == "DoubleGaussian") {
    f = TF1(name.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
            gFitMin, gFitMax);
    f.SetParameters(0.65 * maxY, 0.4976, 0.004, 0.35 * maxY, 0.012);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
  } else if (model == "TripleGaussian") {
    f = TF1(name.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)",
            gFitMin, gFitMax);
    f.SetParameters(0.5 * maxY, 0.4976, 0.003, 0.3 * maxY, 0.009, 0.2 * maxY, 0.03);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
  } else if (model == "QuadGaussian") {
    f = TF1(name.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
            "[5]*exp(-0.5*((x-[1])/[6])^2)+[7]*exp(-0.5*((x-[1])/[8])^2)",
            gFitMin, gFitMax);
    f.SetParameters(0.42 * maxY, 0.4976, 0.003, 0.28 * maxY, 0.008,
                    0.20 * maxY, 0.02, 0.10 * maxY, 0.05);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.04);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
    f.SetParLimits(8, 0.001, 0.20);
  } else if (model == "BifurGauss") {
    f = TF1(name.c_str(),
            "[0]*(x<[1]?exp(-0.5*((x-[1])/[2])^2):exp(-0.5*((x-[1])/[3])^2))",
            gFitMin, gFitMax);
    f.SetParameters(maxY, 0.4976, 0.004, 0.012);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.001, 0.08);
  } else if (model == "LeftTailCB") {
    f = TF1(name.c_str(), CrystalBallLeftTail, gFitMin, gFitMax, 5);
    f.SetParameters(maxY, 0.4976, 0.008, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
  } else if (model == "RightTailCB") {
    f = TF1(name.c_str(), CrystalBallRightTail, gFitMin, gFitMax, 5);
    f.SetParameters(maxY, 0.4976, 0.008, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
  } else if (model == "DoubleSidedCB") {
    f = TF1(name.c_str(), DoubleSidedCrystalBall, gFitMin, gFitMax, 7);
    f.SetParameters(maxY, 0.4976, 0.008, 1.5, 5.0, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
  } else if (model == "GaussPlusLeftCB") {
    f = TF1(name.c_str(), GaussPlusLeftCB, gFitMin, gFitMax, 7);
    f.SetParameters(0.7 * maxY, 0.4976, 0.006, 0.3 * maxY, 0.015, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
  } else if (model == "GaussPlusRightCB") {
    f = TF1(name.c_str(), GaussPlusRightCB, gFitMin, gFitMax, 7);
    f.SetParameters(0.7 * maxY, 0.4976, 0.006, 0.3 * maxY, 0.015, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
  } else if (model == "DoubleGaussPlusLeftCB") {
    f = TF1(name.c_str(), DoubleGaussPlusLeftCB, gFitMin, gFitMax, 9);
    f.SetParameters(0.45 * maxY, 0.4976, 0.004, 0.30 * maxY, 0.010, 0.25 * maxY, 0.020, 1.5, 5.0);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(4, 0.001, 0.08);
    f.SetParLimits(6, 0.001, 0.12);
    f.SetParLimits(7, 0.2, 8.0);
    f.SetParLimits(8, 1.2, 80.0);
  } else if (model == "DoubleSidedCBPlusGauss") {
    f = TF1(name.c_str(), DoubleSidedCBPlusGauss, gFitMin, gFitMax, 9);
    f.SetParameters(0.8 * maxY, 0.4976, 0.008, 1.5, 5.0, 1.5, 5.0, 0.2 * maxY, 0.03);
    f.SetParLimits(1, 0.485, 0.510);
    f.SetParLimits(2, 0.001, 0.05);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
    f.SetParLimits(8, 0.001, 0.15);
  }
  return f;
}

ModelResult runFit(TH1D* h, TF1& f, const std::string& model) {
  ModelResult result;
  result.name = model;
  TH1D* hc = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_" + model + "_tmp").c_str()));
  hc->Fit(&f, "RQ0");
  result.chi2 = f.GetChisquare();
  result.ndf = f.GetNDF();
  if (result.ndf > 0) result.chi2ndf = result.chi2 / result.ndf;
  if (f.GetNpar() > 1) {
    result.mean = f.GetParameter(1);
    result.meanErr = f.GetParError(1);
  }
  if (f.GetNpar() > 2) {
    result.width1 = f.GetParameter(2);
    result.width1Err = f.GetParError(2);
  }
  if (f.GetNpar() > 4) {
    result.width2 = f.GetParameter(4);
    result.width2Err = f.GetParError(4);
  } else if (f.GetNpar() > 3) {
    result.width2 = f.GetParameter(3);
    result.width2Err = f.GetParError(3);
  }
  if (f.GetNpar() > 6) {
    result.width3 = f.GetParameter(6);
    result.width3Err = f.GetParError(6);
  }
  delete hc;
  return result;
}

double getScaleFactor(TH1D* original, TH1D* display) {
  const double originalWidth = original->GetXaxis()->GetBinWidth(1);
  const double displayWidth = display->GetXaxis()->GetBinWidth(1);
  return (originalWidth > 0.0) ? displayWidth / originalWidth : 1.0;
}

void scaleDisplayedFunction(TF1& f, const std::string& model, double scale) {
  if (model == "BifurGauss" || model == "LeftTailCB" || model == "RightTailCB")
    f.SetParameter(0, f.GetParameter(0) * scale);
  else if (model == "DoubleGaussian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
  } else if (model == "TripleGaussian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
  } else if (model == "QuadGaussian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
    f.SetParameter(7, f.GetParameter(7) * scale);
  } else if (model == "DoubleSidedCB") {
    f.SetParameter(0, f.GetParameter(0) * scale);
  } else if (model == "DoubleSidedCBPlusGauss") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(7, f.GetParameter(7) * scale);
  } else if (model == "GaussPlusLeftCB" || model == "GaussPlusRightCB") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
  } else if (model == "DoubleGaussPlusLeftCB") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
    f.SetParameter(5, f.GetParameter(5) * scale);
  }
}

double evalModelBinYield(TF1& f, double xMin, double xMax, double originalBinWidth) {
  return (originalBinWidth > 0.0) ? f.Integral(xMin, xMax) / originalBinWidth : 0.0;
}

void drawFitPage(TH1D* h, const std::vector<std::pair<std::string, TF1>>& models,
                 const std::vector<ModelResult>& results, const ModelResult& best,
                 const std::string& title, const std::string& outputPdf) {
  TH1D* hDisp = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_disp").c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double scale = getScaleFactor(h, hDisp);
  const double originalBinWidth = h->GetXaxis()->GetBinWidth(1);

  TCanvas c("c", "c", 900, 800);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetPad(0.0, 0.32, 1.0, 1.0);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.03);
  hDisp->SetStats(0);
  hDisp->SetLineWidth(2);
  hDisp->SetTitle(title.c_str());
  hDisp->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Candidates / bin");
  hDisp->SetMinimum(0.0);
  hDisp->Draw("E");

  std::vector<int> colors = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kBlack};
  TLegend leg(0.58, 0.56, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "MC histogram", "lep");
  std::vector<TF1*> drawFunctions;
  for (size_t i = 0; i < models.size(); ++i) {
    TF1* fd = new TF1(models[i].second);
    fd->SetName(Form("fdraw_%s_%zu", h->GetName(), i));
    fd->SetTitle(models[i].first.c_str());
    fd->SetNpx(1000);
    scaleDisplayedFunction(*fd, models[i].first, scale);
    fd->SetLineColor(colors[i % colors.size()]);
    fd->SetLineWidth(models[i].first == best.name ? 4 : 2);
    drawFunctions.push_back(fd);
  }
  for (size_t i = 0; i < drawFunctions.size(); ++i) {
    drawFunctions[i]->Draw("same");
    leg.AddEntry(drawFunctions[i], models[i].first.c_str(), "l");
  }
  leg.Draw();

  TPaveText txt(0.58, 0.12, 0.89, 0.44, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  for (const ModelResult& r : results)
    txt.AddText(Form("%s: #chi^{2}/ndf = %.3f", r.name.c_str(), r.chi2ndf));
  txt.AddText(Form("Best: %s", best.name.c_str()));
  txt.Draw();

  c.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.32);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.34);
  TH1D hPull("hPull", "", hDisp->GetNbinsX(), hDisp->GetXaxis()->GetXmin(), hDisp->GetXaxis()->GetXmax());
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

  const TF1* bestFunction = nullptr;
  for (const auto& item : models)
    if (item.first == best.name) bestFunction = &item.second;

  for (int i = 1; i <= hDisp->GetNbinsX(); ++i) {
    const double xMin = hDisp->GetBinLowEdge(i);
    const double xMax = xMin + hDisp->GetBinWidth(i);
    if (xMax < gFitMin || xMin > gFitMax || bestFunction == nullptr) continue;
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0.0) continue;
    const double yfit = evalModelBinYield(*const_cast<TF1*>(bestFunction), xMin, xMax, originalBinWidth);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs(outputPdf.c_str());
  for (TF1* f : drawFunctions)
    delete f;
  delete hDisp;
}

void fitCategory(TH1D* h, const std::string& category, const std::string& outputDir, std::ofstream& csv) {
  const std::vector<std::string> modelNames = {
      "TripleGaussian", "QuadGaussian", "DoubleSidedCB", "DoubleSidedCBPlusGauss"};
  std::vector<std::pair<std::string, TF1>> models;
  std::vector<ModelResult> results;
  for (const std::string& name : modelNames) {
    TF1 f = buildModel(name, "f_" + category + "_" + name, h);
    results.push_back(runFit(h, f, name));
    models.push_back({name, f});
  }

  const ModelResult& best = *std::min_element(
      results.begin(), results.end(),
      [](const ModelResult& a, const ModelResult& b) { return a.chi2ndf < b.chi2ndf; });

  drawFitPage(h, models, results, best,
              category + " K^{0}_{S} signal-only MC fits",
              outputDir + "/kshort_" + category + "_fit_scan.pdf");

  for (const ModelResult& r : results) {
    csv << category << "," << r.name << ","
        << r.chi2 << "," << r.ndf << "," << r.chi2ndf << ","
        << r.mean << "," << r.meanErr << ","
        << r.width1 << "," << r.width1Err << ","
        << r.width2 << "," << r.width2Err << ","
        << r.width3 << "," << r.width3Err << ","
        << (r.name == best.name ? "yes" : "no") << "\n";
  }

  std::cout << category << " best model: " << best.name
            << " with chi2/ndf = " << best.chi2ndf
            << ", mean = " << best.mean
            << ", width1 = " << best.width1 << std::endl;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName = getArgument(argc, argv, "--input", "KShortSignalOnlyHistograms.root");
  const std::string outputDir = getArgument(argc, argv, "--output-dir", "FitResults");
  gFitMin = getDoubleArgument(argc, argv, "--fit-min", gFitMin);
  gFitMax = getDoubleArgument(argc, argv, "--fit-max", gFitMax);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
    return 1;
  }

  TH1D* h1 = nullptr;
  TH1D* h2 = nullptr;
  inputFile.GetObject("hKShortMass1Tag", h1);
  inputFile.GetObject("hKShortMass2Tag", h2);
  if (h1 == nullptr || h2 == nullptr) {
    std::cerr << "Error: required histograms not found in " << inputFileName << std::endl;
    return 1;
  }

  gSystem->mkdir(outputDir.c_str(), true);
  std::ofstream csv(outputDir + "/fit_summary.csv");
  csv << "category,model,chi2,ndf,chi2ndf,mean,meanErr,width1,width1Err,width2,width2Err,width3,width3Err,isBest\n";

  fitCategory(h1, "1tag", outputDir, csv);
  fitCategory(h2, "2tag", outputDir, csv);

  csv.close();
  std::cout << "Wrote fit outputs to " << outputDir << std::endl;
  return 0;
}
