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
double gFitMin = 1.70;
double gFitMax = 2.00;
constexpr int kDisplayRebin = 2;

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

double DoubleSidedCBPlusGauss(double* x, double* p) {
  double dscbp[7] = {p[0], p[1], p[2], p[3], p[4], p[5], p[6]};
  const double cb = DoubleSidedCrystalBall(x, dscbp);
  const double xx = x[0];
  const double g = p[7] * std::exp(-0.5 * std::pow((xx - p[1]) / p[8], 2));
  return cb + g;
}

TF1 buildModel(const std::string& model, const std::string& name, TH1D* h) {
  const double maxY = std::max(50.0, h->GetMaximum());
  TF1 f;
  if (model == "Gaussian") {
    f = TF1(name.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2)", gFitMin, gFitMax);
    f.SetParameters(maxY, 1.865, 0.02);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.002, 0.08);
  } else if (model == "Lorentzian") {
    f = TF1(name.c_str(), "[0] / (1.0 + ((x-[1])/[2])^2)", gFitMin, gFitMax);
    f.SetParameters(maxY, 1.865, 0.02);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.001, 0.08);
  } else if (model == "Voigt") {
    f = TF1(name.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)", gFitMin, gFitMax);
    f.SetParameters(maxY, 1.865, 0.006, 0.02);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.0005, 0.05);
    f.SetParLimits(3, 0.001, 0.08);
  } else if (model == "DoubleGaussian") {
    f = TF1(name.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
            gFitMin, gFitMax);
    f.SetParameters(0.7 * maxY, 1.865, 0.015, 0.3 * maxY, 0.035);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.002, 0.06);
    f.SetParLimits(4, 0.002, 0.10);
  } else if (model == "DoubleSidedCB") {
    f = TF1(name.c_str(), DoubleSidedCrystalBall, gFitMin, gFitMax, 7);
    f.SetParameters(maxY, 1.865, 0.018, 1.5, 5.0, 1.5, 5.0);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.002, 0.08);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
  } else {
    f = TF1(name.c_str(), DoubleSidedCBPlusGauss, gFitMin, gFitMax, 9);
    f.SetParameters(0.8 * maxY, 1.865, 0.018, 1.5, 5.0, 1.5, 5.0, 0.2 * maxY, 0.04);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.002, 0.08);
    f.SetParLimits(3, 0.2, 8.0);
    f.SetParLimits(4, 1.2, 80.0);
    f.SetParLimits(5, 0.2, 8.0);
    f.SetParLimits(6, 1.2, 80.0);
    f.SetParLimits(8, 0.002, 0.12);
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
  if (model == "Gaussian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
  } else if (model == "Lorentzian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
  } else if (model == "Voigt") {
    f.SetParameter(0, f.GetParameter(0) * scale);
  } else if (model == "DoubleGaussian") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(3, f.GetParameter(3) * scale);
  } else if (model == "DoubleSidedCB") {
    f.SetParameter(0, f.GetParameter(0) * scale);
  } else if (model == "DoubleSidedCBPlusGauss") {
    f.SetParameter(0, f.GetParameter(0) * scale);
    f.SetParameter(7, f.GetParameter(7) * scale);
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
  hDisp->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hDisp->GetYaxis()->SetTitle("Candidates / bin");
  hDisp->SetMinimum(0.0);
  hDisp->Draw("E");

  std::vector<int> colors = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
  TLegend leg(0.56, 0.58, 0.89, 0.88);
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

  TPaveText txt(0.14, 0.14, 0.44, 0.34, "NDC");
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
  for (TF1* f : drawFunctions) delete f;
  delete hDisp;
}

void fitCategory(TH1D* h, const std::string& category, const std::string& outputDir, std::ofstream& csv) {
  const std::vector<std::string> modelNames = {
      "Gaussian", "Lorentzian", "Voigt", "DoubleGaussian", "DoubleSidedCB", "DoubleSidedCBPlusGauss"};
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
              category + " D^{0} LooseID signal-only MC fits",
              outputDir + "/d0_looseid_" + category + "_fit_scan.pdf");

  for (const ModelResult& r : results) {
    csv << category << "," << r.name << ","
        << r.chi2 << "," << r.ndf << "," << r.chi2ndf << ","
        << r.mean << "," << r.meanErr << ","
        << r.width1 << "," << r.width1Err << ","
        << r.width2 << "," << r.width2Err << ","
        << (r.name == best.name ? "yes" : "no") << "\n";
  }

  std::cout << category << " best model: " << best.name
            << " with chi2/ndf = " << best.chi2ndf
            << ", mean = " << best.mean
            << ", width1 = " << best.width1 << std::endl;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName = getArgument(argc, argv, "--input", "d0_looseid_combined_assignment_histograms.root");
  const std::string outputDir = getArgument(argc, argv, "--output-dir", "FitResultsCombined");
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
  inputFile.GetObject("hD0MassKaonTag", h1);
  inputFile.GetObject("hD0MassKaonPionTag", h2);
  if (h1 == nullptr || h2 == nullptr) {
    std::cerr << "Error: required histograms not found in " << inputFileName << std::endl;
    return 1;
  }

  gSystem->mkdir(outputDir.c_str(), true);
  std::ofstream csv(outputDir + "/fit_summary.csv");
  csv << "category,model,chi2,ndf,chi2ndf,mean,meanErr,width1,width1Err,width2,width2Err,isBest\n";

  fitCategory(h1, "kaon_tag", outputDir, csv);
  fitCategory(h2, "kaon_pion_tag", outputDir, csv);

  csv.close();
  std::cout << "Wrote fit outputs to " << outputDir << std::endl;
  return 0;
}
