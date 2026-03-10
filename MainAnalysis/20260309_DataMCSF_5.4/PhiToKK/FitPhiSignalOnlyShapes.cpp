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
double gFitMin = 1.000;
double gFitMax = 1.050;
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

double GaussPlusRightTailCB(double* x, double* p) {
  const double xx = x[0];
  const double gNorm = p[0];
  const double mean = p[1];
  const double sigmaG = p[2];
  const double cbNorm = p[3];
  const double sigmaCB = p[4];
  const double alpha = p[5];
  const double n = p[6];

  if (sigmaG <= 0.0 || sigmaCB <= 0.0) return 0.0;

  const double tg = (xx - mean) / sigmaG;
  const double gauss = gNorm * std::exp(-0.5 * tg * tg);

  double cb = 0.0;
  if (alpha > 0.0 && n > 1.0) {
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
}

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

TF1 buildGaussian(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(), "gaus(0)", gFitMin, gFitMax);
  f.SetParNames("A", "mean", "sigma");
  f.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0030);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0005, 0.02);
  return f;
}

TF1 buildVoigt(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3],4)", gFitMin, gFitMax);
  f.SetParNames("A", "mean", "sigma", "gamma");
  f.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0015, 0.0035);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0002, 0.01);
  f.SetParLimits(3, 0.0002, 0.02);
  return f;
}

TF1 buildDoubleGauss(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
        gFitMin, gFitMax);
  f.SetParNames("A1", "mean", "sigma1", "A2", "sigma2");
  f.SetParameters(std::max(60.0, 0.6 * h->GetMaximum()), 1.0195, 0.0020,
                  std::max(30.0, 0.3 * h->GetMaximum()), 0.0045);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(4, 0.0003, 0.03);
  return f;
}

TF1 buildTripleGauss(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)+"
        "[5]*exp(-0.5*((x-[1])/[6])^2)",
        gFitMin, gFitMax);
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

TF1 buildRightTailCB(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(), CrystalBallRightTail, gFitMin, gFitMax, 5);
  f.SetParNames("A", "mean", "sigma", "alpha", "n");
  f.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0028, 1.5, 6.0);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(3, 0.2, 8.0);
  f.SetParLimits(4, 1.2, 80.0);
  return f;
}

TF1 buildDoubleSidedCB(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(), DoubleSidedCrystalBall, gFitMin, gFitMax, 7);
  f.SetParNames("A", "mean", "sigma", "alphaL", "nL", "alphaR", "nR");
  f.SetParameters(std::max(100.0, h->GetMaximum()), 1.0195, 0.0025, 1.2, 5.0, 1.8, 8.0);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(3, 0.2, 8.0);
  f.SetParLimits(4, 1.2, 80.0);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  return f;
}

TF1 buildGaussPlusRightTailCB(const std::string& name, TH1D* h) {
  TF1 f(name.c_str(), GaussPlusRightTailCB, gFitMin, gFitMax, 7);
  f.SetParNames("AG", "mean", "sigmaG", "ACB", "sigmaCB", "alpha", "n");
  f.SetParameters(std::max(60.0, 0.6 * h->GetMaximum()), 1.0195, 0.0016,
                  std::max(30.0, 0.3 * h->GetMaximum()), 0.0035, 1.6, 8.0);
  f.SetParLimits(1, 1.015, 1.024);
  f.SetParLimits(2, 0.0003, 0.02);
  f.SetParLimits(4, 0.0003, 0.03);
  f.SetParLimits(5, 0.2, 8.0);
  f.SetParLimits(6, 1.2, 80.0);
  return f;
}

ModelResult runFit(TH1D* h, TF1& f, const std::string& name) {
  ModelResult result;
  result.name = name;
  TH1D* hc = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_" + name + "_tmp").c_str()));
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
  if (f.GetNpar() > 3) {
    result.width2 = f.GetParameter(3);
    result.width2Err = f.GetParError(3);
  }
  if (f.GetNpar() > 4) {
    result.width2 = f.GetParameter(4);
    result.width2Err = f.GetParError(4);
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
  if (originalWidth <= 0.0) return 1.0;
  return displayWidth / originalWidth;
}

double evalModel(const std::string& model, TF1& gauss, TF1& voigt, TF1& dgauss, TF1& tgauss,
                 TF1& rcb, TF1& dscb, TF1& gcb, double x) {
  if (model == "Gaussian") return gauss.Eval(x);
  if (model == "Voigt") return voigt.Eval(x);
  if (model == "DoubleGaussian") return dgauss.Eval(x);
  if (model == "TripleGaussian") return tgauss.Eval(x);
  if (model == "RightTailCB") return rcb.Eval(x);
  if (model == "DoubleSidedCB") return dscb.Eval(x);
  if (model == "GaussPlusRightTailCB") return gcb.Eval(x);
  return 0.0;
}

double evalModelBinYield(const std::string& model, TF1& gauss, TF1& voigt, TF1& dgauss, TF1& tgauss,
                         TF1& rcb, TF1& dscb, TF1& gcb,
                         double xMin, double xMax, double originalBinWidth) {
  TF1* f = nullptr;
  if (model == "Gaussian") f = &gauss;
  if (model == "Voigt") f = &voigt;
  if (model == "DoubleGaussian") f = &dgauss;
  if (model == "TripleGaussian") f = &tgauss;
  if (model == "RightTailCB") f = &rcb;
  if (model == "DoubleSidedCB") f = &dscb;
  if (model == "GaussPlusRightTailCB") f = &gcb;
  if (f == nullptr) return 0.0;
  if (originalBinWidth <= 0.0) return 0.0;
  return f->Integral(xMin, xMax) / originalBinWidth;
}

void scaleDisplayedFunction(TF1& f, const std::string& model, double displayScale) {
  if (model == "Gaussian" || model == "Voigt" || model == "RightTailCB" || model == "DoubleSidedCB")
    f.SetParameter(0, f.GetParameter(0) * displayScale);
  else if (model == "DoubleGaussian") {
    f.SetParameter(0, f.GetParameter(0) * displayScale);
    f.SetParameter(3, f.GetParameter(3) * displayScale);
  } else if (model == "TripleGaussian") {
    f.SetParameter(0, f.GetParameter(0) * displayScale);
    f.SetParameter(3, f.GetParameter(3) * displayScale);
    f.SetParameter(5, f.GetParameter(5) * displayScale);
  } else if (model == "GaussPlusRightTailCB") {
    f.SetParameter(0, f.GetParameter(0) * displayScale);
    f.SetParameter(3, f.GetParameter(3) * displayScale);
  }
}

void drawFitPage(TH1D* h, TF1& gauss, TF1& voigt, TF1& dgauss, TF1& tgauss,
                 TF1& rcb, TF1& dscb, TF1& gcb,
                 const std::vector<ModelResult>& results, const ModelResult& best,
                 const std::string& title, const std::string& outputPdf) {
  TH1D* hDisp = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + "_disp").c_str()));
  if (kDisplayRebin > 1) hDisp->Rebin(kDisplayRebin);
  const double displayScale = getScaleFactor(h, hDisp);
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
  hDisp->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hDisp->GetYaxis()->SetTitle("Candidates / bin");
  hDisp->Draw("E");

  TF1 gaussDraw(gauss);
  TF1 voigtDraw(voigt);
  TF1 dgaussDraw(dgauss);
  TF1 tgaussDraw(tgauss);
  TF1 rcbDraw(rcb);
  TF1 dscbDraw(dscb);
  TF1 gcbDraw(gcb);

  scaleDisplayedFunction(gaussDraw, "Gaussian", displayScale);
  scaleDisplayedFunction(voigtDraw, "Voigt", displayScale);
  scaleDisplayedFunction(dgaussDraw, "DoubleGaussian", displayScale);
  scaleDisplayedFunction(tgaussDraw, "TripleGaussian", displayScale);
  scaleDisplayedFunction(rcbDraw, "RightTailCB", displayScale);
  scaleDisplayedFunction(dscbDraw, "DoubleSidedCB", displayScale);
  scaleDisplayedFunction(gcbDraw, "GaussPlusRightTailCB", displayScale);

  gaussDraw.SetLineColor(kBlue + 1);
  gaussDraw.SetLineWidth(2);
  gaussDraw.Draw("same");
  voigtDraw.SetLineColor(kRed + 1);
  voigtDraw.SetLineWidth(2);
  voigtDraw.Draw("same");
  dgaussDraw.SetLineColor(kGreen + 2);
  dgaussDraw.SetLineWidth(2);
  dgaussDraw.Draw("same");
  tgaussDraw.SetLineColor(kMagenta + 1);
  tgaussDraw.SetLineWidth(2);
  tgaussDraw.Draw("same");
  rcbDraw.SetLineColor(kOrange + 7);
  rcbDraw.SetLineWidth(2);
  rcbDraw.Draw("same");
  dscbDraw.SetLineColor(kCyan + 2);
  dscbDraw.SetLineWidth(2);
  dscbDraw.Draw("same");
  gcbDraw.SetLineColor(kBlack);
  gcbDraw.SetLineWidth(3);
  gcbDraw.Draw("same");

  TLegend leg(0.58, 0.58, 0.89, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hDisp, "MC histogram", "lep");
  leg.AddEntry(&gaussDraw, "Gaussian", "l");
  leg.AddEntry(&voigtDraw, "Voigt", "l");
  leg.AddEntry(&dgaussDraw, "Double Gaussian", "l");
  leg.AddEntry(&tgaussDraw, "Triple Gaussian", "l");
  leg.AddEntry(&rcbDraw, "Right-tail CB", "l");
  leg.AddEntry(&dscbDraw, "Double-sided CB", "l");
  leg.AddEntry(&gcbDraw, "Gauss + right-tail CB", "l");
  leg.Draw();

  TPaveText txt(0.13, 0.56, 0.48, 0.88, "NDC");
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
    const double x = hDisp->GetBinCenter(i);
    const double xMin = hDisp->GetBinLowEdge(i);
    const double xMax = xMin + hDisp->GetBinWidth(i);
    if (xMax < gFitMin || xMin > gFitMax) continue;
    const double y = hDisp->GetBinContent(i);
    const double ey = hDisp->GetBinError(i);
    if (ey <= 0.0) continue;
    const double yfit =
        evalModelBinYield(best.name, gauss, voigt, dgauss, tgauss, rcb, dscb, gcb, xMin, xMax, originalBinWidth);
    hPull.SetBinContent(i, (y - yfit) / ey);
    hPull.SetBinError(i, 1.0);
  }
  hPull.Draw("E");
  TLine line0(gFitMin, 0.0, gFitMax, 0.0);
  line0.SetLineStyle(2);
  line0.Draw("same");

  c.SaveAs(outputPdf.c_str());
  delete hDisp;
}

void writeSummary(const std::string& path, const std::string& category,
                  const std::vector<ModelResult>& results, const ModelResult& best) {
  std::ofstream out(path, std::ios::app);
  out << "[" << category << "]\n";
  out << "BestModel=" << best.name << "\n";
  for (const ModelResult& r : results) {
    out << r.name << ",chi2=" << r.chi2
        << ",ndf=" << r.ndf
        << ",chi2ndf=" << r.chi2ndf
        << ",mean=" << r.mean
        << ",meanErr=" << r.meanErr
        << ",width1=" << r.width1
        << ",width1Err=" << r.width1Err
        << ",width2=" << r.width2
        << ",width2Err=" << r.width2Err
        << ",width3=" << r.width3
        << ",width3Err=" << r.width3Err
        << "\n";
  }
  out << "\n";
}

void fitCategory(TH1D* h, const std::string& category, const std::string& outputDir, std::ofstream& csv) {
  TF1 gauss = buildGaussian("fGauss_" + category, h);
  TF1 voigt = buildVoigt("fVoigt_" + category, h);
  TF1 dgauss = buildDoubleGauss("fDouble_" + category, h);
  TF1 tgauss = buildTripleGauss("fTriple_" + category, h);
  TF1 rcb = buildRightTailCB("fRightTailCB_" + category, h);
  TF1 dscb = buildDoubleSidedCB("fDoubleSidedCB_" + category, h);
  TF1 gcb = buildGaussPlusRightTailCB("fGaussPlusRightTailCB_" + category, h);

  std::vector<ModelResult> results;
  results.push_back(runFit(h, gauss, "Gaussian"));
  results.push_back(runFit(h, voigt, "Voigt"));
  results.push_back(runFit(h, dgauss, "DoubleGaussian"));
  results.push_back(runFit(h, tgauss, "TripleGaussian"));
  results.push_back(runFit(h, rcb, "RightTailCB"));
  results.push_back(runFit(h, dscb, "DoubleSidedCB"));
  results.push_back(runFit(h, gcb, "GaussPlusRightTailCB"));

  const ModelResult& best = *std::min_element(
      results.begin(), results.end(),
      [](const ModelResult& a, const ModelResult& b) { return a.chi2ndf < b.chi2ndf; });

  drawFitPage(h, gauss, voigt, dgauss, tgauss, rcb, dscb, gcb, results, best,
              category + " signal-only MC fits", outputDir + "/phi_" + category + "_fit_scan.pdf");

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
  const std::string inputFileName = getArgument(argc, argv, "--input", "PhiSignalOnlyHistograms.root");
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
  inputFile.GetObject("hPhiMass1Tag", h1);
  inputFile.GetObject("hPhiMass2Tag", h2);
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
  inputFile.Close();
  std::cout << "Wrote fit outputs to " << outputDir << std::endl;
  return 0;
}
