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
constexpr int kDisplayRebin = 4;
std::string gSignalModel = "Voigt";

struct SignalShape {
  double mean = 1.865;
  double sigma = 0.008;
  double gamma = 0.015;
  double alphaL = 1.5;
  double nL = 5.0;
  double alphaR = 1.5;
  double nR = 5.0;
};

struct FitSummary {
  std::string category;
  std::string model;
  double chi2 = 0.0;
  double ndf = 0.0;
  double chi2ndf = 1e9;
  double signalYield = 0.0;
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

double SignalOnlyShape(double* x, double* p) {
  if (gSignalModel == "DoubleSidedCB") {
    const double xx = x[0];
    const double norm = p[0];
    const double mean = gSignalShape.mean;
    const double sigma = gSignalShape.sigma;
    const double alphaL = gSignalShape.alphaL;
    const double nL = gSignalShape.nL;
    const double alphaR = gSignalShape.alphaR;
    const double nR = gSignalShape.nR;
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
  return p[0] * TMath::Voigt(x[0] - gSignalShape.mean, gSignalShape.sigma, gSignalShape.gamma, 4);
}

double Exp1(double* x, double* p) {
  return p[0] * std::exp(p[1] * x[0]);
}

double Exp2(double* x, double* p) {
  return p[0] * std::exp(p[1] * x[0] + p[2] * x[0] * x[0]);
}

double TotalExp1(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[2] = {p[1], p[2]};
  return SignalOnlyShape(x, sp) + Exp1(x, bp);
}

double TotalExp2(double* x, double* p) {
  double sp[1] = {p[0]};
  double bp[3] = {p[1], p[2], p[3]};
  return SignalOnlyShape(x, sp) + Exp2(x, bp);
}

SignalShape deriveSignalShape(TH1D* h) {
  SignalShape shape;
  if (gSignalModel == "DoubleSidedCB") {
    TF1 f("fSignalShapeD0DSCB",
          [](double* x, double* p) {
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
          },
          gFitMin, gFitMax, 7);
    f.SetParameters(std::max(20.0, 0.8 * h->GetMaximum()), 1.865, 0.018, 1.5, 5.0, 1.5, 5.0);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.002, 0.08);
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
  } else {
    TF1 f("fSignalShapeD0Voigt", "[0]*TMath::Voigt(x-[1],[2],[3],4)", gFitMin, gFitMax);
    f.SetParameters(std::max(20.0, 0.8 * h->GetMaximum()), 1.865, 0.008, 0.015);
    f.SetParLimits(1, 1.82, 1.90);
    f.SetParLimits(2, 0.0005, 0.05);
    f.SetParLimits(3, 0.001, 0.08);
    h->Fit(&f, "RQ0");
    shape.mean = f.GetParameter(1);
    shape.sigma = f.GetParameter(2);
    shape.gamma = f.GetParameter(3);
  }
  return shape;
}

TF1 buildTotalModel(const std::string& model, const std::string& name, TH1D* hSB) {
  const double maxY = std::max(50.0, hSB->GetMaximum());
  TF1 f;
  if (model == "Exp1") {
    f = TF1(name.c_str(), TotalExp1, gFitMin, gFitMax, 3);
    f.SetParNames("S", "N", "b1");
    f.SetParameters(0.10 * maxY, 0.8 * maxY, -2.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e12);
    f.SetParLimits(2, -20.0, 5.0);
  } else {
    f = TF1(name.c_str(), TotalExp2, gFitMin, gFitMax, 4);
    f.SetParNames("S", "N", "b1", "b2");
    f.SetParameters(0.10 * maxY, 0.8 * maxY, -2.0, 0.0);
    f.SetParLimits(0, 0.0, 1e9);
    f.SetParLimits(1, 0.0, 1e12);
    f.SetParLimits(2, -20.0, 5.0);
    f.SetParLimits(3, -10.0, 10.0);
  }
  return f;
}

void scaleDisplayedFunction(TF1& f, const std::string& model, double scale) {
  f.SetParameter(0, f.GetParameter(0) * scale);
  f.SetParameter(1, f.GetParameter(1) * scale);
}

TF1 buildBackgroundDraw(const std::string& model, const std::string& name, TF1& total, double scale) {
  TF1 f;
  if (model == "Exp1") {
    f = TF1(name.c_str(), Exp1, gFitMin, gFitMax, 2);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2));
  } else {
    f = TF1(name.c_str(), Exp2, gFitMin, gFitMax, 3);
    f.SetParameters(total.GetParameter(1) * scale, total.GetParameter(2), total.GetParameter(3));
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
  hDisp->SetTitle((category + " D^{0} LooseID reco-only S+B fit").c_str());
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
  leg.AddEntry(&signalDraw, "Fixed Voigt signal", "l");
  leg.AddEntry(&backgroundDraw, "Background component", "l");
  leg.Draw();

  TPaveText txt(0.14, 0.14, 0.44, 0.30, "NDC");
  txt.SetBorderSize(0);
  txt.SetFillStyle(0);
  txt.AddText(Form("Model: %s", model.c_str()));
  txt.AddText(Form("#chi^{2}/ndf = %.3f", total.GetChisquare() / total.GetNDF()));
  txt.AddText(Form("Signal: %s", gSignalModel.c_str()));
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

  c.SaveAs((outputDir + "/d0_looseid_sb_" + category + "_" + model + "_fit.pdf").c_str());
  delete hDisp;

  FitSummary summary;
  summary.category = category;
  summary.model = model;
  summary.chi2 = total.GetChisquare();
  summary.ndf = total.GetNDF();
  summary.chi2ndf = total.GetChisquare() / total.GetNDF();
  summary.signalYield = signalDraw.Integral(gFitMin, gFitMax) / hSB->GetXaxis()->GetBinWidth(1);
  return summary;
}

FitSummary fitCategory(TH1D* hSignal, TH1D* hSB, const std::string& category,
                       const std::string& outputDir, std::ofstream& csv) {
  gSignalShape = deriveSignalShape(hSignal);

  const std::vector<std::string> models = {"Exp1", "Exp2"};
  std::vector<FitSummary> results;
  for (const std::string& model : models)
    results.push_back(runFit(hSB, category, model, outputDir));

  const FitSummary& best = *std::min_element(results.begin(), results.end(),
                                             [](const FitSummary& a, const FitSummary& b) {
                                               return a.chi2ndf < b.chi2ndf;
                                             });

  for (const FitSummary& r : results) {
    csv << r.category << "," << r.model << "," << r.chi2 << "," << r.ndf << ","
        << r.chi2ndf << "," << r.signalYield << ","
        << (r.model == best.model ? "yes" : "no") << "\n";
  }
  return best;
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string signalInput = getArgument(argc, argv, "--signal-input", "d0_looseid_combined_assignment_histograms.root");
  const std::string sbInput = getArgument(argc, argv, "--sb-input", "d0_looseid_sb_histograms.root");
  const std::string outputDir = getArgument(argc, argv, "--output-dir", "SBFitResults");
  gSignalModel = getArgument(argc, argv, "--signal-model", gSignalModel);
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
  signalFile.GetObject("hD0MassKaonTag", hSignal1);
  signalFile.GetObject("hD0MassKaonPionTag", hSignal2);
  sbFile.GetObject("hD0SBMassKaonTag", hSB1);
  sbFile.GetObject("hD0SBMassKaonPionTag", hSB2);
  if (hSignal1 == nullptr || hSignal2 == nullptr || hSB1 == nullptr || hSB2 == nullptr) {
    std::cerr << "Error: required histograms not found" << std::endl;
    return 1;
  }

  std::ofstream out(outputDir + "/d0_sb_fit_summary.csv");
  out << "category,model,chi2,ndf,chi2ndf,signalYield,isBest\n";

  FitSummary s1 = fitCategory(hSignal1, hSB1, "kaon_tag", outputDir, out);
  FitSummary s2 = fitCategory(hSignal2, hSB2, "kaon_pion_tag", outputDir, out);
  out.close();

  std::cout << "kaon_tag best: " << s1.model << " chi2/ndf = " << s1.chi2ndf << std::endl;
  std::cout << "kaon_pion_tag best: " << s2.model << " chi2/ndf = " << s2.chi2ndf << std::endl;
  std::cout << "Wrote outputs to " << outputDir << std::endl;
  return 0;
}
