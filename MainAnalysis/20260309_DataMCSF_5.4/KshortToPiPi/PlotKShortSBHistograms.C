#include <algorithm>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void SetStyle(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

int PlotKShortSBHistograms(std::string input = "KShortSBHistograms.root",
                           std::string outputDir = "Plots") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open " << input << std::endl;
    return 1;
  }

  TH1D* hAccepted = nullptr;
  TH1D* h1 = nullptr;
  TH1D* h2 = nullptr;
  file.GetObject("hKShortSBMassAccepted", hAccepted);
  file.GetObject("hKShortSBMass1Tag", h1);
  file.GetObject("hKShortSBMass2Tag", h2);
  if (hAccepted == nullptr || h1 == nullptr || h2 == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  TH1D* hAcceptedDraw = static_cast<TH1D*>(hAccepted->Clone("hKShortSBAcceptedDraw"));
  TH1D* h1Draw = static_cast<TH1D*>(h1->Clone("hKShortSB1Draw"));
  TH1D* h2Draw = static_cast<TH1D*>(h2->Clone("hKShortSB2Draw"));
  SetStyle(hAcceptedDraw, kGray + 2);
  SetStyle(h1Draw, kBlue + 1);
  SetStyle(h2Draw, kRed + 1);

  TCanvas c("c", "c", 800, 800);
  TPad p1("p1", "p1", 0.0, 0.42, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.42);
  p1.SetLeftMargin(0.14);
  p1.SetRightMargin(0.05);
  p1.SetBottomMargin(0.03);
  p1.SetLogy(0);
  p2.SetLeftMargin(0.14);
  p2.SetRightMargin(0.05);
  p2.SetTopMargin(0.03);
  p2.SetBottomMargin(0.16);
  p2.SetLogy(0);
  p1.Draw();
  p2.Draw();

  p1.cd();
  hAcceptedDraw->SetTitle("");
  hAcceptedDraw->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  hAcceptedDraw->GetYaxis()->SetTitle("Pairs / bin");
  hAcceptedDraw->GetXaxis()->SetLabelSize(0.0);
  hAcceptedDraw->SetMaximum(hAcceptedDraw->GetMaximum() * 1.20);
  hAcceptedDraw->SetMinimum(0.0);
  hAcceptedDraw->Draw("hist");

  TLegend leg1(0.66, 0.79, 0.88, 0.88);
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.AddEntry(hAcceptedDraw, "Accepted OS pairs", "l");
  leg1.Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.15, 0.88, Form("Accepted = %.0f", hAcceptedDraw->Integral()));

  p2.cd();
  TH1D* h1Linear = static_cast<TH1D*>(h1->Clone("hKShortSB1Linear"));
  TH1D* h2Linear = static_cast<TH1D*>(h2->Clone("hKShortSB2Linear"));
  SetStyle(h1Linear, kBlue + 1);
  SetStyle(h2Linear, kRed + 1);
  h1Linear->SetTitle("");
  h1Linear->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
  h1Linear->GetYaxis()->SetTitle("Pairs / bin");
  h1Linear->SetMinimum(0.0);
  h1Linear->SetMaximum(std::max(h1Linear->GetMaximum(), h2Linear->GetMaximum()) * 1.20);
  h1Linear->Draw("hist");
  h2Linear->Draw("hist same");

  TLegend leg2(0.68, 0.77, 0.88, 0.90);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.AddEntry(h1Linear, "1-tag", "l");
  leg2.AddEntry(h2Linear, "2-tag", "l");
  leg2.Draw();

  latex.SetTextSize(0.030);
  latex.DrawLatex(0.15, 0.88, Form("1-tag = %.0f", h1Linear->Integral()));
  latex.DrawLatex(0.15, 0.82, Form("2-tag = %.0f", h2Linear->Integral()));

  c.SaveAs((outputDir + "/kshort_sb_counts.pdf").c_str());

  delete hAcceptedDraw;
  delete h1Draw;
  delete h2Draw;
  delete h1Linear;
  delete h2Linear;
  file.Close();
  return 0;
}
