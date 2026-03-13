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

int PlotD0SBHistograms(std::string input = "d0_looseid_sb_histograms.root",
                       std::string outputDir = "PlotsSB") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open " << input << std::endl;
    return 1;
  }

  TH1D* hAccepted = nullptr;
  TH1D* hKaonTag = nullptr;
  TH1D* hKaonPionTag = nullptr;
  file.GetObject("hD0SBMassAccepted", hAccepted);
  file.GetObject("hD0SBMassKaonTag", hKaonTag);
  file.GetObject("hD0SBMassKaonPionTag", hKaonPionTag);
  if (hAccepted == nullptr || hKaonTag == nullptr || hKaonPionTag == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  TH1D* hAcceptedDraw = static_cast<TH1D*>(hAccepted->Clone("hD0SBAcceptedDraw"));
  TH1D* hKaonTagDraw = static_cast<TH1D*>(hKaonTag->Clone("hD0SBKaonTagDraw"));
  TH1D* hKaonPionTagDraw = static_cast<TH1D*>(hKaonPionTag->Clone("hD0SBKaonPionTagDraw"));
  SetStyle(hAcceptedDraw, kGray + 2);
  SetStyle(hKaonTagDraw, kBlue + 1);
  SetStyle(hKaonPionTagDraw, kRed + 1);

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

  TLatex latex;
  latex.SetNDC();

  p1.cd();
  hAcceptedDraw->SetTitle("");
  hAcceptedDraw->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hAcceptedDraw->GetYaxis()->SetTitle("Assignments / bin");
  hAcceptedDraw->GetXaxis()->SetLabelSize(0.0);
  hAcceptedDraw->SetMinimum(0.0);
  hAcceptedDraw->SetMaximum(hAcceptedDraw->GetMaximum() * 1.20);
  hAcceptedDraw->Draw("hist");

  TLegend leg1(0.63, 0.80, 0.88, 0.89);
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.AddEntry(hAcceptedDraw, "All assignment fills", "l");
  leg1.Draw();

  latex.SetTextSize(0.032);
  latex.DrawLatex(0.15, 0.88, Form("Accepted fills = %.0f", hAcceptedDraw->Integral()));

  p2.cd();
  hKaonTagDraw->SetTitle("");
  hKaonTagDraw->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hKaonTagDraw->GetYaxis()->SetTitle("Assignments / bin");
  hKaonTagDraw->SetMinimum(0.0);
  hKaonTagDraw->SetMaximum(std::max(hKaonTagDraw->GetMaximum(), hKaonPionTagDraw->GetMaximum()) * 1.20);
  hKaonTagDraw->Draw("hist");
  hKaonPionTagDraw->Draw("hist same");

  TLegend leg2(0.63, 0.77, 0.88, 0.90);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.AddEntry(hKaonTagDraw, "Kaon tag only", "l");
  leg2.AddEntry(hKaonPionTagDraw, "Kaon + pion tags", "l");
  leg2.Draw();

  latex.SetTextSize(0.030);
  latex.DrawLatex(0.15, 0.88, Form("Kaon tag = %.0f", hKaonTagDraw->Integral()));
  latex.DrawLatex(0.15, 0.82, Form("Kaon + pion = %.0f", hKaonPionTagDraw->Integral()));

  c.SaveAs((outputDir + "/d0_looseid_sb_counts.pdf").c_str());

  delete hAcceptedDraw;
  delete hKaonTagDraw;
  delete hKaonPionTagDraw;
  file.Close();
  return 0;
}
