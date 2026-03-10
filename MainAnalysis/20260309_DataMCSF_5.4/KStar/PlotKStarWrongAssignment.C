#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void StyleHistogram(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

void DrawComparison(TH1D* hNominal, TH1D* hSwapped, const std::string& title,
                    const std::string& outputName) {
  TCanvas c("c", "c", 800, 600);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);

  TH1D* hNom = static_cast<TH1D*>(hNominal->Clone((std::string(hNominal->GetName()) + "_draw").c_str()));
  TH1D* hSwap = static_cast<TH1D*>(hSwapped->Clone((std::string(hSwapped->GetName()) + "_draw").c_str()));
  StyleHistogram(hNom, kBlue + 1);
  StyleHistogram(hSwap, kRed + 1);

  hNom->SetTitle(title.c_str());
  hNom->GetXaxis()->SetTitle("m(K#pi) [GeV]");
  hNom->GetYaxis()->SetTitle("Candidates / bin");
  hNom->SetMinimum(0.0);
  hNom->SetMaximum(std::max(hNom->GetMaximum(), hSwap->GetMaximum()) * 1.15);
  hNom->Draw("hist");
  hSwap->Draw("hist same");

  TLegend leg(0.56, 0.72, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hNom, "Nominal assignment", "l");
  leg.AddEntry(hSwap, "Swapped assignment", "l");
  leg.Draw();

  c.SaveAs(outputName.c_str());
  delete hNom;
  delete hSwap;
}

int PlotKStarWrongAssignment(std::string input = "KStarWrongAssignmentHistograms.root",
                             std::string outputDir = "PlotsWrongAssignment") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDir.c_str(), true);

  TFile file(input.c_str(), "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open " << input << std::endl;
    return 1;
  }

  TH1D* hNominalNoPionTag = nullptr;
  TH1D* hSwappedNoPionTag = nullptr;
  TH1D* hNominalWithPionTag = nullptr;
  TH1D* hSwappedWithPionTag = nullptr;
  file.GetObject("hNominalNoPionTag", hNominalNoPionTag);
  file.GetObject("hSwappedNoPionTag", hSwappedNoPionTag);
  file.GetObject("hNominalWithPionTag", hNominalWithPionTag);
  file.GetObject("hSwappedWithPionTag", hSwappedWithPionTag);
  if (hNominalNoPionTag == nullptr || hSwappedNoPionTag == nullptr ||
      hNominalWithPionTag == nullptr || hSwappedWithPionTag == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  DrawComparison(hNominalNoPionTag, hSwappedNoPionTag,
                 "K^{*} wrong-assignment study, no pion tag",
                 outputDir + "/kstar_wrong_assignment_no_pion_tag.pdf");
  DrawComparison(hNominalWithPionTag, hSwappedWithPionTag,
                 "K^{*} wrong-assignment study, with pion tag",
                 outputDir + "/kstar_wrong_assignment_with_pion_tag.pdf");
  return 0;
}
