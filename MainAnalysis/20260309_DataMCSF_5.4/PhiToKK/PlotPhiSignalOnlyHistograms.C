#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void setStyle(TH1D* h, int color) {
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetLineWidth(2);
}

int PlotPhiSignalOnlyHistograms(std::string input = "PhiSignalOnlyHistograms.root",
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
  TH1D* h1Tag = nullptr;
  TH1D* h2Tag = nullptr;
  file.GetObject("hPhiMassAccepted", hAccepted);
  file.GetObject("hPhiMass1Tag", h1Tag);
  file.GetObject("hPhiMass2Tag", h2Tag);
  if (hAccepted == nullptr || h1Tag == nullptr || h2Tag == nullptr) {
    std::cerr << "Missing histograms in " << input << std::endl;
    return 1;
  }

  TH1D* hAcceptedDraw = static_cast<TH1D*>(hAccepted->Clone("hAcceptedDraw"));
  TH1D* h1TagDraw = static_cast<TH1D*>(h1Tag->Clone("h1TagDraw"));
  TH1D* h2TagDraw = static_cast<TH1D*>(h2Tag->Clone("h2TagDraw"));
  setStyle(hAcceptedDraw, kGray + 2);
  setStyle(h1TagDraw, kBlue + 1);
  setStyle(h2TagDraw, kRed + 1);

  TCanvas c1("c1", "c1", 800, 550);
  hAcceptedDraw->SetTitle("");
  hAcceptedDraw->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  hAcceptedDraw->GetYaxis()->SetTitle("Candidates / bin");
  hAcceptedDraw->SetMaximum(hAcceptedDraw->GetMaximum() * 1.25);
  hAcceptedDraw->Draw("hist");
  h1TagDraw->Draw("hist same");
  h2TagDraw->Draw("hist same");

  TLegend leg1(0.62, 0.67, 0.88, 0.88);
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.AddEntry(hAcceptedDraw, "Accepted", "l");
  leg1.AddEntry(h1TagDraw, "1-tag", "l");
  leg1.AddEntry(h2TagDraw, "2-tag", "l");
  leg1.Draw();

  TLatex text1;
  text1.SetNDC();
  text1.SetTextSize(0.032);
  text1.DrawLatex(0.15, 0.88, Form("Accepted = %.0f", hAcceptedDraw->Integral()));
  text1.DrawLatex(0.15, 0.83, Form("1-tag = %.0f", h1TagDraw->Integral()));
  text1.DrawLatex(0.15, 0.78, Form("2-tag = %.0f", h2TagDraw->Integral()));
  c1.SaveAs((outputDir + "/phi_signal_only_counts.pdf").c_str());

  TH1D* h1Shape = static_cast<TH1D*>(h1Tag->Clone("h1Shape"));
  TH1D* h2Shape = static_cast<TH1D*>(h2Tag->Clone("h2Shape"));
  if (h1Shape->Integral() > 0) h1Shape->Scale(1.0 / h1Shape->Integral());
  if (h2Shape->Integral() > 0) h2Shape->Scale(1.0 / h2Shape->Integral());
  setStyle(h1Shape, kBlue + 1);
  setStyle(h2Shape, kRed + 1);

  TCanvas c2("c2", "c2", 800, 550);
  h1Shape->SetTitle("");
  h1Shape->GetXaxis()->SetTitle("m(K^{+}K^{-}) [GeV]");
  h1Shape->GetYaxis()->SetTitle("Unit-normalized candidates / bin");
  h1Shape->SetMaximum(std::max(h1Shape->GetMaximum(), h2Shape->GetMaximum()) * 1.25);
  h1Shape->Draw("hist");
  h2Shape->Draw("hist same");

  TLegend leg2(0.68, 0.76, 0.88, 0.88);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.AddEntry(h1Shape, "1-tag", "l");
  leg2.AddEntry(h2Shape, "2-tag", "l");
  leg2.Draw();
  c2.SaveAs((outputDir + "/phi_signal_only_shapes.pdf").c_str());

  delete hAcceptedDraw;
  delete h1TagDraw;
  delete h2TagDraw;
  delete h1Shape;
  delete h2Shape;
  file.Close();
  return 0;
}
