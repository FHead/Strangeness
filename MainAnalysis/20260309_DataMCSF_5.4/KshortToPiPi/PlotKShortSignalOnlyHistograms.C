void PlotKShortSignalOnlyHistograms(const char* inputFileName = "KShortSignalOnlyHistograms.root",
                                    const char* outputFileName = "Plots/kshort_signal_only_counts.pdf")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile inputFile(inputFileName, "READ");
  TH1D* hAccepted = nullptr;
  TH1D* h1Tag = nullptr;
  TH1D* h2Tag = nullptr;
  inputFile.GetObject("hKShortMassAccepted", hAccepted);
  inputFile.GetObject("hKShortMass1Tag", h1Tag);
  inputFile.GetObject("hKShortMass2Tag", h2Tag);
  if (hAccepted == nullptr || h1Tag == nullptr || h2Tag == nullptr) {
    std::cerr << "Missing required histograms in " << inputFileName << std::endl;
    return;
  }

  gSystem->mkdir("Plots", kTRUE);

  TH1D* hAcceptedDraw = (TH1D*)hAccepted->Clone("hAcceptedDraw");
  TH1D* h1TagDraw = (TH1D*)h1Tag->Clone("h1TagDraw");
  TH1D* h2TagDraw = (TH1D*)h2Tag->Clone("h2TagDraw");

  hAcceptedDraw->SetLineColor(kBlack);
  hAcceptedDraw->SetLineWidth(2);
  h1TagDraw->SetLineColor(kBlue + 1);
  h1TagDraw->SetLineWidth(2);
  h2TagDraw->SetLineColor(kRed + 1);
  h2TagDraw->SetLineWidth(2);

  const double maxY = std::max(hAcceptedDraw->GetMaximum(),
                      std::max(h1TagDraw->GetMaximum(), h2TagDraw->GetMaximum()));

  TCanvas canvas("canvas", "canvas", 900, 700);
  hAcceptedDraw->SetTitle("K^{0}_{S} #rightarrow #pi^{+}#pi^{-} signal-only MC");
  hAcceptedDraw->GetYaxis()->SetTitle("Candidates / bin");
  hAcceptedDraw->SetMaximum(maxY * 1.15);
  hAcceptedDraw->Draw("hist");
  h1TagDraw->Draw("hist same");
  h2TagDraw->Draw("hist same");

  TLegend legend(0.60, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(hAcceptedDraw, "Accepted", "l");
  legend.AddEntry(h1TagDraw, "1-tag", "l");
  legend.AddEntry(h2TagDraw, "2-tag", "l");
  legend.Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.15, 0.86, Form("Accepted = %.0f", hAcceptedDraw->Integral()));
  latex.DrawLatex(0.15, 0.81, Form("1-tag = %.0f", h1TagDraw->Integral()));
  latex.DrawLatex(0.15, 0.76, Form("2-tag = %.0f", h2TagDraw->Integral()));

  canvas.SaveAs(outputFileName);
}
