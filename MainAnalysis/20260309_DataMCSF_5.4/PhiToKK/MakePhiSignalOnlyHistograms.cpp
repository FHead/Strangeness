#include <cmath>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kPhiMassWindowMin = 0.99;
constexpr double kPhiMassWindowMax = 1.06;
constexpr int kPhiMassBins = 280;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMatchAngleMax = 0.01;
constexpr long long kKaonTagThreshold = 2;
constexpr int kMaxReco = 10000;
constexpr int kMaxPhi = 4096;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
};

double buildMass(const TrackKinematics& t1, const TrackKinematics& t2) {
  const double p1sq = t1.px * t1.px + t1.py * t1.py + t1.pz * t1.pz;
  const double p2sq = t2.px * t2.px + t2.py * t2.py + t2.pz * t2.pz;
  const double e1 = std::sqrt(p1sq + kKaonMass * kKaonMass);
  const double e2 = std::sqrt(p2sq + kKaonMass * kKaonMass);

  const double px = t1.px + t2.px;
  const double py = t1.py + t2.py;
  const double pz = t1.pz + t2.pz;
  const double e = e1 + e2;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

bool passAcceptance(const TrackKinematics& t) {
  const double p = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz);
  if (p <= 0.0) return false;
  const double absCosTheta = std::fabs(t.pz / p);
  return (absCosTheta >= kAbsCosMin && absCosTheta <= kAbsCosMax);
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

void printUsage(const char* argv0) {
  std::cout << "Usage: " << argv0
            << " [--input Samples/merged_mc_v2.3.root]"
            << " [--output PhiSignalOnlyHistograms.root]"
            << " [--tree Tree]"
            << " [--mass-min 0.99]"
            << " [--mass-max 1.06]"
            << " [--match-angle-max 0.01]" << std::endl;
}
}  // namespace

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--help") {
      printUsage(argv[0]);
      return 0;
    }
  }

  const std::string inputFileName =
      getArgument(argc, argv, "--input", "../../../../Samples/merged_mc_v2.3.root");
  const std::string outputFileName =
      getArgument(argc, argv, "--output", "PhiSignalOnlyHistograms.root");
  const std::string treeName = getArgument(argc, argv, "--tree", "Tree");
  const double massMin = getDoubleArgument(argc, argv, "--mass-min", kPhiMassWindowMin);
  const double massMax = getDoubleArgument(argc, argv, "--mass-max", kPhiMassWindowMax);
  const double matchAngleMax = getDoubleArgument(argc, argv, "--match-angle-max", kMatchAngleMax);

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
    return 1;
  }

  TTree* tree = nullptr;
  inputFile.GetObject(treeName.c_str(), tree);
  if (tree == nullptr) {
    std::cerr << "Error: cannot find tree '" << treeName << "' in " << inputFileName << std::endl;
    return 1;
  }

  long long nReco = 0;
  double recoPx[kMaxReco] = {0.0};
  double recoPy[kMaxReco] = {0.0};
  double recoPz[kMaxReco] = {0.0};
  double recoCharge[kMaxReco] = {0.0};
  long long recoPIDKaon[kMaxReco] = {0};
  long long recoGoodTrack[kMaxReco] = {0};

  long long nPhi = 0;
  long long phiReco1ID[kMaxPhi] = {0};
  long long phiReco2ID[kMaxPhi] = {0};
  double phiReco1Angle[kMaxPhi] = {0.0};
  double phiReco2Angle[kMaxPhi] = {0.0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  tree->SetBranchAddress("NPhi", &nPhi);
  tree->SetBranchAddress("PhiReco1ID[NPhi]", phiReco1ID);
  tree->SetBranchAddress("PhiReco2ID[NPhi]", phiReco2ID);
  tree->SetBranchAddress("PhiReco1Angle[NPhi]", phiReco1Angle);
  tree->SetBranchAddress("PhiReco2Angle[NPhi]", phiReco2Angle);

  TH1D hMass1Tag("hPhiMass1Tag",
                 "#phi signal-only MC, 1-tag; m(K^{+}K^{-}) [GeV]; Candidates / bin",
                 kPhiMassBins, massMin, massMax);
  TH1D hMass2Tag("hPhiMass2Tag",
                 "#phi signal-only MC, 2-tag; m(K^{+}K^{-}) [GeV]; Candidates / bin",
                 kPhiMassBins, massMin, massMax);
  TH1D hMassAccepted("hPhiMassAccepted",
                     "#phi signal-only MC, accepted; m(K^{+}K^{-}) [GeV]; Candidates / bin",
                     kPhiMassBins, massMin, massMax);

  long long totalPhi = 0;
  long long passValidRecoID = 0;
  long long passMatchAngle = 0;
  long long passGoodTrack = 0;
  long long passAcceptanceBoth = 0;
  long long passOppositeCharge = 0;
  long long count1Tag = 0;
  long long count2Tag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    for (long long iPhi = 0; iPhi < nPhi; ++iPhi) {
      totalPhi++;

      const long long reco1 = phiReco1ID[iPhi];
      const long long reco2 = phiReco2ID[iPhi];
      if (reco1 < 0 || reco2 < 0 || reco1 >= nReco || reco2 >= nReco) continue;
      passValidRecoID++;

      if (phiReco1Angle[iPhi] >= matchAngleMax || phiReco2Angle[iPhi] >= matchAngleMax)
        continue;
      passMatchAngle++;

      if (recoGoodTrack[reco1] != 1 || recoGoodTrack[reco2] != 1) continue;
      passGoodTrack++;

      const TrackKinematics track1{recoPx[reco1], recoPy[reco1], recoPz[reco1]};
      const TrackKinematics track2{recoPx[reco2], recoPy[reco2], recoPz[reco2]};

      if (!passAcceptance(track1) || !passAcceptance(track2)) continue;
      passAcceptanceBoth++;

      if (recoCharge[reco1] * recoCharge[reco2] >= 0) continue;
      passOppositeCharge++;

      const double mass = buildMass(track1, track2);
      hMassAccepted.Fill(mass);

      int nTagged = 0;
      if (recoPIDKaon[reco1] >= kKaonTagThreshold) nTagged++;
      if (recoPIDKaon[reco2] >= kKaonTagThreshold) nTagged++;

      if (nTagged == 1) {
        hMass1Tag.Fill(mass);
        count1Tag++;
      }
      if (nTagged == 2) {
        hMass2Tag.Fill(mass);
        count2Tag++;
      }
    }
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hMass1Tag.Write();
  hMass2Tag.Write();
  hMassAccepted.Write();

  TNamed selection("SelectionSummary",
                   Form("Signal-only phi->KK from MC: valid PhiReco IDs, PhiReco1Angle<%.4f, "
                        "PhiReco2Angle<%.4f, both RecoGoodTrack==1, both 0.15<=|cos(theta)|<=0.675, "
                        "opposite charge, kaon mass hypothesis on both daughters, tag if RecoPIDKaon>=2, "
                        "hist range %.3f-%.3f GeV",
                        matchAngleMax, matchAngleMax, massMin, massMax));
  selection.Write();

  TParameter<long long>("TotalPhiCandidates", totalPhi).Write();
  TParameter<long long>("PassValidRecoID", passValidRecoID).Write();
  TParameter<long long>("PassMatchAngle", passMatchAngle).Write();
  TParameter<long long>("PassGoodTrack", passGoodTrack).Write();
  TParameter<long long>("PassAcceptanceBoth", passAcceptanceBoth).Write();
  TParameter<long long>("PassOppositeCharge", passOppositeCharge).Write();
  TParameter<long long>("Count1Tag", count1Tag).Write();
  TParameter<long long>("Count2Tag", count2Tag).Write();
  TParameter<double>("MassMin", massMin).Write();
  TParameter<double>("MassMax", massMax).Write();
  TParameter<double>("MatchAngleMax", matchAngleMax).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Total phi candidates: " << totalPhi << std::endl;
  std::cout << "  Pass valid reco IDs:  " << passValidRecoID << std::endl;
  std::cout << "  Pass match angles:    " << passMatchAngle << std::endl;
  std::cout << "  Pass good tracks:     " << passGoodTrack << std::endl;
  std::cout << "  Pass acceptance:      " << passAcceptanceBoth << std::endl;
  std::cout << "  Pass opposite charge: " << passOppositeCharge << std::endl;
  std::cout << "  1-tag entries:        " << count1Tag << std::endl;
  std::cout << "  2-tag entries:        " << count2Tag << std::endl;

  return 0;
}
