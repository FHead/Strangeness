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
constexpr double kPionMass = 0.13957039;
constexpr double kMassWindowMin = 0.70;
constexpr double kMassWindowMax = 1.10;
constexpr int kMassBins = 320;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMatchAngleMax = 0.01;
constexpr long long kKaonTagThreshold = 2;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;
constexpr int kMaxKStar = 4096;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
};

double buildMass(const TrackKinematics& assumedKaon, const TrackKinematics& assumedPion) {
  const double pK2 = assumedKaon.px * assumedKaon.px + assumedKaon.py * assumedKaon.py + assumedKaon.pz * assumedKaon.pz;
  const double pPi2 = assumedPion.px * assumedPion.px + assumedPion.py * assumedPion.py + assumedPion.pz * assumedPion.pz;
  const double eK = std::sqrt(pK2 + kKaonMass * kKaonMass);
  const double ePi = std::sqrt(pPi2 + kPionMass * kPionMass);
  const double px = assumedKaon.px + assumedPion.px;
  const double py = assumedKaon.py + assumedPion.py;
  const double pz = assumedKaon.pz + assumedPion.pz;
  const double e = eK + ePi;
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
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName =
      getArgument(argc, argv, "--input", "../../../../Samples/merged_mc_v2.3.root");
  const std::string outputFileName =
      getArgument(argc, argv, "--output", "KStarWrongAssignmentHistograms.root");
  const std::string treeName = getArgument(argc, argv, "--tree", "Tree");
  const double massMin = getDoubleArgument(argc, argv, "--mass-min", kMassWindowMin);
  const double massMax = getDoubleArgument(argc, argv, "--mass-max", kMassWindowMax);
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
  long long recoPIDPion[kMaxReco] = {0};
  long long recoGoodTrack[kMaxReco] = {0};

  long long nKStar = 0;
  long long reco1ID[kMaxKStar] = {0};
  long long reco2ID[kMaxKStar] = {0};
  double reco1Angle[kMaxKStar] = {0.0};
  double reco2Angle[kMaxKStar] = {0.0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  tree->SetBranchAddress("NKStar", &nKStar);
  tree->SetBranchAddress("KStarReco1ID[NKStar]", reco1ID);
  tree->SetBranchAddress("KStarReco2ID[NKStar]", reco2ID);
  tree->SetBranchAddress("KStarReco1Angle[NKStar]", reco1Angle);
  tree->SetBranchAddress("KStarReco2Angle[NKStar]", reco2Angle);

  TH1D hNominalNoPionTag("hNominalNoPionTag",
                         "Nominal assignment, kaon-tag only; m(K#pi) [GeV]; Candidates / bin",
                         kMassBins, massMin, massMax);
  TH1D hSwappedNoPionTag("hSwappedNoPionTag",
                         "Swapped assignment, kaon-tag only; m(K#pi) [GeV]; Candidates / bin",
                         kMassBins, massMin, massMax);
  TH1D hNominalWithPionTag("hNominalWithPionTag",
                           "Nominal assignment, kaon+pion tag; m(K#pi) [GeV]; Candidates / bin",
                           kMassBins, massMin, massMax);
  TH1D hSwappedWithPionTag("hSwappedWithPionTag",
                           "Swapped assignment, kaon+pion tag; m(K#pi) [GeV]; Candidates / bin",
                           kMassBins, massMin, massMax);

  long long totalCandidates = 0;
  long long passSelection = 0;
  long long countNominalNoPionTag = 0;
  long long countSwappedNoPionTag = 0;
  long long countNominalWithPionTag = 0;
  long long countSwappedWithPionTag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    for (long long i = 0; i < nKStar; ++i) {
      totalCandidates++;

      const long long i1 = reco1ID[i];
      const long long i2 = reco2ID[i];
      if (i1 < 0 || i2 < 0 || i1 >= nReco || i2 >= nReco) continue;
      if (reco1Angle[i] >= matchAngleMax || reco2Angle[i] >= matchAngleMax) continue;
      if (recoGoodTrack[i1] != 1 || recoGoodTrack[i2] != 1) continue;

      const TrackKinematics d1{recoPx[i1], recoPy[i1], recoPz[i1]};
      const TrackKinematics d2{recoPx[i2], recoPy[i2], recoPz[i2]};
      if (!passAcceptance(d1) || !passAcceptance(d2)) continue;
      if (recoCharge[i1] * recoCharge[i2] >= 0) continue;
      passSelection++;

      if (recoPIDKaon[i1] >= kKaonTagThreshold) {
        hNominalNoPionTag.Fill(buildMass(d1, d2));
        countNominalNoPionTag++;
        if (recoPIDPion[i2] >= kPionTagThreshold) {
          hNominalWithPionTag.Fill(buildMass(d1, d2));
          countNominalWithPionTag++;
        }
      }

      if (recoPIDKaon[i2] >= kKaonTagThreshold) {
        hSwappedNoPionTag.Fill(buildMass(d2, d1));
        countSwappedNoPionTag++;
        if (recoPIDPion[i1] >= kPionTagThreshold) {
          hSwappedWithPionTag.Fill(buildMass(d2, d1));
          countSwappedWithPionTag++;
        }
      }
    }
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hNominalNoPionTag.Write();
  hSwappedNoPionTag.Write();
  hNominalWithPionTag.Write();
  hSwappedWithPionTag.Write();

  TNamed selection("SelectionSummary",
                   Form("KStar wrong-assignment study: valid KStarReco IDs, KStarReco1Angle<%.4f, "
                        "KStarReco2Angle<%.4f, both RecoGoodTrack==1, both 0.15<=|cos(theta)|<=0.675, "
                        "opposite charge, compare nominal daughter1->K daughter2->pi assignment against "
                        "swapped daughter2->K daughter1->pi assignment, with kaon-tag and optional pion-tag selections, "
                        "hist range %.3f-%.3f GeV",
                        matchAngleMax, matchAngleMax, massMin, massMax));
  selection.Write();
  TParameter<long long>("TotalKStarCandidates", totalCandidates).Write();
  TParameter<long long>("PassSelection", passSelection).Write();
  TParameter<long long>("CountNominalNoPionTag", countNominalNoPionTag).Write();
  TParameter<long long>("CountSwappedNoPionTag", countSwappedNoPionTag).Write();
  TParameter<long long>("CountNominalWithPionTag", countNominalWithPionTag).Write();
  TParameter<long long>("CountSwappedWithPionTag", countSwappedWithPionTag).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Pass common selection:       " << passSelection << std::endl;
  std::cout << "  Nominal, no pion tag:        " << countNominalNoPionTag << std::endl;
  std::cout << "  Swapped, no pion tag:        " << countSwappedNoPionTag << std::endl;
  std::cout << "  Nominal, with pion tag:      " << countNominalWithPionTag << std::endl;
  std::cout << "  Swapped, with pion tag:      " << countSwappedWithPionTag << std::endl;
  return 0;
}
