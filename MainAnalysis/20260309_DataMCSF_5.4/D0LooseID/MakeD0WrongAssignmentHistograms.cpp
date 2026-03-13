#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kPionMass = 0.13957039;
constexpr double kMassWindowMin = 1.70;
constexpr double kMassWindowMax = 2.00;
constexpr int kMassBins = 320;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMatchAngleMax = 0.01;
constexpr double kTruthEnergyMatchMax = 0.025;
constexpr double kTruthMomentumAngleMax = 0.025;
constexpr long long kKaonTagThreshold = 1;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;
constexpr int kMaxD0 = 4096;
constexpr int kMaxGen = 20000;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
};

double angleBetween(double px1, double py1, double pz1,
                    double px2, double py2, double pz2) {
  const double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
  const double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
  if (p1 <= 0.0 || p2 <= 0.0) return std::numeric_limits<double>::infinity();
  double c = (px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return std::acos(c);
}

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
      getArgument(argc, argv, "--output", "d0_looseid_wrong_assignment_histograms.root");
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

  long long nD0 = 0;
  long long reco1ID[kMaxD0] = {0};
  long long reco2ID[kMaxD0] = {0};
  double reco1Angle[kMaxD0] = {0.0};
  double reco2Angle[kMaxD0] = {0.0};
  long long nGen = 0;
  double genPx[kMaxGen] = {0.0};
  double genPy[kMaxGen] = {0.0};
  double genPz[kMaxGen] = {0.0};
  double genE[kMaxGen] = {0.0};
  long long genID[kMaxGen] = {0};
  long long genStatus[kMaxGen] = {0};
  long long genMatchIndex[kMaxGen] = {0};
  double genMatchAngle[kMaxGen] = {0.0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  tree->SetBranchAddress("ND0", &nD0);
  tree->SetBranchAddress("D0Reco1ID[ND0]", reco1ID);
  tree->SetBranchAddress("D0Reco2ID[ND0]", reco2ID);
  tree->SetBranchAddress("D0Reco1Angle[ND0]", reco1Angle);
  tree->SetBranchAddress("D0Reco2Angle[ND0]", reco2Angle);
  tree->SetBranchAddress("NGen", &nGen);
  tree->SetBranchAddress("GenPx", genPx);
  tree->SetBranchAddress("GenPy", genPy);
  tree->SetBranchAddress("GenPz", genPz);
  tree->SetBranchAddress("GenE", genE);
  tree->SetBranchAddress("GenID", genID);
  tree->SetBranchAddress("GenStatus", genStatus);
  tree->SetBranchAddress("GenMatchIndex", genMatchIndex);
  tree->SetBranchAddress("GenMatchAngle", genMatchAngle);

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

    if (nD0 > 0) {
      for (long long i = 0; i < nD0; ++i) {
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
    } else {
      std::vector<long long> kPlus, kMinus, piPlus, piMinus;
      for (long long i = 0; i < nGen; ++i) {
        if (genStatus[i] != 1) continue;
        const long long reco = genMatchIndex[i];
        if (reco < 0 || reco >= nReco) continue;
        if (genMatchAngle[i] >= matchAngleMax) continue;
        if (recoGoodTrack[reco] != 1) continue;
        const TrackKinematics t{recoPx[reco], recoPy[reco], recoPz[reco]};
        if (!passAcceptance(t)) continue;
        if (genID[i] == 321 && recoCharge[reco] > 0) kPlus.push_back(i);
        if (genID[i] == -321 && recoCharge[reco] < 0) kMinus.push_back(i);
        if (genID[i] == 211 && recoCharge[reco] > 0) piPlus.push_back(i);
        if (genID[i] == -211 && recoCharge[reco] < 0) piMinus.push_back(i);
      }

      for (long long iD0 = 0; iD0 < nGen; ++iD0) {
        if (genID[iD0] != 421 && genID[iD0] != -421) continue;
        totalCandidates++;
        const bool isD0 = (genID[iD0] == 421);
        const std::vector<long long>& kaons = isD0 ? kMinus : kPlus;
        const std::vector<long long>& pions = isD0 ? piPlus : piMinus;

        long long bestK = -1;
        long long bestPi = -1;
        double bestAbsDeltaE = std::numeric_limits<double>::infinity();
        double bestAngle = std::numeric_limits<double>::infinity();

        for (long long iK : kaons) {
          for (long long iPi : pions) {
            const double absDeltaE = std::fabs(genE[iK] + genE[iPi] - genE[iD0]);
            if (absDeltaE >= bestAbsDeltaE) continue;
            const double sx = genPx[iK] + genPx[iPi];
            const double sy = genPy[iK] + genPy[iPi];
            const double sz = genPz[iK] + genPz[iPi];
            const double angle = angleBetween(sx, sy, sz, genPx[iD0], genPy[iD0], genPz[iD0]);
            bestAbsDeltaE = absDeltaE;
            bestAngle = angle;
            bestK = iK;
            bestPi = iPi;
          }
        }
        if (bestK < 0 || bestPi < 0) continue;
        if (bestAbsDeltaE >= kTruthEnergyMatchMax) continue;
        if (bestAngle >= kTruthMomentumAngleMax) continue;
        const long long i1 = genMatchIndex[bestK];
        const long long i2 = genMatchIndex[bestPi];
        if (i1 < 0 || i2 < 0 || i1 >= nReco || i2 >= nReco) continue;
        if (recoCharge[i1] * recoCharge[i2] >= 0) continue;
        passSelection++;

        const TrackKinematics d1{recoPx[i1], recoPy[i1], recoPz[i1]};
        const TrackKinematics d2{recoPx[i2], recoPy[i2], recoPz[i2]};
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
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hNominalNoPionTag.Write();
  hSwappedNoPionTag.Write();
  hNominalWithPionTag.Write();
  hSwappedWithPionTag.Write();

  TNamed selection("SelectionSummary",
                   Form("D0 LooseID wrong-assignment study: use D0* branches when available, otherwise for each truth "
                        "D0/D0bar choose the final-state Kpi pair minimizing |E(Kpi)-E(D0)|, require "
                        "|E(Kpi)-E(D0)|<%.3f and angle(Kpi,D0)<%.3f, then require GenMatchAngle<%.4f on both legs; "
                        "compare nominal daughter1->K daughter2->pi assignment against swapped daughter2->K "
                        "daughter1->pi assignment, with RecoPIDKaon>=1 on the assumed kaon and optional pion-tag "
                        "selections, hist range %.3f-%.3f GeV",
                        kTruthEnergyMatchMax, kTruthMomentumAngleMax, matchAngleMax, massMin, massMax));
  selection.Write();
  TParameter<long long>("TotalD0Candidates", totalCandidates).Write();
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
