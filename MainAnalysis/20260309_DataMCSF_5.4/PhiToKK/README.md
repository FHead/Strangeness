# Phi -> KK analysis

This subfolder holds the `phi -> KK` mass-spectrum fit workflow used to derive kaon efficiency.

Current contents:
- `MakePhiSignalOnlyHistograms.cpp`: builds signal-only MC `K^{+}K^{-}` mass histograms for the 1-tag and 2-tag categories.
- `MakePhiSBHistograms.cpp`: builds reco-only same-event `K^{+}K^{-}` mass histograms.
- `FitPhiSignalOnlyShapes.cpp`: scans candidate signal-only line shapes.
- `FitPhiSB.cpp`: fits reco-only MC signal-plus-background spectra.
- `FitPhiSBData.cpp`: fits reco-only data signal-plus-background spectra.
- `RunPhiSystematics.sh`: runs the current set of `phi` systematic variations.
- `Systematics/`: output area for systematic-variation runs.
- `makefile`: builds the standalone executables used in this folder.
