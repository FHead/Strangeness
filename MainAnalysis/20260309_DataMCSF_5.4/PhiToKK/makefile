ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

default: ExecuteMakePhiSignalOnlyHistograms ExecuteFitPhiSignalOnlyShapes ExecuteMakePhiSBHistograms ExecuteFitPhiSB ExecuteFitPhiSBData ExecuteEvaluatePhiScaleFactors

ExecuteMakePhiSignalOnlyHistograms: MakePhiSignalOnlyHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakePhiSignalOnlyHistograms.cpp \
		-o ExecuteMakePhiSignalOnlyHistograms \
		$(ROOTLIBS)

ExecuteFitPhiSignalOnlyShapes: FitPhiSignalOnlyShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitPhiSignalOnlyShapes.cpp \
		-o ExecuteFitPhiSignalOnlyShapes \
		$(ROOTLIBS)

ExecuteMakePhiSBHistograms: MakePhiSBHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakePhiSBHistograms.cpp \
		-o ExecuteMakePhiSBHistograms \
		$(ROOTLIBS)

ExecuteFitPhiSB: FitPhiSB.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitPhiSB.cpp \
		-o ExecuteFitPhiSB \
		$(ROOTLIBS)

ExecuteFitPhiSBData: FitPhiSBData.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitPhiSBData.cpp \
		-o ExecuteFitPhiSBData \
		$(ROOTLIBS)

ExecuteEvaluatePhiScaleFactors: EvaluatePhiScaleFactors.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		EvaluatePhiScaleFactors.cpp \
		-o ExecuteEvaluatePhiScaleFactors \
		$(ROOTLIBS)
