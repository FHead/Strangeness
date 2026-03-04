ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

ExecuteKtoPiAnalysis: KtoPiAnalysis.cpp
	 g++ -O3 -I. -I$(ProjectBase)/CommonCode/include -I./include \
	    $(ROOTCFLAGS) \
	    KtoPiAnalysis.cpp \
	    $(ProjectBase)/CommonCode/library/StrangenessMessenger.o \
	    -o ExecuteKtoPiAnalysis \
	    $(ROOTLIBS)
