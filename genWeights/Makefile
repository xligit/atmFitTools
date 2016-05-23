# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(T2KREWEIGHT)/make/Make.include

all: genWeightsFromSK_NEUTFSI genWeightsFromSK_BANFF14c genWeights_2015c WeightOut genXsecSpline genWeights_sfrpa_rfg

genWeightsFromSK_NEUTFSI: FORCE
	@echo -e "\n\n***** Compiling genWeightsFromSK_NEUTFSI\n"
	$(CXX) $(CXXFLAGS) -c genWeightsFromSK_NEUTFSI.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) genWeightsFromSK_NEUTFSI.o $(LIBRARIES) -o genWeightsFromSK_NEUTFSI.exe

genWeightsFromSK_BANFF14c: FORCE
	@echo -e "\n\n***** Compiling genWeightsFromSK_BANFF14c\n"
	$(CXX) $(CXXFLAGS) -c genWeightsFromSK_BANFF14c.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) genWeightsFromSK_BANFF14c.o $(LIBRARIES) -o genWeightsFromSK_BANFF14c.exe

genWeights_2015c: FORCE
	@echo -e "\n\n***** Compiling genWeights_2015c\n"
	$(CXX) $(CXXFLAGS) -c genWeights_2015c.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) genWeights_2015c.o $(LIBRARIES) -o genWeights_2015c.exe

WeightOut:
	$(CXX) $(CXXFLAGS) -c WeightOut.cc $(INCLUDES)
	$(LD) $(LDFLAGS) WeightOut.o $(LIBRARIES) -o WeightOut.exe

genXsecSpline: FORCE
	@echo -e "\n\n*** Compiling genXsecSpline\n"
	$(CXX) $(CXXFLAGS) -c genXsecSpline.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) genXsecSpline.o $(LIBRARIES) -o genXsecSpline.exe

genWeights_sfrpa_rfg: FORCE
	@echo -e "\n\n*** Compiling genWeights_sfrpa_rfg\n"
	$(CXX) $(CXXFLAGS) -c genWeights_sfrpa_rfg.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) genWeights_sfrpa_rfg.o $(LIBRARIES) -o genWeights_sfrpa_rfg.exe

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) *.exe

FORCE:

# DO NOT DELETE
