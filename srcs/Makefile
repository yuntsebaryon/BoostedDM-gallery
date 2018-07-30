#Makefile for gallery c++ programs.
#Note, being all-incllusive here: you can cut out libraries/includes you do not need
#you can also change the flags if you want too (Werror, pedantic, etc.)

CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CANVAS_ROOT_IO_INC) \
         -I $(CETLIB_INC) \
         -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(LARCOREALG_INC) \
         -I $(NUSIMDATA_INC) \
         -I $(ROOT_INC)

CXXFLAGS=-std=c++17 -Wall -Werror -pedantic
CXX=clang++
LDFLAGS=$$(root-config --libs) \
        -l GenVector \
        -L $(CANVAS_LIB) -l canvas \
        -L $(CANVAS_ROOT_IO_LIB) -l canvas_root_io \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase \
        -L $(LARCOREALG_LIB) -l larcorealg_Geometry

all: MCTruthDistributions MCParticleDistributions CheckMCTruthMCParticle MCDistributions

MCTruthDistributions: MCTruthDistributions.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

MCParticleDistributions: MCParticleDistributions.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

CheckMCTruthMCParticle: CheckMCTruthMCParticle.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

MCDistributions: MCDistributions.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

clean: rm *.o MCTruthDistributions MCParticleDistributions CheckMCTruthMCParticle MCDistributions
