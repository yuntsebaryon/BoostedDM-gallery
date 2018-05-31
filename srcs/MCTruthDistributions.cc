#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"


// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

  
using CounterMap_t = std::map< int, unsigned int >;

void ResetCounters( CounterMap_t& Multiplicity ) {

    Multiplicity[0] = 0; // All particles, including visible and invisible, but not GENIE artifacts, nor DMs, Ar40, Ar39, Cl39.
    Multiplicity[2212] = 0; // protons
    Multiplicity[2112] = 0; // neutrons
    Multiplicity[211] = 0;  // charged pions
    Multiplicity[111] = 0;  // neutral pions
    Multiplicity[300] = 0;  // Mesons including K0(311), K+(321), but not charged and neutral pions
    Multiplicity[3000] = 0; // Baryons including Lambda(3122), Sigma-(3112), Sigma+(3222), but not protons and neutrons
    Multiplicity[1000180400] = 0; // argon 40
    Multiplicity[1000180390] = 0; // argon 39
    Multiplicity[1000170390] = 0; // Cl 39
    Multiplicity[2000010000] = 0; // dark matter
    Multiplicity[2000000000] = 0; // GENIE artifacts
    
}
CounterMap_t* InitTreeGen( TTree* pTree ) {

    CounterMap_t* pMultipGen = new CounterMap_t;
    CounterMap_t& MultipGen = (*pMultipGen);
    
    ResetCounters( MultipGen );

    pTree->Branch( "nParticlesGen", &MultipGen[0], "nParticleGen/i" );
    pTree->Branch( "nProtonsGen", &MultipGen[2212], "nProtonsGen/i" );
    pTree->Branch( "nNeutronsGen", &MultipGen[2112], "nNeutronsGen/i" );
    pTree->Branch( "nPionsGen", &MultipGen[211], "nPionsGen/i" );
    pTree->Branch( "nPi0sGen", &MultipGen[111], "nPi0sGen/i" );
    pTree->Branch( "nMesonsGen", &MultipGen[300], "nMesonsGen/i" );
    pTree->Branch( "nBaryonsGen", &MultipGen[3000], "nBaryonsGen/i" );
    pTree->Branch( "nAr40Gen", &MultipGen[1000180400], "nAr40Gen/i" );
    pTree->Branch( "nAr39Gen", &MultipGen[1000180390], "nAr39Gen/i" );
    pTree->Branch( "nCl39Gen", &MultipGen[1000170390], "nCl39Gen/i" );
    pTree->Branch( "nDMsGen", &MultipGen[2000010000], "nDMsGen/i" );
    pTree->Branch( "nGENIEGen", &MultipGen[2000000000], "nGENIEGen/i" );
    return pMultipGen;
}

CounterMap_t* InitTreeG4( TTree* pTree ) {
    
    CounterMap_t* pMultipG4 = new CounterMap_t;
    CounterMap_t& MultipG4 = (*pMultipG4);
    
    ResetCounters( MultipG4 );

    pTree->Branch( "nParticlesG4", &MultipG4[0], "nParticlesG4/i" );
    pTree->Branch( "nProtonsG4", &MultipG4[2212], "nProtonsG4/i" );
    pTree->Branch( "nNeutronsG4", &MultipG4[2112], "nNeutronsG4/i" );
    pTree->Branch( "nPionsG4", &MultipG4[211], "nPionsG4/i" );
    pTree->Branch( "nPi0sG4", &MultipG4[111], "nPi0sG4/i" );
    pTree->Branch( "nMesonsG4", &MultipG4[300], "nMesonsG4/i" );
    pTree->Branch( "nBaryonsG4", &MultipG4[3000], "nBaryonsG4/i" );
    pTree->Branch( "nAr40G4", &MultipG4[1000180400], "nAr40G4/i" );
    pTree->Branch( "nAr39G4", &MultipG4[1000180390], "nAr39G4/i" );
    pTree->Branch( "nCl39G4", &MultipG4[1000170390], "nCl39G4/i" );
    pTree->Branch( "nDMsG4", &MultipG4[2000010000], "nDMsG4/i" );
    pTree->Branch( "nGENIEG4", &MultipG4[2000000000], "nGENIEG4/i" );
    return pMultipG4;
}

int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "gsimple";
    std::string G4Label = "largeant";
    art::InputTag MCTruthTag { GenLabel };
    art::InputTag MCParticleTag { G4Label };

    TFile *fOut = new TFile( "MCDistributions.root", "RECREATE" );
    TTree *fTree = new TTree( "MCParticles", "Primary MC Particles" );
   
    CounterMap_t *pMultiplicityGen = InitTreeGen( fTree );
    CounterMap_t& MultiplicityGen = *pMultiplicityGen;
    CounterMap_t *pMultiplicityG4 = InitTreeG4( fTree );
    CounterMap_t& MultiplicityG4 = *pMultiplicityG4;
    
    for ( gallery::Event ev( Filenames ); !ev.atEnd(); ev.next() ) {

        
        for ( auto& multPair: MultiplicityGen ) multPair.second = 0;
        
        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;
        art::FindMany< simb::MCParticle, sim::GeneratedParticleInfo > G4MCParticlesAssn( MCTruthHandle, ev, MCParticleTag );

        // std::map< int, double > LeadingMomentum;
        // std::map< int, double > LeadingEnergy;

        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {

            simb::MCTruth MCTruthObj = MCTruthObjs[iMCTruth];
            int nParticles = MCTruthObj.NParticles();

            for ( int iParticle = 0; iParticle < nParticles; ++iParticle ) {

                const simb::MCParticle& MCParticleObj = MCTruthObj.GetParticle( iParticle );
                int pdgCode = MCParticleObj.PdgCode();
                // std::cout << "MCParticle: " << iParticle << ", PDGCode: " << pdgCode << std::endl;
                // double momentum = MCParticleObj.P();
                // double energy = MCParticleObj.E();

                if ( pdgCode > 2000000000 && pdgCode != 2000010000 ) { 
                    MultiplicityGen[2000000000]++;
                    // if ( momentum > LeadingMomentum[2000000001] ) LeadingMomentum[2000000001] = momentum;
                    // if ( energy > LeadingEnergy[2000000001] ) LeadingEnergy[2000000001] = energy;
               
                } else if ( abs(pdgCode) > 300 && abs(pdgCode) < 400 ) {
                    MultiplicityGen[300]++;
                    MultiplicityGen[0]++;
                } else if ( abs(pdgCode) > 3000 && abs(pdgCode) < 4000 ) {
                    MultiplicityGen[3000]++;
                    MultiplicityGen[0]++;
                } else {
                    MultiplicityGen[abs(pdgCode)]++;
                    if ( pdgCode == 2000010000 || pdgCode == 1000180400 || pdgCode == 1000180390 || pdgCode == 1000170390 ) continue; 
                    MultiplicityGen[0]++;
                    // if ( momentum > LeadingMomentum[abs(pdgCode)] ) LeadingMomentum[abs(pdgCode)] = momentum;
                    // if ( energy > LeadingEnergy[abs(pdgCode)] ) LeadingEnergy[abs(pdgCode)] = energy;
                }
               
            } // Loop over MCParticles for each MCTruth
            
            std::vector< simb::MCParticle const* > const& G4MCParticles = G4MCParticlesAssn.at( iMCTruth );
            
            for ( size_t iMCParticle = 0; iMCParticle < G4MCParticles.size(); ++iMCParticle ) {
                int pdgCode = G4MCParticles[iMCParticle]->PdgCode();
                
                if ( pdgCode > 2000000000 && pdgCode != 2000010000 ) { 
                    MultiplicityG4[2000000000]++;
               
                } else if ( abs(pdgCode) > 300 && abs(pdgCode) < 400 ) {
                    MultiplicityG4[300]++;
                    MultiplicityG4[0]++;
                } else if ( abs(pdgCode) > 3000 && abs(pdgCode) < 4000 ) {
                    MultiplicityG4[3000]++;
                    MultiplicityG4[0]++;
                } else {
                    MultiplicityG4[abs(pdgCode)]++;
                    if ( pdgCode == 2000010000 || pdgCode == 1000180400 || pdgCode == 1000180390 || pdgCode == 1000170390 ) continue; 
                    MultiplicityG4[0]++;

                }
            }
        } // Loop over MCTruth 

        /* for ( std::map< int, int >::iterator it = MultiplicityGen.begin(); it != MultiplicityGen.end(); ++it ) {
            hMultiplicity[it->first]->Fill( it->second );
            if ( LeadingMomentum[it->first] > 0. ) hLeadingMomentum[it->first]->Fill( LeadingMomentum[it->first] );
            if ( LeadingEnergy[it->first] > 0. ) hLeadingEnergy[it->first]->Fill( LeadingEnergy[it->first] );
        } */
        fTree->Fill();
    } // End of an event

    fOut->Write();
    return 0;
}
