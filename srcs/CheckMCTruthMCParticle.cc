#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"


// LArSoft, nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "dk2nu";
    std::string G4Label = "largeant";
    art::InputTag MCTruthTag { GenLabel };
    art::InputTag MCParticleTag { G4Label };

    for ( gallery::Event ev( Filenames ); !ev.atEnd(); ev.next() ) {

        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;
        
        // Find the MCParticles produced by LArG4 and associated to the MCTruth
        art::FindMany< simb::MCParticle, sim::GeneratedParticleInfo > G4MCParticlesAssn( MCTruthHandle, ev, MCParticleTag );

        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {
            
            simb::MCTruth MCTruthObj = MCTruthObjs[iMCTruth];
            int nParticles = MCTruthObj.NParticles();
            TLorentzVector EventMomentum, iP, fP;
            
            // Check the momentum and energy conservation in an interaction
            for ( int iParticle = 0; iParticle < nParticles; ++iParticle ) {

                const simb::MCParticle& MCParticleObj = MCTruthObj.GetParticle( iParticle );
                std::cout << "GENIE Particle[" << iParticle << "]: PdgCode: " << MCParticleObj.PdgCode() << ", StatudCode: " << MCParticleObj.StatusCode() << std::endl;
                const TLorentzVector& iMomentum = MCParticleObj.Momentum();
                // Initial state particles
                if ( MCParticleObj.StatusCode() == 0 ) iP += iMomentum;
                // Stable final state particles and nuclear remnants
                else if (  MCParticleObj.StatusCode() == 1 || MCParticleObj.StatusCode() == 15 ) fP += iMomentum;
            } // Loop over MCParticles for each MCTruth
            
            EventMomentum = fP - iP;
            if ( std::abs( EventMomentum.M2() ) > 1e-10 ) {
                std::cout << "MCTruth " << iMCTruth << " doesn't conserve the momentum and energy.  The total 4-momentum is ( " << EventMomentum.Px() << ", " << EventMomentum.Py() << ", " << EventMomentum.Pz() << ", " << EventMomentum.E() << " )." << std::endl;
            }
            
            // Find all the (stable final state) MCParticles associated to the MCTruth object
            std::vector< simb::MCParticle const* > const& G4MCParticles = G4MCParticlesAssn.at( iMCTruth );
            std::vector< sim::GeneratedParticleInfo const* > const& G4MCParticleInfo = G4MCParticlesAssn.data( iMCTruth );
            for ( size_t iMCParticle = 0; iMCParticle < G4MCParticles.size(); ++iMCParticle ) {

                const simb::MCParticle* thisMCParticle = G4MCParticles[iMCParticle];
                int pdgCode = thisMCParticle->PdgCode();
                int statusCode = thisMCParticle->StatusCode();
                int trackID = thisMCParticle->TrackId();
                if ( statusCode != 1 ) {
                    std::cout << "MCParticle[" << trackID << "] has the status code " << statusCode << "!" << std::endl;
                    continue;
                }
                
                // Sanity check: whether the MCParticle propagated to LArG4 is identical to that in GENIE (or other MC generators).
                // Check PdgCode, StatusCode, and 4-momentum.
                const sim::GeneratedParticleInfo* thisMCParticleInfo = G4MCParticleInfo[iMCParticle];
                const simb::MCParticle& thisMCParticleInGenerator = MCTruthObj.GetParticle( thisMCParticleInfo->generatedParticleIndex() );
                // Check if there is a corresponding GENIE index
                bool hasGeneratedParticleIndex = thisMCParticleInfo->hasGeneratedParticleIndex();
                if ( !hasGeneratedParticleIndex ) continue;
                int refpdgCode = thisMCParticleInGenerator.PdgCode();
                if ( pdgCode != refpdgCode ) {
                    std::cout << "MCParticle[" << trackID << "] has PDG code difference : GENIE(" << refpdgCode << "), GEANT(" << pdgCode << ")!" << std::endl;
                }
                int refstatusCode = thisMCParticleInGenerator.StatusCode();
                if ( statusCode != refstatusCode ) {
                    std::cout << "MCParticle[" << trackID << "] has Status code difference : GENIE(" << refstatusCode << "), GEANT(" << statusCode << ")!" << std::endl;
                }
                

                if ( thisMCParticle->Process() != "primary" ) {
                    std::cout << "MCParticle[" << trackID << "] has the process " << thisMCParticle->Process() << "!" << std::endl;
                    continue;
                }
                
                // std::cout << "MCParticle [" << trackID << "] has index in GENIE: " <<  thisMCParticleInfo->generatedParticleIndex() << std::endl;
                // std::cout << "Reference MCParticle 4-momentum: (" << thisMCParticleInGenerator.Px() << ", " << thisMCParticleInGenerator.Py() << ", " << thisMCParticleInGenerator.Pz() << ", " << thisMCParticleInGenerator.E() << ")" << std::endl;
                const TLorentzVector& Momentum = thisMCParticle->Momentum();
                const TLorentzVector& refMomentum = thisMCParticleInGenerator.Momentum();
                TLorentzVector diffMomentum = Momentum - refMomentum;
                if ( std::abs( diffMomentum.M2() ) > 1e-10 ) {
                    std::cout << "MCParticle[" << trackID << "] has 4-momentum difference: GENIE(" << refMomentum.Px() << ", " << refMomentum.Py() << ", " << refMomentum.Pz() << ", " << refMomentum.E() << "), GEANT(" << Momentum.Px() << ", " << Momentum.Py() << ", " << Momentum.Pz() << ", " << Momentum.E() << ")!" << std::endl;
                }
                // End of sanity checks ----------------------------------------------------------------------
                
            } // Loop over MCParticles
        } // Loop over MCTruth
    }
    return 0;
}
