#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>

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
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"


using CounterMap_t = std::map< int, unsigned int >;
using KinematicMap_t = std::map< int, std::vector< double > >;

void SetParticleTypes( std::vector< int >& ParticleTypes ) {

    ParticleTypes.push_back( 0 );  // All particles, including visible and invisible, but not GENIE artifacts, nor DMs, Ar40, Ar39, Cl39.
    ParticleTypes.push_back( 2212 );  // protons
    ParticleTypes.push_back( 2112 );  // neutrons
    ParticleTypes.push_back( 211  );  // charged pions
    ParticleTypes.push_back( 111  );  // neutral pions
    ParticleTypes.push_back( 300  );  // Mesons including K0(311), K+(321), but not charged and neutral pions
    ParticleTypes.push_back( 3000 );  // Baryons including Lambda(3122), Sigma-(3112), Sigma+(3222), but not protons and neutrons
    ParticleTypes.push_back( 1000180400 );  // argon 40
    ParticleTypes.push_back( 1000180390 );  // argon 39
    ParticleTypes.push_back( 1000170390 );  // Cl 39
    ParticleTypes.push_back( 2000010000 );  // dark matter
    ParticleTypes.push_back( 2000000000 );  // GENIE artifacts

}

void ResetCounters( CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E ) {

    std::vector< int > ParticleTypes;
    SetParticleTypes( ParticleTypes );

    for ( size_t iType = 0; iType < ParticleTypes.size(); ++iType ) {
        int type = ParticleTypes[iType];
        Multiplicity[type] = 0;
        Px[type].clear();
        Py[type].clear();
        Pz[type].clear();
        P[type].clear();
        E[type].clear();
    }
}

void InitTree( TTree* pTree, CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, double& EventPx, double& EventPy, double& EventPz, double& EventE, double& angle ) {

    ResetCounters( Multiplicity, Px, Py, Pz, P, E );

    pTree->Branch( "EventPx", &EventPx );
    pTree->Branch( "EventPy", &EventPy );
    pTree->Branch( "EventPz", &EventPz );
    pTree->Branch( "EventE", &EventE );
    pTree->Branch( "Angle", &angle );
    pTree->Branch( "nParticles", &Multiplicity[0], "nParticles/i" );
    pTree->Branch( "nProtons", &Multiplicity[2212], "nProtons/i" );
    pTree->Branch( "nNeutrons", &Multiplicity[2112], "nNeutrons/i" );
    pTree->Branch( "nPions", &Multiplicity[211], "nPions/i" );
    pTree->Branch( "nPi0s", &Multiplicity[111], "nPi0s/i" );
    pTree->Branch( "nMesons", &Multiplicity[300], "nMesons/i" );
    pTree->Branch( "nBaryons", &Multiplicity[3000], "nBaryons/i" );
    pTree->Branch( "nAr40", &Multiplicity[1000180400], "nAr40/i" );
    pTree->Branch( "nAr39", &Multiplicity[1000180390], "nAr39/i" );
    pTree->Branch( "nCl39", &Multiplicity[1000170390], "nCl39/i" );
    pTree->Branch( "nInDMs", &Multiplicity[2000010000], "nDMs/i" );
    pTree->Branch( "nGENIE", &Multiplicity[2000000000], "nGENIE/i" );
    pTree->Branch( "ProtonPx", &Px[2212] );
    pTree->Branch( "NeutronPx", &Px[2112] );
    pTree->Branch( "PionPx", &Px[211] );
    pTree->Branch( "Pi0Px", &Px[111] );
    pTree->Branch( "MesonPx", &Px[300] );
    pTree->Branch( "BaryonPx", &Px[3000] );
    pTree->Branch( "Ar40Px", &Px[1000180400] );
    pTree->Branch( "Ar39Px", &Px[1000180390] );
    pTree->Branch( "Cl39Px", &Px[1000170390] );
    pTree->Branch( "InDMPx", &Px[2000010000] );
    pTree->Branch( "GENIEPx", &Px[2000000000] );
    pTree->Branch( "ProtonPy", &Py[2212] );
    pTree->Branch( "NeutronPy", &Py[2112] );
    pTree->Branch( "PionPy", &Py[211] );
    pTree->Branch( "Pi0Py", &Py[111] );
    pTree->Branch( "MesonPy", &Py[300] );
    pTree->Branch( "BaryonPy", &Py[3000] );
    pTree->Branch( "Ar40Py", &Py[1000180400] );
    pTree->Branch( "Ar39Py", &Py[1000180390] );
    pTree->Branch( "Cl39Py", &Py[1000170390] );
    pTree->Branch( "InDMPy", &Py[2000010000] );
    pTree->Branch( "GENIEPy", &Py[2000000000] );
    pTree->Branch( "ProtonPz", &Pz[2212] );
    pTree->Branch( "NeutronPz", &Pz[2112] );
    pTree->Branch( "PionPz", &Pz[211] );
    pTree->Branch( "Pi0Pz", &Pz[111] );
    pTree->Branch( "MesonPz", &Pz[300] );
    pTree->Branch( "BaryonPz", &Pz[3000] );
    pTree->Branch( "Ar40Pz", &Pz[1000180400] );
    pTree->Branch( "Ar39Pz", &Pz[1000180390] );
    pTree->Branch( "Cl39Pz", &Pz[1000170390] );
    pTree->Branch( "InDMPz", &Pz[2000010000] );
    pTree->Branch( "GENIEPz", &Pz[2000000000] );
    pTree->Branch( "ProtonP", &P[2212] );
    pTree->Branch( "NeutronP", &P[2112] );
    pTree->Branch( "PionP", &P[211] );
    pTree->Branch( "Pi0P", &P[111] );
    pTree->Branch( "MesonP", &P[300] );
    pTree->Branch( "BaryonP", &P[3000] );
    pTree->Branch( "Ar40P", &P[1000180400] );
    pTree->Branch( "Ar39P", &P[1000180390] );
    pTree->Branch( "Cl39P", &P[1000170390] );
    pTree->Branch( "InDMP", &P[2000010000] );
    pTree->Branch( "GENIEP", &P[2000000000] );
    pTree->Branch( "ProtonE", &E[2212] );
    pTree->Branch( "NeutronE", &E[2112] );
    pTree->Branch( "PionE", &E[211] );
    pTree->Branch( "Pi0E", &E[111] );
    pTree->Branch( "MesonE", &E[300] );
    pTree->Branch( "BaryonE", &E[3000] );
    pTree->Branch( "Ar40E", &E[1000180400] );
    pTree->Branch( "Ar39E", &E[1000180390] );
    pTree->Branch( "Cl39E", &E[1000170390] );
    pTree->Branch( "InDME", &E[2000010000] );
    pTree->Branch( "GENIEE", &E[2000000000] );
}


int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "gsimple";
    std::string G4Label = "largeant";
    art::InputTag MCTruthTag { GenLabel };
    art::InputTag MCParticleTag { G4Label };

    TFile *fOut = new TFile( "MCParticleDistributionsTEST.root", "RECREATE" );
    TTree *fTree = new TTree( "MCParticles", "Primary MC Particles" );

    CounterMap_t Multiplicity;
    KinematicMap_t Px, Py, Pz, P, E;
    double EventPx, EventPy, EventPz, EventE, angle, counter = 0, summed = 0;
    InitTree( fTree, Multiplicity, Px, Py, Pz, P, E, EventPx, EventPy, EventPz, EventE, angle );


    for ( gallery::Event ev( Filenames ); !ev.atEnd(); ev.next() ) {

        // for ( auto& multPair: Multiplicity ) multPair.second = 0;
        ResetCounters( Multiplicity, Px, Py, Pz, P, E );
        EventPx = 0;
        EventPy = 0;
        EventPz = 0;
        EventE = 0;
        
        // default value of angle is -1
        angle = -1;

        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;
        // Find the MCParticles produced by LArG4 and associated to the MCTruth
        art::FindMany< simb::MCParticle, sim::GeneratedParticleInfo > G4MCParticlesAssn( MCTruthHandle, ev, MCParticleTag );

        // std::map< int, double > LeadingMomentum;
        // std::map< int, double > LeadingEnergy;
        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {

            simb::MCTruth MCTruthObj = MCTruthObjs[iMCTruth];
            // int nParticles = MCTruthObj.NParticles();

            // The incident DM particle is stored in the neutrino container in the MCTruth
            simb::MCNeutrino const& InDMObj = MCTruthObj.GetNeutrino();
            simb::MCParticle const& InDM    = InDMObj.Nu();
            Multiplicity[2000010000] = 1;
            auto& InDMPx = Px[2000010000]; auto& InDMPy = Py[2000010000]; auto& InDMPz = Pz[2000010000];
            auto& InDMP = P[2000010000]; auto& InDME = E[2000010000];
            InDMPx.resize( 1, 0. ); InDMPy.resize( 1, 0. ); InDMPz.resize( 1, 0. );
            InDMP.resize( 1, 0. ); InDME.resize( 1, 0. );
            InDMPx[0] = InDM.Px(); InDMPy[0] = InDM.Py(); InDMPz[0] = InDM.Pz();
            InDMP[0] = InDM.P(); InDME[0] = InDM.E();

            // Find all the (stable final state) MCParticles associated to the MCTruth object
            std::vector< simb::MCParticle const* > const& G4MCParticles = G4MCParticlesAssn.at( iMCTruth );
            // std::vector< sim::GeneratedParticleInfo const* > const& G4MCParticleInfo = G4MCParticlesAssn.data( iMCTruth );

            // keep track of leading momentum particle over all MCParticles to calculate angle with incident DM later
            // double leading_index;
            double leading_P = 0;
            double leading_Px = 0;
            double leading_Py = 0;
            double leading_Pz = 0;
            // const simb::MCParticle* leading_particle;

            for ( size_t iMCParticle = 0; iMCParticle < G4MCParticles.size(); ++iMCParticle ) {
                const simb::MCParticle* thisMCParticle = G4MCParticles[iMCParticle];
                int pdgCode = thisMCParticle->PdgCode();
                auto& nAllParticles = Multiplicity[0];

                // if the particle is proton, neutron, charged pion, pi0, meson, or baryon, add P and E to running sum
                if ( abs(pdgCode) == 2212 || abs(pdgCode) == 2112 || abs(pdgCode) == 211 || abs(pdgCode) == 111 || abs(pdgCode) == 300 || abs(pdgCode) == 3000 ) {
                    EventPx += thisMCParticle->Px();
                    EventPy += thisMCParticle->Py();
                    EventPz += thisMCParticle->Pz();
                    EventE += thisMCParticle->E();

                    // check to see if this MCParticle's P is leading
                    if (thisMCParticle->P() > leading_P) {
                      leading_P = thisMCParticle->P();
                      leading_Px = thisMCParticle->Px();
                      leading_Py = thisMCParticle->Py();
                      leading_Pz = thisMCParticle->Pz();
                      // leading_particle = thisMCParticle;
                    }
                }

                if ( pdgCode > 2000000000 && pdgCode != 2000010000 ) {
                    auto& nGENIE = Multiplicity[2000000000];
                    auto& GENIEPx = Px[2000000000];
                    if ( GENIEPx.size() <= nGENIE ) GENIEPx.resize( nGENIE + 1, 0. );
                    GENIEPx[nGENIE] = thisMCParticle->Px();
                    auto& GENIEPy = Py[2000000000];
                    if ( GENIEPy.size() <= nGENIE ) GENIEPy.resize( nGENIE + 1, 0. );
                    GENIEPy[nGENIE] = thisMCParticle->Py();
                    auto& GENIEPz = Pz[2000000000];
                    if ( GENIEPz.size() <= nGENIE ) GENIEPz.resize( nGENIE + 1, 0. );
                    GENIEPz[nGENIE] = thisMCParticle->Pz();
                    auto& GENIEP = P[2000000000];
                    if ( GENIEP.size() <= nGENIE ) GENIEP.resize( nGENIE + 1, 0. );
                    GENIEP[nGENIE] = thisMCParticle->P();
                    auto& GENIEE = E[2000000000];
                    if ( GENIEE.size() <= nGENIE ) GENIEE.resize( nGENIE + 1, 0. );
                    GENIEE[nGENIE] = thisMCParticle->E();

                    ++nGENIE;
                } else if ( abs(pdgCode) > 300 && abs(pdgCode) < 400 ) {
                    auto& nMesons = Multiplicity[300];
                    auto& MesonPx = Px[300];
                    if ( MesonPx.size() <= nMesons ) MesonPx.resize( nMesons + 1, 0. );
                    MesonPx[nMesons] = thisMCParticle->Px();
                    auto& MesonPy = Py[300];
                    if ( MesonPy.size() <= nMesons ) MesonPy.resize( nMesons + 1, 0. );
                    MesonPy[nMesons] = thisMCParticle->Py();
                    auto& MesonPz = Pz[300];
                    if ( MesonPz.size() <= nMesons ) MesonPz.resize( nMesons + 1, 0. );
                    MesonPz[nMesons] = thisMCParticle->Pz();
                    auto& MesonP = P[300];
                    if ( MesonP.size() <= nMesons ) MesonP.resize( nMesons + 1, 0. );
                    MesonP[nMesons] = thisMCParticle->P();
                    auto& MesonE = E[300];
                    if ( MesonE.size() <= nMesons ) MesonE.resize( nMesons + 1, 0. );
                    MesonE[nMesons] = thisMCParticle->E();
                    ++nMesons;
                    ++nAllParticles;
                } else if ( abs(pdgCode) > 3000 && abs(pdgCode) < 4000 ) {
                    auto& nBaryons = Multiplicity[3000];
                    auto& BaryonPx = Px[3000];
                    if ( BaryonPx.size() <= nBaryons ) BaryonPx.resize( nBaryons + 1, 0. );
                    BaryonPx[nBaryons] = thisMCParticle->Px();
                    auto& BaryonPy = Py[3000];
                    if ( BaryonPy.size() <= nBaryons ) BaryonPy.resize( nBaryons + 1, 0. );
                    BaryonPy[nBaryons] = thisMCParticle->Py();
                    auto& BaryonPz = Pz[3000];
                    if ( BaryonPz.size() <= nBaryons ) BaryonPz.resize( nBaryons + 1, 0. );
                    BaryonPz[nBaryons] = thisMCParticle->Pz();
                    auto& BaryonP = P[3000];
                    if ( BaryonP.size() <= nBaryons ) BaryonP.resize( nBaryons + 1, 0. );
                    BaryonP[nBaryons] = thisMCParticle->P();
                    auto& BaryonE = E[3000];
                    if ( BaryonE.size() <= nBaryons ) BaryonE.resize( nBaryons + 1, 0. );
                    BaryonE[nBaryons] = thisMCParticle->E();
                    ++nBaryons;
                    ++nAllParticles;
                } else {
                    auto& nParticles = Multiplicity[abs(pdgCode)];
                    auto& ParticlePx = Px[abs(pdgCode)];
                    if ( ParticlePx.size() <= nParticles ) ParticlePx.resize( nParticles + 1, 0. );
                    ParticlePx[nParticles] = thisMCParticle->Px();
                    auto& ParticlePy = Py[abs(pdgCode)];
                    if ( ParticlePy.size() <= nParticles ) ParticlePy.resize( nParticles + 1, 0. );
                    ParticlePy[nParticles] = thisMCParticle->Py();
                    auto& ParticlePz = Pz[abs(pdgCode)];
                    if ( ParticlePz.size() <= nParticles ) ParticlePz.resize( nParticles + 1, 0. );
                    ParticlePz[nParticles] = thisMCParticle->Pz();
                    auto& ParticleP = P[abs(pdgCode)];
                    if ( ParticleP.size() <= nParticles ) ParticleP.resize( nParticles + 1, 0. );
                    ParticleP[nParticles] = thisMCParticle->P();
                    auto& ParticleE = E[abs(pdgCode)];
                    if ( ParticleE.size() <= nParticles ) ParticleE.resize( nParticles + 1, 0. );
                    ParticleE[nParticles] = thisMCParticle->E();
                    ++nParticles;
                    if ( pdgCode == 2000010000 || pdgCode == 1000180400 || pdgCode == 1000180390 || pdgCode == 1000170390 ) continue;
                    ++nAllParticles;

                }
            } // Loop over MCParticles

            // now that we have leading particle, calculate angle between leading particle and incident DM
            // but we only fill "Angle" if we have well-defined leading momentum modulus
            // const simb::MCParticle* leading_particle = G4MCParticles[leading_index];
            if (leading_P > 0) {
              double numerator = leading_Px * InDM.Px() + leading_Py * InDM.Py() + leading_Pz * InDM.Pz();
              double denominator = leading_P * InDM.P();
              angle = acos(numerator/denominator);
              summed += angle;
              counter += 1;
            }
        } // Loop over MCTruth
        /* for ( std::map< int, int >::iterator it = MultiplicityGen.begin(); it != MultiplicityGen.end(); ++it ) {
            hMultiplicity[it->first]->Fill( it->second );
            if ( LeadingMomentum[it->first] > 0. ) hLeadingMomentum[it->first]->Fill( LeadingMomentum[it->first] );
            if ( LeadingEnergy[it->first] > 0. ) hLeadingEnergy[it->first]->Fill( LeadingEnergy[it->first] );
        } */
        fTree->Fill();
    } // End of an event
    std::cout << summed/10000 << std::endl;
    std::cout << summed/counter << std::endl;
    std::cout << counter << std::endl;

    fOut->Write();
    return 0;
}
