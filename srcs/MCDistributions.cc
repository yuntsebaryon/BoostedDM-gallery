/***
 * @file    MCDistributions.cc
 * @brief   Create a TTree filling with the 4-momenta of the initial and final state particles of the dark matter particle-argon interaction
 * @date    August 5, 2018
 * @author  Yun-Tse Tsai (yuntse@slac.stanford.edu), Dane Stocks (dstocks@stanford.edu)
 * @version v1.0
 *
 * Currently we assume there is only one interaction, which is a DM-Argon interaction, in each event.
 * Therefore there should be only one MCTruth object in each event.
 * We will have to reconsider the scope of the TTree entries (per event or per interaction, how to count particles in two interactions nearby, etc.) when we start to consider multiple interactions (multiple DM-argon interactions, or one DM-argon interaction and multiple atmospheric neutrino interactions in an event, or both).
 * When we look for the leading particle or the leading proton, we are comparing the modulus 3-momentum.
 *
 ***/

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
// #include <cassert>

// some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"
//#include "Math/GenVector/LorentzVector.h"

// "art" includes (canvas, and gallery)
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
#include "larcorealg/Geometry/geo_vectors_utils.h"


using CounterMap_t = std::map< int, unsigned int >;
using KinematicMap_t = std::map< int, std::vector< double > >;

using Momentum4_t = ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D< double > >;

void SetParticleTypes( std::vector< int >& ParticleTypes ) {

    ParticleTypes.push_back( 0 );  // All particles, including visible and invisible
    ParticleTypes.push_back( 2212 );  // protons
    ParticleTypes.push_back( 2112 );  // neutrons
    ParticleTypes.push_back( 211  );  // charged pions
    ParticleTypes.push_back( 111  );  // neutral pions
    ParticleTypes.push_back( 300  );  // Mesons including K0(311), K+(321), KL(130), but not charged and neutral pions
    ParticleTypes.push_back( 3000 );  // Baryons including Lambda(3122), Sigma-(3112), Sigma+(3222), but not protons and neutrons
    ParticleTypes.push_back( 1000180400 );  // argon 40
    ParticleTypes.push_back( 1000180390 );  // argon 39
    ParticleTypes.push_back( 1000170390 );  // Cl 39
    ParticleTypes.push_back( -2000010000 );  // Incident dark matter
    ParticleTypes.push_back( 2000010000 ); // Outgoing dark matter
    ParticleTypes.push_back( 2000000000 );  // GENIE artifacts
    ParticleTypes.push_back( 2000000400 ); // All the visible particles and neutrons ( not GENIE artifacts, nor DMs, Ar40, Ar39, Cl39 )
    ParticleTypes.push_back( 2000000401 ); // All the visible particles but not neutrons
    ParticleTypes.push_back( 2000000410 ); // Leading particle among all the visible particles and neutrons
    ParticleTypes.push_back( 2000000411 ); // Leading particle among all the visible particles but not neutrons
    ParticleTypes.push_back( 2000000412 ); // Leading proton
}

void ResetCounters( CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle ) {

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
        Angle[type].clear();
    }
}

void ResizeKinematics( std::vector< double >& Px, std::vector< double >& Py, std::vector< double >& Pz, std::vector< double >& P, std::vector< double >& E, std::vector< double >& Angle, size_t n ) {

    Px.resize( n, 0. );
    Py.resize( n, 0. );
    Pz.resize( n, 0. );
    P.resize( n, 0. );
    E.resize( n, 0. );
    // The angle between the vector sum of all the visible final state particles and the incident DM particle
    // Set the default value to -1.
    Angle.resize( n, -1. );

}

// CalculateAngle accepts two vectors (in our case, 3-vectors) and calculates the angle between them
double CalculateAngle( std::vector< double > vec1, std::vector< double > vec2 ) {
  if ( vec1.size() != vec2.size() ) {
    std::cout << "Vectors are not the same length! Aborting..." << std::endl;
    return 0.;
  }
  else {
    double numerator = 0, denominator = 0, vec1mag = 0, vec2mag = 0, angle;
    for (int i = 0; i < vec1.size(); ++i) {
      numerator += vec1[i] * vec2[i];
      vec1mag += pow(vec1[i], 2);
      vec2mag += pow(vec2[i], 2);
    }
    vec1mag = pow(vec1mag, 0.5);
    vec2mag = pow(vec2mag, 0.5);
    denominator = vec1mag * vec2mag;
    angle = acos(numerator/denominator);
    return angle;
  }
}

void InitTree( TTree* pTree, CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle, simb::Origin_t& EventOrigin, int& Mode, int& InteractionType, int& CCNC ) {

    ResetCounters( Multiplicity, Px, Py, Pz, P, E, Angle );

    // TODO: Fill all the information into TTree

    // Event origin: neutrino, cosmics, etc.
    pTree->Branch( "EventOrigin", &EventOrigin );
    // Interaction mode
    pTree->Branch( "Mode", &Mode );
    // Interaction type
    pTree->Branch( "InteractionType", &InteractionType );
    // CC or NC
    pTree->Branch( "CCNC", &CCNC );

    // Use the index "0" for the event-wide four-momentum
    pTree->Branch( "EventPx", &Px[0] );
    pTree->Branch( "EventPy", &Py[0] );
    pTree->Branch( "EventPz", &Pz[0] );
    pTree->Branch( "EventP", &P[0] );
    pTree->Branch( "EventE", &E[0] );

    pTree->Branch( "VisiblePx", &Px[2000000400] );
    pTree->Branch( "VisibleNoNPx", &Px[2000000401] );
    pTree->Branch( "LeadingParticlePx", &Px[2000000410] );
    pTree->Branch( "LeadingParticleNoNPx", &Px[2000000411] );
    pTree->Branch( "LeadingProtonPx", &Px[2000000412] );

    pTree->Branch( "VisiblePy", &Py[2000000400] );
    pTree->Branch( "VisibleNoNPy", &Py[2000000401] );
    pTree->Branch( "LeadingParticlePy", &Py[2000000410] );
    pTree->Branch( "LeadingParticleNoNPy", &Py[2000000411] );
    pTree->Branch( "LeadingProtonPy", &Py[2000000412] );

    pTree->Branch( "VisiblePz", &Pz[2000000400] );
    pTree->Branch( "VisibleNoNPz", &Pz[2000000401] );
    pTree->Branch( "LeadingParticlePz", &Pz[2000000410] );
    pTree->Branch( "LeadingParticleNoNPz", &Pz[2000000411] );
    pTree->Branch( "LeadingProtonPz", &Pz[2000000412] );

    pTree->Branch( "VisibleP", &P[2000000400] );
    pTree->Branch( "VisibleNoNP", &P[2000000401] );
    pTree->Branch( "LeadingParticleP", &P[2000000410] );
    pTree->Branch( "LeadingParticleNoNP", &P[2000000411] );
    pTree->Branch( "LeadingProtonP", &P[2000000412] );

    pTree->Branch( "VisibleE", &E[2000000400] );
    pTree->Branch( "VisibleNoNE", &E[2000000401] );
    pTree->Branch( "LeadingParticleE", &E[2000000410] );
    pTree->Branch( "LeadingParticleNoNE", &E[2000000411] );
    pTree->Branch( "LeadingProtonE", &E[2000000412] );

    pTree->Branch( "EventAngle", &Angle[0] );
    pTree->Branch( "VisibleAngle", &Angle[2000000400] );
    pTree->Branch( "VisibleNoNAngle", &Angle[2000000401] );
    pTree->Branch( "LeadingParticleAngle", &Angle[2000000410] );
    pTree->Branch( "LeadingParticleNoNAngle", &Angle[2000000411] );
    pTree->Branch( "LeadingProtonAngle", &Angle[2000000412] );


    // Event-wide particle multiplicity
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
    pTree->Branch( "nInDMs", &Multiplicity[-2000010000], "nInDMs/i" );
    pTree->Branch( "nOutDMs", &Multiplicity[2000010000], "nOutDMs/i" );
    pTree->Branch( "nGENIE", &Multiplicity[2000000000], "nGENIE/i" );

    // Variables for each final state particle
    pTree->Branch( "ProtonPx", &Px[2212] );
    pTree->Branch( "NeutronPx", &Px[2112] );
    pTree->Branch( "PionPx", &Px[211] );
    pTree->Branch( "Pi0Px", &Px[111] );
    pTree->Branch( "MesonPx", &Px[300] );
    pTree->Branch( "BaryonPx", &Px[3000] );
    pTree->Branch( "Ar40Px", &Px[1000180400] );
    pTree->Branch( "Ar39Px", &Px[1000180390] );
    pTree->Branch( "Cl39Px", &Px[1000170390] );
    pTree->Branch( "InDMPx", &Px[-2000010000] );
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
    pTree->Branch( "InDMPy", &Py[-2000010000] );
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
    pTree->Branch( "InDMPz", &Pz[-2000010000] );
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
    pTree->Branch( "InDMP", &P[-2000010000] );
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
    pTree->Branch( "InDME", &E[-2000010000] );
    pTree->Branch( "GENIEE", &E[2000000000] );
}


int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );

    std::string GenLabel = "gsimple";
    std::string G4Label = "largeant";
    art::InputTag MCTruthTag { GenLabel };
    art::InputTag MCParticleTag { G4Label };

    TFile *fOut = new TFile( "MCParticleDistributions.root", "RECREATE" );
    TTree *fTree = new TTree( "MCParticles", "Primary MC Particles" );

    CounterMap_t Multiplicity;
    KinematicMap_t Px, Py, Pz, P, E, Angle;
    simb::Origin_t EventOrigin;
    int Mode, InteractionType, CCNC;

    InitTree( fTree, Multiplicity, Px, Py, Pz, P, E, Angle, EventOrigin, Mode, InteractionType, CCNC );


    for ( gallery::Event ev( Filenames ); !ev.atEnd(); ev.next() ) {

        // for ( auto& multPair: Multiplicity ) multPair.second = 0;
        ResetCounters( Multiplicity, Px, Py, Pz, P, E, Angle );

        std::cout << "Processing "
                  << "Run " << ev.eventAuxiliary().run() << ", "
                  << "Event " << ev.eventAuxiliary().event() << std::endl;

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;

        // Find the MCParticles produced by LArG4 and associated to the MCTruth
        art::FindMany< simb::MCParticle, sim::GeneratedParticleInfo > G4MCParticlesAssn( MCTruthHandle, ev, MCParticleTag );

        // Access the element of the map directly
        auto& nAllParticles = Multiplicity[0];
        auto& nInDM = Multiplicity[-2000010000];
        auto& nVisible = Multiplicity[2000000400];
        auto& nVisibleNoN = Multiplicity[2000000401];

        // Initialize the incident DM particle.
        // Currently allow only one DM-Argon interaction in an event
        auto& InDMPx = Px[-2000010000]; auto& InDMPy = Py[-2000010000]; auto& InDMPz = Pz[-2000010000];
        auto& InDMP = P[-2000010000]; auto& InDME = E[-2000010000]; auto& InDMAngle = Angle[-2000010000];
        ResizeKinematics( InDMPx, InDMPy, InDMPz, InDMP, InDME, InDMAngle, 1 );

        auto& nGENIE = Multiplicity[2000000000];
        auto& GENIEPx = Px[2000000000]; auto& GENIEPy = Py[2000000000];
        auto& GENIEPz = Pz[2000000000]; auto& GENIEP = P[2000000000];
        auto& GENIEE = E[2000000000]; auto& GENIEAngle = Angle[2000000000];

        auto& nMesons = Multiplicity[300];
        auto& MesonPx = Px[300]; auto& MesonPy = Py[300]; auto& MesonPz = Pz[300];
        auto& MesonP = P[300]; auto& MesonE = E[300]; auto& MesonAngle = Angle[300];

        auto& nBaryons = Multiplicity[3000];
        auto& BaryonPx = Px[3000]; auto& BaryonPy = Py[3000]; auto& BaryonPz = Pz[3000];
        auto& BaryonP = P[3000]; auto& BaryonE = E[3000]; auto& BaryonAngle = Angle[3000];

        double LeadingParticleMax = 0., LeadingParticleNoNMax = 0., LeadingProtonMax = 0.;

        // The default constructor gives the initial values 0.
        Momentum4_t Event, Visible, VisibleNoN, LeadingParticle, LeadingParticleNoN, LeadingProton;

        for ( size_t iMCTruth = 0; iMCTruth < MCTruthObjs.size(); ++iMCTruth ) {

            // Assume only one incident DM particle, or one interaction in an event

            simb::MCTruth const& MCTruthObj = MCTruthObjs[iMCTruth];

            // The incident DM particle is stored in the neutrino container in the MCTruth
            simb::MCNeutrino const& InDMObj = MCTruthObj.GetNeutrino();
            simb::MCParticle const& InDM    = InDMObj.Nu();

            // Skip the loop in the case that the incident particle is not a DM particle
            if ( InDM.PdgCode() != 2000010000 ) {
                std::cout << "The incident particle in the MCTruth " << iMCTruth << " is a neutrino, skip the interaction..." << std::endl;
                continue;
            }

            // Currently deal with only one DM-Argon interaction in an event
            if ( nInDM > 0 ) {
                std::cout << "Found more than 1 DM-Argon interaction in the event!" << std::endl;
                continue;
            }

            // Interaction information
            Mode = InDMObj.Mode();
            InteractionType = InDMObj.InteractionType();
            CCNC = InDMObj.CCNC();

            // 4-momentum of the incident DM particle
            const auto InDMMom = geo::vect::convertTo< Momentum4_t >( InDM.Momentum() );
            Event += InDMMom;

            InDMPx[0] = InDM.Px(); InDMPy[0] = InDM.Py(); InDMPz[0] = InDM.Pz();
            InDMP[0] = InDM.P(); InDME[0] = InDM.E(); InDMAngle[0] = 0.;

            // Now look for the initial argon 4-momentum
            for ( size_t iMCParticle = 0; iMCParticle < MCTruthObj.NParticles(); ++iMCParticle ) {

                const simb::MCParticle& thisMCParticle = MCTruthObj.GetParticle( iMCParticle );
                if ( thisMCParticle.StatusCode() == 1 ) continue;
                // The incident DM particle has been filled
                if ( thisMCParticle.StatusCode() == 0 ) {
                    if ( thisMCParticle.PdgCode() == 2000010000 ) continue;
                    Event += geo::vect::convertTo< Momentum4_t >( thisMCParticle.Momentum() );
                    break;
                }
            }

            // Find all the (stable final state) MCParticles associated to the MCTruth object
            std::vector< simb::MCParticle const* > const& G4MCParticles = G4MCParticlesAssn.at( iMCTruth );

            for ( size_t iMCParticle = 0; iMCParticle < G4MCParticles.size(); ++iMCParticle ) {

                const simb::MCParticle* thisMCParticle = G4MCParticles[iMCParticle];
                int pdgCode = thisMCParticle->PdgCode();
                int statusCode = thisMCParticle->StatusCode();

                if ( statusCode != 1 ) {
                    std::cout << "MCParticle[" << thisMCParticle->TrackId() << "] has the status code " << statusCode << "!" << std::endl;
                }

                const auto Momentum = geo::vect::convertTo< Momentum4_t >( thisMCParticle->Momentum() );
                Event += Momentum;

                if ( pdgCode > 2000000000 && pdgCode < 2000000301 ) {

                    if ( GENIEPx.size() <= nGENIE )
                        ResizeKinematics( GENIEPx, GENIEPy, GENIEPz, GENIEP, GENIEE, GENIEAngle, nGENIE + 1 );

                    GENIEPx[nGENIE] = thisMCParticle->Px();
                    GENIEPy[nGENIE] = thisMCParticle->Py();
                    GENIEPz[nGENIE] = thisMCParticle->Pz();
                    GENIEP[nGENIE] = thisMCParticle->P();
                    GENIEE[nGENIE] = thisMCParticle->E();

                    ++nGENIE;
                    ++nAllParticles;

                } else if ( ( abs(pdgCode) > 300 && abs(pdgCode) < 400 ) || abs(pdgCode) == 130 ) {

                    if ( MesonPx.size() <= nMesons )
                        ResizeKinematics( MesonPx, MesonPy, MesonPz, MesonP, MesonE, MesonAngle, nMesons + 1 );

                    MesonPx[nMesons] = thisMCParticle->Px();
                    MesonPy[nMesons] = thisMCParticle->Py();
                    MesonPz[nMesons] = thisMCParticle->Pz();
                    MesonP[nMesons] = thisMCParticle->P();
                    MesonE[nMesons] = thisMCParticle->E();

                    Visible += Momentum;
                    VisibleNoN += Momentum;
                    if ( Momentum.P() > LeadingParticleMax ) {
                        LeadingParticle = Momentum;
                        LeadingParticleMax = Momentum.P();
                    }
                    if ( Momentum.P() > LeadingParticleNoNMax ) {
                        LeadingParticleNoN = Momentum;
                        LeadingParticleNoNMax = Momentum.P();
                    }

                    ++nMesons;
                    ++nVisible;
                    ++nVisibleNoN;
                    ++nAllParticles;

                } else if ( abs(pdgCode) > 3000 && abs(pdgCode) < 4000 ) {

                    if ( BaryonPx.size() <= nBaryons )
                        ResizeKinematics( BaryonPx, BaryonPy, BaryonPz, BaryonP, BaryonE, BaryonAngle, nBaryons + 1 );

                    BaryonPx[nBaryons] = thisMCParticle->Px();
                    BaryonPy[nBaryons] = thisMCParticle->Py();
                    BaryonPz[nBaryons] = thisMCParticle->Pz();
                    BaryonP[nBaryons] = thisMCParticle->P();
                    BaryonE[nBaryons] = thisMCParticle->E();

                    Visible += Momentum;
                    VisibleNoN += Momentum;
                    if ( Momentum.P() > LeadingParticleMax ) {
                        LeadingParticle = Momentum;
                        LeadingParticleMax = Momentum.P();
                    }
                    if ( Momentum.P() > LeadingParticleNoNMax ) {
                        LeadingParticleNoN = Momentum;
                        LeadingParticleNoNMax = Momentum.P();
                    }

                    ++nBaryons;
                    ++nVisible;
                    ++nVisibleNoN;
                    ++nAllParticles;

                } else {

                    auto& nParticles = Multiplicity[abs(pdgCode)];
                    auto& ParticlePx = Px[abs(pdgCode)]; auto& ParticlePy = Py[abs(pdgCode)];
                    auto& ParticlePz = Pz[abs(pdgCode)]; auto& ParticleP = P[abs(pdgCode)];
                    auto& ParticleE = E[abs(pdgCode)]; auto& ParticleAngle = Angle[abs(pdgCode)];

                    if ( ParticlePx.size() <= nParticles )
                        ResizeKinematics( ParticlePx, ParticlePy, ParticlePz, ParticleP, ParticleE, ParticleAngle, nParticles + 1 );

                    ParticlePx[nParticles] = thisMCParticle->Px();
                    ParticlePy[nParticles] = thisMCParticle->Py();
                    ParticlePz[nParticles] = thisMCParticle->Pz();
                    ParticleP[nParticles] = thisMCParticle->P();
                    ParticleE[nParticles] = thisMCParticle->E();

                    ++nParticles;
                    ++nAllParticles;

                    if ( pdgCode == 2000010000 ) continue;
                    Visible += Momentum;
                    if ( Momentum.P() > LeadingParticleMax ) {
                        LeadingParticle = Momentum;
                        LeadingParticleMax = Momentum.P();
                    }
                    ++nVisible;

                    if ( pdgCode == 2112 ) continue;
                    VisibleNoN += Momentum;
                    if ( Momentum.P() > LeadingParticleNoNMax ) {
                        LeadingParticleNoN = Momentum;
                        LeadingParticleNoNMax = Momentum.P();
                    }
                    ++nVisibleNoN;

                    if ( pdgCode == 2212 && Momentum.P() > LeadingProtonMax ) {
                        LeadingProton = Momentum;
                        LeadingProtonMax = Momentum.P();
                    }
                }
            } // Loop over MCParticles
            // calculate the relevant angles between event-wide values: InDM and EventP, VisibleP, VisibleNoNP, etc.
            // unsure how to deal with ROOT::Math::DisplacementVector3D<Cartesian3D<Scalar> > objects right now so converting
            // to vector -- fix later. (Dane)
            std::vector< double > InDMMom3Vec = { InDMMom.Px(), InDMMom.Py(), InDMMom.Pz() };
            std::vector< double > Event3Vec = { Event.Px(), Event.Py(), Event.Pz() };
            std::vector< double > Visible3Vec = { Visible.Px(), Visible.Py(), Visible.Pz() };
            std::vector< double > VisibleNoN3Vec = { VisibleNoN.Px(), VisibleNoN.Py(), VisibleNoN.Pz() };
            std::vector< double > LeadingParticle3Vec = { LeadingParticle.Px(), LeadingParticle.Py(), LeadingParticle.Pz() };
            std::vector< double > LeadingParticleNoN3Vec = { LeadingParticleNoN.Px(), LeadingParticleNoN.Py(), LeadingParticleNoN.Pz() };
            std::vector< double > LeadingProton3Vec = { LeadingProton.Px(), LeadingProton.Py(), LeadingProton.Pz() };

            Angle[0][0]          = CalculateAngle( InDMMom3Vec, Event3Vec );
            Angle[2000000400][0] = CalculateAngle( InDMMom3Vec, Visible3Vec );
            Angle[2000000401][0] = CalculateAngle( InDMMom3Vec, VisibleNoN3Vec );
            Angle[2000000410][0] = CalculateAngle( InDMMom3Vec, LeadingParticle3Vec );
            Angle[2000000411][0] = CalculateAngle( InDMMom3Vec, LeadingParticleNoN3Vec );
            Angle[2000000412][0] = CalculateAngle( InDMMom3Vec, LeadingProton3Vec );
        } // Loop over MCTruth

        // Fill the event-wide information
        // TODO: finish the filling
        Px[0][0] = Event.Px();
        Py[0][0] = Event.Py();
        Pz[0][0] = Event.Pz();
        P[0][0]  = Event.P();
        E[0][0]  = Event.E();

        Px[2000000400][0] = Visible.Px();
        Px[2000000401][0] = VisibleNoN.Px();
        Px[2000000410][0] = LeadingParticle.Px();
        Px[2000000411][0] = LeadingParticleNoN.Px();
        Px[2000000412][0] = LeadingProton.Px();

        Py[2000000400][0] = Visible.Py();
        Py[2000000401][0] = VisibleNoN.Py();
        Py[2000000410][0] = LeadingParticle.Py();
        Py[2000000411][0] = LeadingParticleNoN.Py();
        Py[2000000412][0] = LeadingProton.Py();

        Pz[2000000400][0] = Visible.Pz();
        Pz[2000000401][0] = VisibleNoN.Pz();
        Pz[2000000410][0] = LeadingParticle.Pz();
        Pz[2000000411][0] = LeadingParticleNoN.Pz();
        Pz[2000000412][0] = LeadingProton.Pz();

        P[2000000400][0] = Visible.P();
        P[2000000401][0] = VisibleNoN.P();
        P[2000000410][0] = LeadingParticle.P();
        P[2000000411][0] = LeadingParticleNoN.P();
        P[2000000412][0] = LeadingProton.P();

        E[2000000400][0] = Visible.E();
        E[2000000401][0] = VisibleNoN.E();
        E[2000000410][0] = LeadingParticle.E();
        E[2000000411][0] = LeadingParticleNoN.E();
        E[2000000412][0] = LeadingProton.E();

        fTree->Fill();
    } // End of an event

    fOut->Write();
    return 0;
}
