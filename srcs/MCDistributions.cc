/***
 * @file    MCDistributions.cc
 * @brief   Create a TTree filling with the 4-momenta of the initial and final state particles of the dark matter particle-argon interaction
 * @date    August 7, 2018
 * @author  Yun-Tse Tsai (yuntse@slac.stanford.edu), Dane Stocks (dstocks@stanford.edu)
 * @update  May 10, 2019.  Allow to deal with both the signal and background events
 * @version v2.0
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
#include <fstream>
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
#include "Math/GenVector/DisplacementVector3D.h"

// CLHEP libraries
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

// "art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
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

// BDM libraries
#include "boosteddmanalysis/DataObjects/SmearedMCParticle.h"


using CounterMap_t = std::map< int, unsigned int >;
using KinematicMap_t = std::map< int, std::vector< double > >;

using Momentum4_t = ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D< double > >;
using Momentum3_t = ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D< double > >;

void SetParticleTypes( std::vector< int >& ParticleTypes ) {

    // For MCParticle information
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
    // For smeared MC particles
    ParticleTypes.push_back( -2212 );  // protons
    ParticleTypes.push_back( -2112 );  // neutrons
    ParticleTypes.push_back( -211  );  // charged pions
    ParticleTypes.push_back( -13  );   // muons
    ParticleTypes.push_back( -11  );   // electrons
    ParticleTypes.push_back( -22  );   // photons
    ParticleTypes.push_back( -2000000400 ); // All the smeared visible particles and smeared neutrons
    ParticleTypes.push_back( -2000000401 ); // All the smeared visible particles but not neutrons
    ParticleTypes.push_back( -2000000402 ); // All the smeared reconstructable particles and reconstructable neutrons
    ParticleTypes.push_back( -2000000403 ); // All the smeared reconstructable particles but not neutrons
    ParticleTypes.push_back( -2000000410 ); // Leading smeared particle among all the smeared visible particles and smeared neutrons
    ParticleTypes.push_back( -2000000411 ); // Leading smeared particle among all the smeared visible particles but not neutrons
    ParticleTypes.push_back( -2000000412 ); // Leading smeared particle among all the smeared reconstructable particles and reconstructable neutrons
    ParticleTypes.push_back( -2000000413 ); // Leading smeared particle among all the smeared reconstructable particles but not neutrons
    ParticleTypes.push_back( -2000000414 ); // Leading smeared proton
    ParticleTypes.push_back( -2000000415 ); // Leading smeared reconstructable proton
    
    
}

void ResetCounters( CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle, double& SunPx, double& SunPy, double& SunPz, double& VtxX, double& VtxY, double& VtxZ, double& T0 ) {

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
    SunPx = 0.; SunPy = 0.; SunPz = 0.;
    VtxX = 0.; VtxY = 0.; VtxZ = 0.; T0 = 0.;
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

double CalculateAngle( Momentum3_t v1, Momentum3_t v2 ) {
  
  double v1dotv2 = v1.Dot( v2 );
  double cosTheta = v1dotv2 / ( v1.R() * v2.R() );
  return acos( cosTheta );
}

void InitTree( TTree* pTree, CounterMap_t& Multiplicity, KinematicMap_t& Px, KinematicMap_t& Py, KinematicMap_t& Pz, KinematicMap_t& P, KinematicMap_t& E, KinematicMap_t& Angle, int& run, int& eventNo, simb::Origin_t& EventOrigin, int& Mode, int& InteractionType, int& CCNC, int& isIn10kton, double& SunPx, double& SunPy, double& SunPz, double& VtxX, double& VtxY, double& VtxZ, double& T0 ) {

    ResetCounters( Multiplicity, Px, Py, Pz, P, E, Angle, SunPx, SunPy, SunPz, VtxX, VtxY, VtxZ, T0 );

    pTree->Branch( "Run", &run );
    pTree->Branch( "Event", &eventNo );
    
    // Event origin: neutrino, cosmics, etc.
    pTree->Branch( "EventOrigin", &EventOrigin );
    // Interaction mode
    pTree->Branch( "Mode", &Mode );
    // Interaction type
    pTree->Branch( "InteractionType", &InteractionType );
    // CC or NC
    pTree->Branch( "CCNC", &CCNC );
    
    // Whether the event is in the 10k ton volumn
    pTree->Branch( "isIn10kton", &isIn10kton );

    // Sun position
    pTree->Branch( "SunPx", &SunPx );
    pTree->Branch( "SunPy", &SunPy );
    pTree->Branch( "SunPz", &SunPz );

    // Use the index "0" for the event-wide four-momentum
    pTree->Branch( "EventPx", &Px[0] );
    pTree->Branch( "EventPy", &Py[0] );
    pTree->Branch( "EventPz", &Pz[0] );
    pTree->Branch( "EventP", &P[0] );
    pTree->Branch( "EventE", &E[0] );
    pTree->Branch( "EventAngle", &Angle[0] );

    // Interaction vertex position and time
    pTree->Branch( "VertexX", &VtxX );
    pTree->Branch( "VertexY", &VtxY );
    pTree->Branch( "VertexZ", &VtxZ );
    pTree->Branch( "T0", &T0 );
    
    pTree->Branch( "VisiblePx", &Px[2000000400] );
    pTree->Branch( "VisiblePy", &Py[2000000400] );
    pTree->Branch( "VisiblePz", &Pz[2000000400] );
    pTree->Branch( "VisibleP", &P[2000000400] );
    pTree->Branch( "VisibleE", &E[2000000400] );
    pTree->Branch( "VisibleAngle", &Angle[2000000400] );

    pTree->Branch( "VisibleNoNPx", &Px[2000000401] );
    pTree->Branch( "VisibleNoNPy", &Py[2000000401] );
    pTree->Branch( "VisibleNoNPz", &Pz[2000000401] );
    pTree->Branch( "VisibleNoNP", &P[2000000401] );
    pTree->Branch( "VisibleNoNE", &E[2000000401] );
    pTree->Branch( "VisibleNoNAngle", &Angle[2000000401] );

    pTree->Branch( "LeadingParticlePx", &Px[2000000410] );
    pTree->Branch( "LeadingParticlePy", &Py[2000000410] );
    pTree->Branch( "LeadingParticlePz", &Pz[2000000410] );
    pTree->Branch( "LeadingParticleP", &P[2000000410] );
    pTree->Branch( "LeadingParticleE", &E[2000000410] );
    pTree->Branch( "LeadingParticleAngle", &Angle[2000000410] );

    pTree->Branch( "LeadingParticleNoNPx", &Px[2000000411] );
    pTree->Branch( "LeadingParticleNoNPy", &Py[2000000411] );
    pTree->Branch( "LeadingParticleNoNPz", &Pz[2000000411] );
    pTree->Branch( "LeadingParticleNoNP", &P[2000000411] );
    pTree->Branch( "LeadingParticleNoNE", &E[2000000411] );
    pTree->Branch( "LeadingParticleNoNAngle", &Angle[2000000411] );

    pTree->Branch( "LeadingProtonPx", &Px[2000000412] );
    pTree->Branch( "LeadingProtonPy", &Py[2000000412] );
    pTree->Branch( "LeadingProtonPz", &Pz[2000000412] );
    pTree->Branch( "LeadingProtonP", &P[2000000412] );
    pTree->Branch( "LeadingProtonE", &E[2000000412] );
    pTree->Branch( "LeadingProtonAngle", &Angle[2000000412] );

    pTree->Branch( "SmearedVisiblePx", &Px[-2000000400] );
    pTree->Branch( "SmearedVisiblePy", &Py[-2000000400] );
    pTree->Branch( "SmearedVisiblePz", &Pz[-2000000400] );
    pTree->Branch( "SmearedVisibleP", &P[-2000000400] );
    pTree->Branch( "SmearedVisibleE", &E[-2000000400] );
    pTree->Branch( "SmearedVisibleAngle", &Angle[-2000000400] );
    
    pTree->Branch( "SmearedVisibleNonPx", &Px[-2000000401] );
    pTree->Branch( "SmearedVisibleNonPy", &Py[-2000000401] );
    pTree->Branch( "SmearedVisibleNonPz", &Pz[-2000000401] );
    pTree->Branch( "SmearedVisibleNonP", &P[-2000000401] );
    pTree->Branch( "SmearedVisibleNonE", &E[-2000000401] );
    pTree->Branch( "SmearedVisibleNonAngle", &Angle[-2000000401] );

    pTree->Branch( "SmearedReconstructablePx", &Px[-2000000402] );
    pTree->Branch( "SmearedReconstructablePy", &Py[-2000000402] );
    pTree->Branch( "SmearedReconstructablePz", &Pz[-2000000402] );
    pTree->Branch( "SmearedReconstructableP", &P[-2000000402] );
    pTree->Branch( "SmearedReconstructableE", &E[-2000000402] );
    pTree->Branch( "SmearedReconstructableAngle", &Angle[-2000000402] );

    pTree->Branch( "SmearedReconstructableNoNPx", &Px[-2000000403] );
    pTree->Branch( "SmearedReconstructableNoNPy", &Py[-2000000403] );
    pTree->Branch( "SmearedReconstructableNoNPz", &Pz[-2000000403] );
    pTree->Branch( "SmearedReconstructableNoNP", &P[-2000000403] );
    pTree->Branch( "SmearedReconstructableNoNE", &E[-2000000403] );
    pTree->Branch( "SmearedReconstructableNoNAngle", &Angle[-2000000403] );
    
    pTree->Branch( "LeadingSmearedPx", &Px[-2000000410] );
    pTree->Branch( "LeadingSmearedPy", &Py[-2000000410] );
    pTree->Branch( "LeadingSmearedPz", &Pz[-2000000410] );
    pTree->Branch( "LeadingSmearedP", &P[-2000000410] );
    pTree->Branch( "LeadingSmearedE", &E[-2000000410] );
    pTree->Branch( "LeadingSmearedAngle", &Angle[-2000000410] );
    
    pTree->Branch( "LeadingSmearedNoNPx", &Px[-2000000411] );
    pTree->Branch( "LeadingSmearedNoNPy", &Py[-2000000411] );
    pTree->Branch( "LeadingSmearedNoNPz", &Pz[-2000000411] );
    pTree->Branch( "LeadingSmearedNoNP", &P[-2000000411] );
    pTree->Branch( "LeadingSmearedNoNE", &E[-2000000411] );
    pTree->Branch( "LeadingSmearedNoNAngle", &Angle[-2000000411] );

    pTree->Branch( "LeadingSmearedReconstructablePx", &Px[-2000000412] );
    pTree->Branch( "LeadingSmearedReconstructablePy", &Py[-2000000412] );
    pTree->Branch( "LeadingSmearedReconstructablePz", &Pz[-2000000412] );
    pTree->Branch( "LeadingSmearedReconstructableP", &P[-2000000412] );
    pTree->Branch( "LeadingSmearedReconstructableE", &E[-2000000412] );
    pTree->Branch( "LeadingSmearedReconstructableAngle", &Angle[-2000000412] );
    
    pTree->Branch( "LeadingSmearedReconstructableNoNPx", &Px[-2000000413] );
    pTree->Branch( "LeadingSmearedReconstructableNoNPy", &Py[-2000000413] );
    pTree->Branch( "LeadingSmearedReconstructableNoNPz", &Pz[-2000000413] );
    pTree->Branch( "LeadingSmearedReconstructableNoNP", &P[-2000000413] );
    pTree->Branch( "LeadingSmearedReconstructableNoNE", &E[-2000000413] );
    pTree->Branch( "LeadingSmearedReconstructableNoNAngle", &Angle[-2000000413] );

    pTree->Branch( "LeadingSmearedProtonPx", &Px[-2000000414] );
    pTree->Branch( "LeadingSmearedProtonPy", &Py[-2000000414] );
    pTree->Branch( "LeadingSmearedProtonPz", &Pz[-2000000414] );
    pTree->Branch( "LeadingSmearedProtonP", &P[-2000000414] );
    pTree->Branch( "LeadingSmearedProtonE", &E[-2000000414] );
    pTree->Branch( "LeadingSmearedProtonAngle", &Angle[-2000000414] );
    
    pTree->Branch( "LeadingSmearedReconstructableProtonPx", &Px[-2000000415] );
    pTree->Branch( "LeadingSmearedReconstructableProtonPy", &Py[-2000000415] );
    pTree->Branch( "LeadingSmearedReconstructableProtonPz", &Pz[-2000000415] );
    pTree->Branch( "LeadingSmearedReconstructableProtonP", &P[-2000000415] );
    pTree->Branch( "LeadingSmearedReconstructableProtonE", &E[-2000000415] );
    pTree->Branch( "LeadingSmearedReconstructableProtonAngle", &Angle[-2000000415] );
    
    
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
    pTree->Branch( "nInParticles", &Multiplicity[-2000010000], "nInPaticles/i" );
    pTree->Branch( "nOutParticles", &Multiplicity[2000010000], "nOutParticles/i" );
    pTree->Branch( "nGENIE", &Multiplicity[2000000000], "nGENIE/i" );
    pTree->Branch( "nVisible", &Multiplicity[2000000400], "nVisible/i" );
    pTree->Branch( "nVisibleNoN", &Multiplicity[2000000401], "nVisibleNoN/i" );
    pTree->Branch( "nSelProtons", &Multiplicity[-2212], "nSelProtons/i" );
    pTree->Branch( "nSelNeutrons", &Multiplicity[-2112], "nSelNeutrons/i" );
    pTree->Branch( "nSelPions", &Multiplicity[-211], "nSelPions/i" );
    pTree->Branch( "nSelMuons", &Multiplicity[-13], "nSelMuons/i" );
    pTree->Branch( "nSelElectrons", &Multiplicity[-11], "nSelElectrons/i" );
    pTree->Branch( "nSelPhotons", &Multiplicity[-22], "nSelPhotons/i" );
    pTree->Branch( "nSmearedVisible", &Multiplicity[-2000000400], "nSmearedVisiable/i" );
    pTree->Branch( "nSmearedVisibleNoN", &Multiplicity[-2000000401], "nSmearedVisibleNoN/i" );
    pTree->Branch( "nSmearedReconstructable", &Multiplicity[-2000000402], "nSmearedReconstructable/i" );
    pTree->Branch( "nSmearedReconstructableNoN", &Multiplicity[-2000000403], "nSmearedReconstructableNoN/i" );
    
    

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
    pTree->Branch( "InParticlePx", &Px[-2000010000] );
    pTree->Branch( "OutParticlePx", &Px[2000010000] );
    pTree->Branch( "GENIEPx", &Px[2000000000] );
    pTree->Branch( "SmearedProtonPx", &Px[-2212] );
    pTree->Branch( "SmearedNeutronPx", &Px[-2112] );
    pTree->Branch( "SmearedPionPx", &Px[-211] );
    pTree->Branch( "SmearedMuonPx", &Px[-13] );
    pTree->Branch( "SmearedElectronPx", &Px[-11] );
    pTree->Branch( "SmearedPhotonPx", &Px[-22] );    
    pTree->Branch( "ProtonPy", &Py[2212] );
    pTree->Branch( "NeutronPy", &Py[2112] );
    pTree->Branch( "PionPy", &Py[211] );
    pTree->Branch( "Pi0Py", &Py[111] );
    pTree->Branch( "MesonPy", &Py[300] );
    pTree->Branch( "BaryonPy", &Py[3000] );
    pTree->Branch( "Ar40Py", &Py[1000180400] );
    pTree->Branch( "Ar39Py", &Py[1000180390] );
    pTree->Branch( "Cl39Py", &Py[1000170390] );
    pTree->Branch( "InParticlePy", &Py[-2000010000] );
    pTree->Branch( "OutParticlePy", &Py[2000010000] );
    pTree->Branch( "GENIEPy", &Py[2000000000] );
    pTree->Branch( "SmearedProtonPy", &Py[-2212] );
    pTree->Branch( "SmearedNeutronPy", &Py[-2112] );
    pTree->Branch( "SmearedPionPy", &Py[-211] );
    pTree->Branch( "SmearedMuonPy", &Py[-13] );
    pTree->Branch( "SmearedElectronPy", &Py[-11] );
    pTree->Branch( "SmearedPhotonPy", &Py[-22] );   
    pTree->Branch( "ProtonPz", &Pz[2212] );
    pTree->Branch( "NeutronPz", &Pz[2112] );
    pTree->Branch( "PionPz", &Pz[211] );
    pTree->Branch( "Pi0Pz", &Pz[111] );
    pTree->Branch( "MesonPz", &Pz[300] );
    pTree->Branch( "BaryonPz", &Pz[3000] );
    pTree->Branch( "Ar40Pz", &Pz[1000180400] );
    pTree->Branch( "Ar39Pz", &Pz[1000180390] );
    pTree->Branch( "Cl39Pz", &Pz[1000170390] );
    pTree->Branch( "InParticlePz", &Pz[-2000010000] );
    pTree->Branch( "OutParticlePz", &Pz[2000010000] );
    pTree->Branch( "GENIEPz", &Pz[2000000000] );
    pTree->Branch( "SmearedProtonPz", &Pz[-2212] );
    pTree->Branch( "SmearedNeutronPz", &Pz[-2112] );
    pTree->Branch( "SmearedPionPz", &Pz[-211] );
    pTree->Branch( "SmearedMuonPz", &Pz[-13] );
    pTree->Branch( "SmearedElectronPz", &Pz[-11] );
    pTree->Branch( "SmearedPhotonPz", &Pz[-22] );   
    pTree->Branch( "ProtonP", &P[2212] );
    pTree->Branch( "NeutronP", &P[2112] );
    pTree->Branch( "PionP", &P[211] );
    pTree->Branch( "Pi0P", &P[111] );
    pTree->Branch( "MesonP", &P[300] );
    pTree->Branch( "BaryonP", &P[3000] );
    pTree->Branch( "Ar40P", &P[1000180400] );
    pTree->Branch( "Ar39P", &P[1000180390] );
    pTree->Branch( "Cl39P", &P[1000170390] );
    pTree->Branch( "InParticleP", &P[-2000010000] );
    pTree->Branch( "OutParticleP", &P[2000010000] );
    pTree->Branch( "GENIEP", &P[2000000000] );
    pTree->Branch( "SmearedProtonP", &P[-2212] );
    pTree->Branch( "SmearedNeutronP", &P[-2112] );
    pTree->Branch( "SmearedPionP", &P[-211] );
    pTree->Branch( "SmearedMuonP", &P[-13] );
    pTree->Branch( "SmearedElectronP", &P[-11] );
    pTree->Branch( "SmearedPhotonP", &P[-22] );   
    pTree->Branch( "ProtonE", &E[2212] );
    pTree->Branch( "NeutronE", &E[2112] );
    pTree->Branch( "PionE", &E[211] );
    pTree->Branch( "Pi0E", &E[111] );
    pTree->Branch( "MesonE", &E[300] );
    pTree->Branch( "BaryonE", &E[3000] );
    pTree->Branch( "Ar40E", &E[1000180400] );
    pTree->Branch( "Ar39E", &E[1000180390] );
    pTree->Branch( "Cl39E", &E[1000170390] );
    pTree->Branch( "InParticleE", &E[-2000010000] );
    pTree->Branch( "OutParticleE", &E[2000010000] );
    pTree->Branch( "GENIEE", &E[2000000000] );
    pTree->Branch( "SmearedProtonE", &E[-2212] );
    pTree->Branch( "SmearedNeutronE", &E[-2112] );
    pTree->Branch( "SmearedPionE", &E[-211] );
    pTree->Branch( "SmearedMuonE", &E[-13] );
    pTree->Branch( "SmearedElectronE", &E[-11] );
    pTree->Branch( "SmearedPhotonE", &E[-22] );   
    pTree->Branch( "ProtonAngle", &Angle[2212] );
    pTree->Branch( "NeutronAngle", &Angle[2112] );
    pTree->Branch( "PionAngle", &Angle[211] );
    pTree->Branch( "Pi0Angle", &Angle[111] );
    pTree->Branch( "MesonAngle", &Angle[300] );
    pTree->Branch( "BaryonAngle", &Angle[3000] );
    pTree->Branch( "InParticleAngle", &Angle[-2000010000] );
    pTree->Branch( "OutParticleAngle", &Angle[2000010000] );
    pTree->Branch( "GENIEAngle", &Angle[2000000000] );
    pTree->Branch( "SmearedProtonAngle", &Angle[-2212] );
    pTree->Branch( "SmearedNeutronAngle", &Angle[-2112] );
    pTree->Branch( "SmearedPionAngle", &Angle[-211] );
    pTree->Branch( "SmearedMuonAngle", &Angle[-13] );
    pTree->Branch( "SmearedElectronAngle", &Angle[-11] );
    pTree->Branch( "SmearedPhotonAngle", &Angle[-22] );   
}

std::vector< Momentum3_t > LoadSunPositions( std::string SunFilename ) {

  std::ifstream SunFile( SunFilename );
  if ( !SunFile.is_open() ) throw std::runtime_error( "Can't open file '" + SunFilename + "'" );
  
  double x, y, z;
  unsigned int iLine = 0;
  std::vector< Momentum3_t > SunPositions;
  
  while ( SunFile ) {
    ++iLine;
    SunFile >> x >> y >> z;
    if ( SunFile.eof() ) break;
    if ( SunFile.fail() ) throw std::runtime_error( "Can't read x, y or z at line " + std::to_string( iLine ) );
    Momentum3_t SunPos( x, y, z );
    SunPositions.push_back( SunPos );
  }
  std::cout << "Read " << iLine << " lines from '" << SunFilename << "'" << std::endl;
  
  return SunPositions;
}

int main( int argc, char ** argv ) {

    std::vector< std::string > Filenames;
    Filenames.push_back( argv[1] );
    
    std::string outFilename( argv[2] );

    bool isBDM = true;
    if ( argc > 3 ) {
      std::string sampleType( argv[3] );
      if ( sampleType.find("background") || sampleType.find("neutrino") ) isBDM = false;
    }
    
    if ( argc < 4 && !isBDM ) {
      std::cerr << "Need the sun position input for the atmospheric neutrino events!" << std::endl;
      return 1;
    }
    
    std::string SunFilename;
    std::vector< Momentum3_t > SunPositions;
    if ( !isBDM ) {
      SunFilename = std::string( argv[4] );
      // std::cout << "Sun position file: " << SunFilename << std::endl;
      SunPositions = LoadSunPositions( SunFilename );
      // std::cout << "Number of sun positions: " << SunPositions.size();
      // std::cout << "1st: ( " << SunPositions[0].X() << ", " << SunPositions[0].Y() << ", " << SunPositions[0].Z() << " )" << std::endl;
      // std::cout << "2nd: ( " << SunPositions[1].X() << ", " << SunPositions[1].Y() << ", " << SunPositions[1].Z() << " )" << std::endl;
      // std::cout << "3rd: ( " << SunPositions[2].X() << ", " << SunPositions[2].Y() << ", " << SunPositions[2].Z() << " )" << std::endl;
      // std::cout << "4th: ( " << SunPositions[3].X() << ", " << SunPositions[3].Y() << ", " << SunPositions[3].Z() << " )" << std::endl;
    }
    
    std::string GenLabel = "dk2nu";
    if ( !isBDM ) GenLabel = "generator";

    std::string G4Label = "largeant";
    // std::string SelMCLabel = "selector";
    std::string SmearedMCLabel = "smear";
    art::InputTag MCTruthTag { GenLabel };
    art::InputTag MCParticleTag { G4Label };
    // art::InputTag SelMCParticleTag { SelMCLabel };
    art::InputTag SmearedMCTag { SmearedMCLabel };

    
    TFile *fOut = new TFile( outFilename.c_str(), "RECREATE" );
    TTree *fTree = new TTree( "MCParticles", "Primary MC particles and smeared stable final state particles" );

    CounterMap_t Multiplicity;
    KinematicMap_t Px, Py, Pz, P, E, Angle;
    simb::Origin_t EventOrigin;
    int run, eventNo, Mode, InteractionType, CCNC, isIn10kton;
    double SunPx, SunPy, SunPz, VtxX, VtxY, VtxZ, T0;
    int nIn10kton = 0;
    
    InitTree( fTree, Multiplicity, Px, Py, Pz, P, E, Angle, run, eventNo, EventOrigin, Mode, InteractionType, CCNC, isIn10kton, SunPx, SunPy, SunPz, VtxX, VtxY, VtxZ, T0 );

    // Random engine for picking an entry for the sun position
    int RandMax = SunPositions.size();
    CLHEP::HepJamesRandom engine;
    CLHEP::RandFlat flat( engine, 0, RandMax );
    
    for ( gallery::Event ev( Filenames ); !ev.atEnd(); ev.next() ) {

        // for ( auto& multPair: Multiplicity ) multPair.second = 0;
        ResetCounters( Multiplicity, Px, Py, Pz, P, E, Angle, SunPx, SunPy, SunPz, VtxX, VtxY, VtxZ, T0 );
        run = ev.eventAuxiliary().run();
        eventNo = ev.eventAuxiliary().event();
        
        std::cout << "Processing "
                  << "Run " << run << ", "
                  << "Event " << eventNo << std::endl;

        // In the background events, the sun position should be the same in the entire events.
        // It should be the same for the signal events as well, but we obtain the sun position of the signal BDM events from the incoming direction of the BDM, and therefore will be filled later.
        Momentum3_t SunPosition;
        if ( !isBDM ) {
          double randnumber = flat();
          int iSunPosition = floor( randnumber );
          SunPosition = SunPositions[iSunPosition];
        
        }

        auto const& MCTruthHandle = ev.getValidHandle< std::vector< simb::MCTruth > >( MCTruthTag );
        auto const& MCTruthObjs = *MCTruthHandle;

        // Find the MCParticles produced by LArG4 and associated to the MCTruth
        art::FindMany< simb::MCParticle, sim::GeneratedParticleInfo > G4MCParticlesAssn( MCTruthHandle, ev, MCParticleTag );

        // Access the selected MC particles: protons, muons, charged pions, electrons, photons, and possibly neutrons (depending on the configuration)
        // auto const& SelMCHandle = ev.getValidHandle< std::vector< art::Ptr< simb::MCParticle > > >( SelMCParticleTag );
        // auto const& SelMCObj    = *SelMCHandle;

        // Access the smeared selected MC particles
        auto const& SmearedMCHandle = ev.getValidHandle< std::vector< bdm::SmearedMCParticle > >( SmearedMCTag );
        auto const& SmearedMCParts   = *SmearedMCHandle;


        // Access the element of the map directly
        auto& nAllParticles = Multiplicity[0];
        auto& nInParticles = Multiplicity[-2000010000];
        auto& nVisible = Multiplicity[2000000400];
        auto& nVisibleNoN = Multiplicity[2000000401];

        // Initialize the incident DM particle.
        // Currently allow only one DM-Argon interaction in an event
        auto& InParticlePx = Px[-2000010000]; auto& InParticlePy = Py[-2000010000]; auto& InParticlePz = Pz[-2000010000];
        auto& InParticleP = P[-2000010000]; auto& InParticleE = E[-2000010000]; auto& InParticleAngle = Angle[-2000010000];
        ResizeKinematics( InParticlePx, InParticlePy, InParticlePz, InParticleP, InParticleE, InParticleAngle, 1 );

        auto& OutParticlePx = Px[2000010000]; auto& OutParticlePy = Py[2000010000]; auto& OutParticlePz = Pz[2000010000];
        auto& OutParticleP = P[2000010000]; auto& OutParticleE = E[2000010000]; auto& OutParticleAngle = Angle[2000010000];
        ResizeKinematics( OutParticlePx, OutParticlePy, OutParticlePz, OutParticleP, OutParticleE, OutParticleAngle, 1 );

        auto& Ar40Px = Px[1000180400]; auto& Ar40Py = Py[1000180400]; auto& Ar40Pz = Pz[1000180400];
        auto& Ar40P = P[1000180400]; auto& Ar40E = E[1000180400]; auto& Ar40Angle = Angle[1000180400];
        ResizeKinematics( Ar40Px, Ar40Py, Ar40Pz, Ar40P, Ar40E, Ar40Angle, 1 );
        
        // Access event-wide variables directly for later use
        auto& EventPx = Px[0]; auto& EventPy = Py[0]; auto& EventPz = Pz[0];
        auto& EventP = P[0]; auto& EventE = E[0]; auto& EventAngle = Angle[0];

        auto& VisiblePx = Px[2000000400]; auto& VisiblePy = Py[2000000400]; auto& VisiblePz = Pz[2000000400];
        auto& VisibleP = P[2000000400]; auto& VisibleE = E[2000000400]; auto& VisibleAngle = Angle[2000000400];

        auto& VisibleNoNPx = Px[2000000401]; auto& VisibleNoNPy = Py[2000000401]; auto& VisibleNoNPz = Pz[2000000401];
        auto& VisibleNoNP = P[2000000401]; auto& VisibleNoNE = E[2000000401]; auto& VisibleNoNAngle = Angle[2000000401];

        auto& LeadingParticlePx = Px[2000000410]; auto& LeadingParticlePy = Py[2000000410]; auto& LeadingParticlePz = Pz[2000000410];
        auto& LeadingParticleP = P[2000000410]; auto& LeadingParticleE = E[2000000410]; auto& LeadingParticleAngle = Angle[2000000410];

        auto& LeadingParticleNoNPx = Px[2000000411]; auto& LeadingParticleNoNPy = Py[2000000411]; auto& LeadingParticleNoNPz = Pz[2000000411];
        auto& LeadingParticleNoNP = P[2000000411]; auto& LeadingParticleNoNE = E[2000000411]; auto& LeadingParticleNoNAngle = Angle[2000000411];

        auto& LeadingProtonPx = Px[2000000412]; auto& LeadingProtonPy = Py[2000000412]; auto& LeadingProtonPz = Pz[2000000412];
        auto& LeadingProtonP = P[2000000412]; auto& LeadingProtonE = E[2000000412]; auto& LeadingProtonAngle = Angle[2000000412];

        // Access particle varaibles
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
            // Check to see if event originates from a beam neutrino
            EventOrigin = MCTruthObj.Origin();
            //std::cout << "EventOrigin is " << EventOrigin << std::endl;
            if ( EventOrigin != simb::kBeamNeutrino ) {
              std::cout << "EventOrigin is not from beam neutrino! . . . "
                        << "MCTruth origin is " << EventOrigin << std::endl;
            }

            // The incident DM particle is stored in the neutrino container in the MCTruth
            simb::MCNeutrino const& InParticleObj = MCTruthObj.GetNeutrino();
            simb::MCParticle const& InParticle    = InParticleObj.Nu();

            // Skip the loop in the case that the incident particle is not a DM particle
            if ( isBDM && InParticle.PdgCode() != 2000010000 ) {
                std::cout << "The incident particle in the MCTruth " << iMCTruth << " is a neutrino, skip the interaction..." << std::endl;
                continue;
            }

            // Currently deal with only one DM-Argon interaction in an event
            if ( nInParticles > 0 ) {
                std::cout << "Found more than 1 DM-Argon interaction in the event!" << std::endl;
                continue;
            }

            // Interaction information
            Mode = InParticleObj.Mode();
            InteractionType = InParticleObj.InteractionType();
            CCNC = InParticleObj.CCNC();

            // Position of the vertex
            const auto Position = geo::vect::convertTo< Momentum4_t >( InParticle.Position() );
            VtxX = Position.X(); VtxY = Position.Y(); VtxZ = Position.Z(); T0 = Position.T();
            if ( abs( VtxX ) <= 745.744+1.e-12 && abs( VtxY ) <= 600.019+1.e-12 && VtxZ >= -0.87625-1.e-12 && VtxZ <= 5808.87+1.e-12 ) {
              isIn10kton = 1;
              ++nIn10kton;
            } else isIn10kton = 0;
            
            // 4-momentum of the incident DM particle
            const auto InParticleMom = geo::vect::convertTo< Momentum4_t >( InParticle.Momentum() );
            Event += InParticleMom;

            Momentum3_t InParticleMom3( InParticleMom.Px(), InParticleMom.Py(), InParticleMom.Pz() );
            InParticlePx[0] = InParticle.Px(); InParticlePy[0] = InParticle.Py(); InParticlePz[0] = InParticle.Pz();
            InParticleP[0] = InParticle.P(); InParticleE[0] = InParticle.E(); InParticleAngle[0] = 0.;
            if ( isBDM ) SunPosition = InParticleMom3;
            SunPx = SunPosition.X(); SunPy = SunPosition.Y(); SunPz = SunPosition.Z();

            // Now look for the initial argon 4-momentum and assign outgoing (final state) DM kinematic values
            for ( size_t iMCParticle = 0; iMCParticle < MCTruthObj.NParticles(); ++iMCParticle ) {

                const simb::MCParticle& thisMCParticle = MCTruthObj.GetParticle( iMCParticle );
                if ( thisMCParticle.StatusCode() == 1 && thisMCParticle.PdgCode() == 2000010000 ) {
                    const auto Momentum = geo::vect::convertTo< Momentum4_t >( thisMCParticle.Momentum() );
                    OutParticlePx[0] = Momentum.Px(); OutParticlePy[0] = Momentum.Py(); OutParticlePz[0] = Momentum.Pz();
                    OutParticleP[0] = Momentum.P(); OutParticleE[0] = Momentum.E();
                } else if ( !isBDM && thisMCParticle.StatusCode() == 1 && ( abs(thisMCParticle.PdgCode() ) == 12 ||abs(thisMCParticle.PdgCode() ) == 14 ) ) {
                    const auto Momentum = geo::vect::convertTo< Momentum4_t >( thisMCParticle.Momentum() );
                    OutParticlePx[0] = Momentum.Px(); OutParticlePy[0] = Momentum.Py(); OutParticlePz[0] = Momentum.Pz();
                    OutParticleP[0] = Momentum.P(); OutParticleE[0] = Momentum.E();
                } else if ( thisMCParticle.StatusCode() == 1 ) continue;
                // The incident DM particle has been filled
                if ( thisMCParticle.StatusCode() == 0 ) {
                    if ( thisMCParticle.PdgCode() == 2000010000 || abs( thisMCParticle.PdgCode() ) == 12 || abs( thisMCParticle.PdgCode() ) == 14 ) continue;
                    Event += geo::vect::convertTo< Momentum4_t >( thisMCParticle.Momentum() );
                    
                    // (Dane) adding in block to assign initial Ar40 Px, Py, Pz, E, etc.
                    Ar40Px[0] = thisMCParticle.Px(); Ar40Py[0] = thisMCParticle.Py(); Ar40Pz[0] = thisMCParticle.Pz();
                    Ar40P[0] = thisMCParticle.P(); Ar40E[0] = thisMCParticle.E();
                    // break;
                }
            }

            Momentum3_t OutParticleMom3( OutParticlePx[0], OutParticlePy[0], OutParticlePz[0] );
            OutParticleAngle[0] = CalculateAngle( SunPosition, OutParticleMom3 );

            
            // Find all the (stable final state) MCParticles associated to the MCTruth object
            std::vector< simb::MCParticle const* > const& G4MCParticles = G4MCParticlesAssn.at( iMCTruth );

            for ( size_t iMCParticle = 0; iMCParticle < G4MCParticles.size(); ++iMCParticle ) {

                const simb::MCParticle* thisMCParticle = G4MCParticles[iMCParticle];
                // Only consider primary particles
                if ( thisMCParticle->Process() != "primary" ) continue;
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

                    Momentum3_t GENIEMom3( GENIEPx[nGENIE], GENIEPy[nGENIE], GENIEPz[nGENIE] );
                    GENIEAngle[nGENIE] = CalculateAngle( SunPosition, GENIEMom3 );

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

                    Momentum3_t MesonMom3( MesonPx[nMesons], MesonPy[nMesons], MesonPz[nMesons] );
                    MesonAngle[nMesons] = CalculateAngle( SunPosition, MesonMom3 );

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

                    Momentum3_t BaryonMom3( BaryonPx[nBaryons], BaryonPy[nBaryons], BaryonPz[nBaryons] );
                    BaryonAngle[nBaryons] = CalculateAngle( SunPosition, BaryonMom3 );

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

                    Momentum3_t ParticleMom3( ParticlePx[nParticles], ParticlePy[nParticles], ParticlePz[nParticles] );
                    ParticleAngle[nParticles] = CalculateAngle( SunPosition, ParticleMom3 );

                    ++nParticles;
                    ++nAllParticles;

                    if ( pdgCode == 2000010000 || abs( pdgCode ) == 12 || abs( pdgCode ) == 14 ) continue;
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
            
        } // Loop over MCTruth

        // calculate the relevant angles between event-wide values: InDM and EventP, VisibleP, VisibleNoNP, etc.
        // unsure how to deal with ROOT::Math::DisplacementVector3D<Cartesian3D<Scalar> > objects right now so converting
        // to vector -- fix later. (Dane)
        Momentum3_t EventMom3( Event.Px(), Event.Py(), Event.Pz() );
        Momentum3_t VisibleMom3( Visible.Px(), Visible.Py(), Visible.Pz() );
        Momentum3_t VisibleNoNMom3( VisibleNoN.Px(), VisibleNoN.Py(), VisibleNoN.Pz() );
        Momentum3_t LeadingParticleMom3( LeadingParticle.Px(), LeadingParticle.Py(), LeadingParticle.Pz() );
        Momentum3_t LeadingParticleNoNMom3( LeadingParticleNoN.Px(), LeadingParticleNoN.Py(), LeadingParticleNoN.Pz() );
        Momentum3_t LeadingProtonMom3( LeadingProton.Px(), LeadingProton.Py(), LeadingProton.Pz() );

        ResizeKinematics( EventPx, EventPy, EventPz, EventP, EventE, EventAngle, 1 );
        ResizeKinematics( VisiblePx, VisiblePy, VisiblePz, VisibleP, VisibleE, VisibleAngle, 1 );
        ResizeKinematics( VisibleNoNPx, VisibleNoNPy, VisibleNoNPz, VisibleNoNP, VisibleNoNE, VisibleNoNAngle, 1 );
        ResizeKinematics( LeadingParticlePx, LeadingParticlePy, LeadingParticlePz, LeadingParticleP, LeadingParticleE, LeadingParticleAngle, 1 );
        ResizeKinematics( LeadingParticleNoNPx, LeadingParticleNoNPy, LeadingParticleNoNPz, LeadingParticleNoNP, LeadingParticleNoNE, LeadingParticleNoNAngle, 1 );
        ResizeKinematics( LeadingProtonPx, LeadingProtonPy, LeadingProtonPz, LeadingProtonP, LeadingProtonE, LeadingProtonAngle, 1 );

        EventAngle[0]              = CalculateAngle( SunPosition, EventMom3 );
        VisibleAngle[0]            = CalculateAngle( SunPosition, VisibleMom3 );
        VisibleNoNAngle[0]         = CalculateAngle( SunPosition, VisibleNoNMom3 );
        LeadingParticleAngle[0]    = CalculateAngle( SunPosition, LeadingParticleMom3 );
        LeadingParticleNoNAngle[0] = CalculateAngle( SunPosition, LeadingParticleNoNMom3 );
        LeadingProtonAngle[0]      = CalculateAngle( SunPosition, LeadingProtonMom3 );
        
        // Fill the event-wide information
        EventPx[0] = Event.Px();
        EventPy[0] = Event.Py();
        EventPz[0] = Event.Pz();
        EventP[0]  = Event.P();
        EventE[0]  = Event.E();

        VisiblePx[0] = Visible.Px();
        VisiblePy[0] = Visible.Py();
        VisiblePz[0] = Visible.Pz();
        VisibleP[0]  = Visible.P();
        VisibleE[0]  = Visible.E();

        VisibleNoNPx[0] = VisibleNoN.Px();
        VisibleNoNPy[0] = VisibleNoN.Py();
        VisibleNoNPz[0] = VisibleNoN.Pz();
        VisibleNoNP[0]  = VisibleNoN.P();
        VisibleNoNE[0]  = VisibleNoN.E();

        LeadingParticlePx[0] = LeadingParticle.Px();
        LeadingParticlePy[0] = LeadingParticle.Py();
        LeadingParticlePz[0] = LeadingParticle.Pz();
        LeadingParticleP[0] = LeadingParticle.P();
        LeadingParticleE[0] = LeadingParticle.E();

        LeadingParticleNoNPx[0] = LeadingParticleNoN.Px();
        LeadingParticleNoNPy[0] = LeadingParticleNoN.Py();
        LeadingParticleNoNPz[0] = LeadingParticleNoN.Pz();
        LeadingParticleNoNP[0] = LeadingParticleNoN.P();
        LeadingParticleNoNE[0] = LeadingParticleNoN.E();

        LeadingProtonPx[0] = LeadingProton.Px();
        LeadingProtonPy[0] = LeadingProton.Py();
        LeadingProtonPz[0] = LeadingProton.Pz();
        LeadingProtonP[0] = LeadingProton.P();
        LeadingProtonE[0] = LeadingProton.E();

        
        // ----------------------------------------------------------------------------------------------------
        // Smeared MC particles (cheated reconstruction)
        // ----------------------------------------------------------------------------------------------------
        
        auto& nSmearedVisible       = Multiplicity[-2000000400];
        auto& nSmearedVisibleNoN    = Multiplicity[-2000000401];
        auto& nSmearedRecontable    = Multiplicity[-2000000402];
        auto& nSmearedRecontableNoN = Multiplicity[-2000000403];
        Momentum4_t SmearedVisible, SmearedVisibleNoN, SmearedRecontable, SmearedRecontableNoN;
        Momentum4_t LeadingSmeared, LeadingSmearedNoN, LeadingSmearedRecontable, LeadingSmearedRecontableNoN;
        Momentum4_t LeadingSmearedProton, LeadingSmearedRecontableProton;
        double LeadingSmearedMax = 0., LeadingSmearedNoNMax = 0.;
        double LeadingSmearedRecontableMax = 0., LeadingSmearedRecontableNoNMax = 0.;
        double LeadingSmearedProtonMax = 0., LeadingSmearedRecontableProtonMax = 0.;
        
        for ( size_t iSmearedMC = 0; iSmearedMC < SmearedMCParts.size(); ++iSmearedMC ) {
          
            bdm::SmearedMCParticle const& SmearedMCPart = SmearedMCParts[iSmearedMC];
            int pdgCode = SmearedMCPart.particleId();
            int particleCode = -1*abs(pdgCode);

            auto& nParticles = Multiplicity[particleCode];
            auto& ParticlePx = Px[particleCode]; auto& ParticlePy = Py[particleCode];
            auto& ParticlePz = Pz[particleCode]; auto& ParticleP = P[particleCode];
            auto& ParticleE = E[particleCode]; auto& ParticleAngle = Angle[particleCode];

            if ( ParticlePx.size() <= nParticles )
                ResizeKinematics( ParticlePx, ParticlePy, ParticlePz, ParticleP, ParticleE, ParticleAngle, nParticles + 1 );

            ParticlePx[nParticles] = SmearedMCPart.momentum().X();
            ParticlePy[nParticles] = SmearedMCPart.momentum().Y();
            ParticlePz[nParticles] = SmearedMCPart.momentum().Z();
            ParticleP[nParticles] = SmearedMCPart.momentum().R();
            ParticleE[nParticles] = SmearedMCPart.E();
            Momentum4_t thisMomentum4( ParticlePx[nParticles], ParticlePy[nParticles], ParticlePz[nParticles], ParticleE[nParticles] );

            Momentum3_t ParticleMom3( ParticlePx[nParticles], ParticlePy[nParticles], ParticlePz[nParticles] );
            ParticleAngle[nParticles] = CalculateAngle( SunPosition, ParticleMom3 );

            ++nParticles;
            ++nSmearedVisible;
            SmearedVisible += thisMomentum4;
            if ( thisMomentum4.P() > LeadingSmearedMax ) {
                LeadingSmeared = thisMomentum4;
                LeadingSmearedMax = thisMomentum4.P();
            }
            
            if ( abs(pdgCode) == 2212 && thisMomentum4.P() > LeadingSmearedProtonMax ) {
                LeadingSmearedProton = thisMomentum4;
                LeadingSmearedProtonMax = thisMomentum4.P();
            }
            
            // Reconstructable
            if ( SmearedMCPart.isValid() ) {
                SmearedRecontable += thisMomentum4;
                ++nSmearedRecontable;
                if ( thisMomentum4.P() > LeadingSmearedRecontableMax ) {
                    LeadingSmearedRecontable = thisMomentum4;
                    LeadingSmearedRecontableMax = thisMomentum4.P();
                }
                if ( abs(pdgCode) == 2212 && thisMomentum4.P() > LeadingSmearedRecontableProtonMax ) {
                    LeadingSmearedRecontableProton = thisMomentum4;
                    LeadingSmearedRecontableProtonMax = thisMomentum4.P();
                }
            }
            
            // No neutrons
            if ( abs(pdgCode) != 2112 ) {
                SmearedVisibleNoN += thisMomentum4;
                ++nSmearedVisibleNoN;
                if ( thisMomentum4.P() > LeadingSmearedNoNMax ) {
                    LeadingSmearedNoN = thisMomentum4;
                    LeadingSmearedNoNMax = thisMomentum4.P();
                }
            }
            
            // No neutron and reconstrutable
            if ( abs(pdgCode) != 2112 && SmearedMCPart.isValid() ) {
                SmearedRecontableNoN += thisMomentum4;
                ++nSmearedRecontableNoN;
                if ( thisMomentum4.P() > LeadingSmearedRecontableNoNMax ) {
                    LeadingSmearedRecontableNoN = thisMomentum4;
                    LeadingSmearedRecontableNoNMax = thisMomentum4.P();
                }
            }
            
        } // Loop over smeared MC particles, std::vector< bdm::SmearedMCParticle >
        

        // Initialize the map
        auto& SmearedVisiblePx = Px[-2000000400]; auto& SmearedVisiblePy = Py[-2000000400]; auto& SmearedVisiblePz = Pz[-2000000400];
        auto& SmearedVisibleP = P[-2000000400]; auto& SmearedVisibleE = E[-2000000400]; auto& SmearedVisibleAngle = Angle[-2000000400];
        auto& SmearedVisibleNoNPx = Px[-2000000401]; auto& SmearedVisibleNoNPy = Py[-2000000401]; auto& SmearedVisibleNoNPz = Pz[-2000000401];
        auto& SmearedVisibleNoNP = P[-2000000401]; auto& SmearedVisibleNoNE = E[-2000000401]; auto& SmearedVisibleNoNAngle = Angle[-2000000401];
        auto& SmearedRecontablePx = Px[-2000000402]; auto& SmearedRecontablePy = Py[-2000000402]; auto& SmearedRecontablePz = Pz[-2000000402];
        auto& SmearedRecontableP = P[-2000000402]; auto& SmearedRecontableE = E[-2000000402]; auto& SmearedRecontableAngle = Angle[-2000000402];
        auto& SmearedRecontableNoNPx = Px[-2000000403]; auto& SmearedRecontableNoNPy = Py[-2000000403]; auto& SmearedRecontableNoNPz = Pz[-2000000403];
        auto& SmearedRecontableNoNP = P[-2000000403]; auto& SmearedRecontableNoNE = E[-2000000403]; auto& SmearedRecontableNoNAngle = Angle[-2000000403];
        auto& LeadingSmearedPx = Px[-2000000410]; auto& LeadingSmearedPy = Py[-2000000410]; auto& LeadingSmearedPz = Pz[-2000000410];
        auto& LeadingSmearedP = P[-2000000410]; auto& LeadingSmearedE = E[-2000000410]; auto& LeadingSmearedAngle = Angle[-2000000410];
        auto& LeadingSmearedNoNPx = Px[-2000000411]; auto& LeadingSmearedNoNPy = Py[-2000000411]; auto& LeadingSmearedNoNPz = Pz[-2000000411];
        auto& LeadingSmearedNoNP = P[-2000000411]; auto& LeadingSmearedNoNE = E[-2000000411]; auto& LeadingSmearedNoNAngle = Angle[-2000000411];
        auto& LeadingSmearedRecontablePx = Px[-2000000412]; auto& LeadingSmearedRecontablePy = Py[-2000000412]; auto& LeadingSmearedRecontablePz = Pz[-2000000412];
        auto& LeadingSmearedRecontableP = P[-2000000412]; auto& LeadingSmearedRecontableE = E[-2000000412]; auto& LeadingSmearedRecontableAngle = Angle[-2000000412];
        auto& LeadingSmearedRecontableNoNPx = Px[-2000000413]; auto& LeadingSmearedRecontableNoNPy = Py[-2000000413]; auto& LeadingSmearedRecontableNoNPz = Pz[-2000000413];
        auto& LeadingSmearedRecontableNoNP = P[-2000000413]; auto& LeadingSmearedRecontableNoNE = E[-2000000413]; auto& LeadingSmearedRecontableNoNAngle = Angle[-2000000413];
        auto& LeadingSmearedProtonPx = Px[-2000000414]; auto& LeadingSmearedProtonPy = Py[-2000000414]; auto& LeadingSmearedProtonPz = Pz[-2000000414];
        auto& LeadingSmearedProtonP = P[-2000000414]; auto& LeadingSmearedProtonE = E[-2000000414]; auto& LeadingSmearedProtonAngle = Angle[-2000000414];
        auto& LeadingSmearedRecontableProtonPx = Px[-2000000415]; auto& LeadingSmearedRecontableProtonPy = Py[-2000000415]; auto& LeadingSmearedRecontableProtonPz = Pz[-2000000415];
        auto& LeadingSmearedRecontableProtonP = P[-2000000415]; auto& LeadingSmearedRecontableProtonE = E[-2000000415]; auto& LeadingSmearedRecontableProtonAngle = Angle[-2000000415];
        
        
        ResizeKinematics( SmearedVisiblePx, SmearedVisiblePy, SmearedVisiblePz, SmearedVisibleP, SmearedVisibleE, SmearedVisibleAngle, 1 );
        ResizeKinematics( SmearedVisibleNoNPx, SmearedVisibleNoNPy, SmearedVisibleNoNPz, SmearedVisibleNoNP, SmearedVisibleNoNE, SmearedVisibleNoNAngle, 1 );
        ResizeKinematics( SmearedRecontablePx, SmearedRecontablePy, SmearedRecontablePz, SmearedRecontableP, SmearedRecontableE, SmearedRecontableAngle, 1 );
        ResizeKinematics( SmearedRecontableNoNPx, SmearedRecontableNoNPy, SmearedRecontableNoNPz, SmearedRecontableNoNP, SmearedRecontableNoNE, SmearedRecontableNoNAngle, 1 );
        ResizeKinematics( LeadingSmearedPx, LeadingSmearedPy, LeadingSmearedPz, LeadingSmearedP, LeadingSmearedE, LeadingSmearedAngle, 1 );
        ResizeKinematics( LeadingSmearedNoNPx, LeadingSmearedNoNPy, LeadingSmearedNoNPz, LeadingSmearedNoNP, LeadingSmearedNoNE, LeadingSmearedNoNAngle, 1 );
        ResizeKinematics( LeadingSmearedRecontablePx, LeadingSmearedRecontablePy, LeadingSmearedRecontablePz, LeadingSmearedRecontableP, LeadingSmearedRecontableE, LeadingSmearedRecontableAngle, 1 );
        ResizeKinematics( LeadingSmearedRecontableNoNPx, LeadingSmearedRecontableNoNPy, LeadingSmearedRecontableNoNPz, LeadingSmearedRecontableNoNP, LeadingSmearedRecontableNoNE, LeadingSmearedRecontableNoNAngle, 1 );
        ResizeKinematics( LeadingSmearedProtonPx, LeadingSmearedProtonPy, LeadingSmearedProtonPz, LeadingSmearedProtonP, LeadingSmearedProtonE, LeadingSmearedProtonAngle, 1 );
        ResizeKinematics( LeadingSmearedRecontableProtonPx, LeadingSmearedRecontableProtonPy, LeadingSmearedRecontableProtonPz, LeadingSmearedRecontableProtonP, LeadingSmearedRecontableProtonE, LeadingSmearedRecontableProtonAngle, 1 );

        
        // Calculate the angle
        Momentum3_t SmearedVisibleMom3( SmearedVisible.X(), SmearedVisible.Y(), SmearedVisible.Z() );
        Momentum3_t SmearedVisibleNoNMom3( SmearedVisibleNoN.X(), SmearedVisibleNoN.Y(), SmearedVisibleNoN.Z() );
        Momentum3_t SmearedRecontableMom3( SmearedRecontable.X(), SmearedRecontable.Y(), SmearedRecontable.Z() );
        Momentum3_t SmearedRecontableNoNMom3( SmearedRecontableNoN.X(), SmearedRecontableNoN.Y(), SmearedRecontableNoN.Z() );
        Momentum3_t LeadingSmearedMom3( LeadingSmeared.X(), LeadingSmeared.Y(), LeadingSmeared.Z() );
        Momentum3_t LeadingSmearedNoNMom3( LeadingSmearedNoN.X(), LeadingSmearedNoN.Y(), LeadingSmearedNoN.Z() );
        Momentum3_t LeadingSmearedRecontableMom3( LeadingSmearedRecontable.X(), LeadingSmearedRecontable.Y(), LeadingSmearedRecontableNoN.Z() );
        Momentum3_t LeadingSmearedRecontableNoNMom3( LeadingSmearedRecontableNoN.X(), LeadingSmearedRecontableNoN.Y(), LeadingSmearedRecontable.Z() );
        Momentum3_t LeadingSmearedProtonMom3( LeadingSmearedProton.X(), LeadingSmearedProton.Y(), LeadingSmearedProton.Z() );
        Momentum3_t LeadingSmearedRecontableProtonMom3( LeadingSmearedRecontableProton.X(), LeadingSmearedRecontableProton.Y(), LeadingSmearedRecontableProton.Z() );
        
        SmearedVisibleAngle[0] = CalculateAngle( SunPosition, SmearedVisibleMom3 );
        SmearedVisibleNoNAngle[0] = CalculateAngle( SunPosition, SmearedVisibleNoNMom3 );
        SmearedRecontableAngle[0] = CalculateAngle( SunPosition, SmearedRecontableMom3 );
        SmearedRecontableNoNAngle[0] = CalculateAngle( SunPosition, SmearedRecontableNoNMom3 );
        LeadingSmearedAngle[0] = CalculateAngle( SunPosition, LeadingSmearedMom3 );
        LeadingSmearedNoNAngle[0] = CalculateAngle( SunPosition, LeadingSmearedNoNMom3 );
        LeadingSmearedRecontableAngle[0] = CalculateAngle( SunPosition, LeadingSmearedRecontableMom3 );
        LeadingSmearedRecontableNoNAngle[0] = CalculateAngle( SunPosition, LeadingSmearedRecontableNoNMom3 );
        LeadingSmearedProtonAngle[0] = CalculateAngle( SunPosition, LeadingSmearedProtonMom3 );
        LeadingSmearedRecontableProtonAngle[0] = CalculateAngle( SunPosition, LeadingSmearedRecontableProtonMom3 );

        
        SmearedVisiblePx[0] = SmearedVisible.X();
        SmearedVisiblePy[0] = SmearedVisible.Y();
        SmearedVisiblePz[0] = SmearedVisible.Z();
        SmearedVisibleP[0] = SmearedVisible.P();
        SmearedVisibleE[0] = SmearedVisible.E();
        
        SmearedVisibleNoNPx[0] = SmearedVisibleNoN.X();
        SmearedVisibleNoNPy[0] = SmearedVisibleNoN.Y();
        SmearedVisibleNoNPz[0] = SmearedVisibleNoN.Z();
        SmearedVisibleNoNP[0] = SmearedVisibleNoN.P();
        SmearedVisibleNoNE[0] = SmearedVisibleNoN.E();
        
        SmearedRecontablePx[0] = SmearedRecontable.X();
        SmearedRecontablePy[0] = SmearedRecontable.Y();
        SmearedRecontablePz[0] = SmearedRecontable.Z();
        SmearedRecontableP[0] = SmearedRecontable.P();
        SmearedRecontableE[0] = SmearedRecontable.E();
        
        SmearedRecontableNoNPx[0] = SmearedRecontableNoN.X();
        SmearedRecontableNoNPy[0] = SmearedRecontableNoN.Y();
        SmearedRecontableNoNPz[0] = SmearedRecontableNoN.Z();
        SmearedRecontableNoNP[0] = SmearedRecontableNoN.P();
        SmearedRecontableNoNE[0] = SmearedRecontableNoN.E();
        
        LeadingSmearedPx[0] = LeadingSmeared.X();
        LeadingSmearedPy[0] = LeadingSmeared.Y();
        LeadingSmearedPz[0] = LeadingSmeared.Z();
        LeadingSmearedP[0] = LeadingSmeared.P();
        LeadingSmearedE[0] = LeadingSmeared.E();
        
        LeadingSmearedNoNPx[0] = LeadingSmearedNoN.X();
        LeadingSmearedNoNPy[0] = LeadingSmearedNoN.Y();
        LeadingSmearedNoNPz[0] = LeadingSmearedNoN.Z();
        LeadingSmearedNoNP[0] = LeadingSmearedNoN.P();
        LeadingSmearedNoNE[0] = LeadingSmearedNoN.E();
        
        LeadingSmearedRecontablePx[0] = LeadingSmearedRecontable.X();
        LeadingSmearedRecontablePy[0] = LeadingSmearedRecontable.Y();
        LeadingSmearedRecontablePz[0] = LeadingSmearedRecontable.Z();
        LeadingSmearedRecontableP[0] = LeadingSmearedRecontable.P();
        LeadingSmearedRecontableE[0] = LeadingSmearedRecontable.E();
        
        LeadingSmearedRecontableNoNPx[0] = LeadingSmearedRecontableNoN.X();
        LeadingSmearedRecontableNoNPy[0] = LeadingSmearedRecontableNoN.Y();
        LeadingSmearedRecontableNoNPz[0] = LeadingSmearedRecontableNoN.Z();
        LeadingSmearedRecontableNoNP[0] = LeadingSmearedRecontableNoN.P();
        LeadingSmearedRecontableNoNE[0] = LeadingSmearedRecontableNoN.E();
        
        LeadingSmearedProtonPx[0] = LeadingSmearedProton.X();
        LeadingSmearedProtonPy[0] = LeadingSmearedProton.Y();
        LeadingSmearedProtonPz[0] = LeadingSmearedProton.Z();
        LeadingSmearedProtonP[0] = LeadingSmearedProton.P();
        LeadingSmearedProtonE[0] = LeadingSmearedProton.E();
        
        LeadingSmearedRecontableProtonPx[0] = LeadingSmearedRecontableProton.X();
        LeadingSmearedRecontableProtonPy[0] = LeadingSmearedRecontableProton.Y();
        LeadingSmearedRecontableProtonPz[0] = LeadingSmearedRecontableProton.Z();
        LeadingSmearedRecontableProtonP[0] = LeadingSmearedRecontableProton.P();
        LeadingSmearedRecontableProtonE[0] = LeadingSmearedRecontableProton.E();
        
        
        
        
        fTree->Fill();
    } // End of an event

    fOut->Write();
    
    std::cout << "Number of events in 10k ton volume: " << nIn10kton << std::endl;

    return 0;
}
