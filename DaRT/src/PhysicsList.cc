//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 70268 2013-05-28 14:17:50Z maire $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
// #include "G4EmDNAPhysics_option2.hh"
// #include "G4EmDNAPhysicsActivator.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4EmPenelopePhysics.hh"
// #include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// #include "G4HadronElasticPhysicsHP.hh"
// #include "G4HadronPhysicsFTFP_BERT_HP.hh"
// #include "G4HadronPhysicsQGSP_BIC_HP.hh"
// #include "G4HadronInelasticQBBC.hh"
// #include "G4HadronPhysicsINCLXX.hh"
// #include "G4IonElasticPhysics.hh"
// #include "G4IonPhysics.hh"
// #include "G4IonINCLXXPhysics.hh"

// particles

// #include "G4BosonConstructor.hh"
// #include "G4LeptonConstructor.hh"
// #include "G4MesonConstructor.hh"
// #include "G4BosonConstructor.hh"
// #include "G4BaryonConstructor.hh"
// #include "G4IonConstructor.hh"
// #include "G4ShortLivedConstructor.hh"

#include "G4StepLimiterPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
:G4VModularPhysicsList()
{
  G4int verb = 1;
  SetVerboseLevel(verb);
  SetDefaultCutValue(1.0*nanometer);

  //add new units for radioActive decays
  //
//   new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
  // EM physics
  RegisterPhysics(new G4EmPenelopePhysics());

  G4EmParameters* param = G4EmParameters::Instance();
  param->SetAugerCascade(true);

  G4ProductionCutsTable::GetProductionCutsTable()->
    SetEnergyRange(100*eV, 1*GeV);

  // Decay
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  RegisterPhysics(new G4StepLimiterPhysics());
            
  // Hadron Elastic scattering
  // RegisterPhysics( new G4HadronElasticPhysicsHP(verb) );
  
  // // Hadron Inelastic physics
  // RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));
  
  // // Ion Elastic scattering
  // RegisterPhysics( new G4IonElasticPhysics(verb));
      
  // // Ion Inelastic physics
  // RegisterPhysics( new G4IonPhysics(verb));
    
  // // Gamma-Nuclear Physics
  // G4EmExtraPhysics* gnuc = new G4EmExtraPhysics(verb);
  // gnuc->ElectroNuclear(false);
  // gnuc->MuonNuclear(false);
  // RegisterPhysics(gnuc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void PhysicsList::ConstructParticle()
// {
//   G4BosonConstructor  pBosonConstructor;
//   pBosonConstructor.ConstructParticle();

//   G4LeptonConstructor pLeptonConstructor;
//   pLeptonConstructor.ConstructParticle();

//   G4MesonConstructor pMesonConstructor;
//   pMesonConstructor.ConstructParticle();

//   G4BaryonConstructor pBaryonConstructor;
//   pBaryonConstructor.ConstructParticle();

//   G4IonConstructor pIonConstructor;
//   pIonConstructor.ConstructParticle();

//   G4ShortLivedConstructor pShortLivedConstructor;
//   pShortLivedConstructor.ConstructParticle();  
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void PhysicsList::SetCuts()
// {
  
//   SetCutValue(10*um, "proton");
//   SetCutValue(10*um, "e-");
//   SetCutValue(10*um, "e+");
//   SetCutValue(10*um, "gamma");  
//   SetCutValue(1*nm, "Pb212");
//   SetCutValue(1*um, "alpha");



// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......