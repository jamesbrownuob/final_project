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
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include <memory>
#include "G4ThreeVector.hh"

#include "G4SingleParticleSource.hh"
#include "G4GeneralParticleSourceMessenger.hh"
#include "G4GeneralParticleSourceData.hh"
#include "G4ParticleGun.hh"

#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4String.hh"

class G4GeneralParticleSource;

class PrimaryGeneratorAction
    : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(G4String PS_data);
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event *event) override;
    G4String primaryName;
    G4int part1_CopyNum{-1};
    G4int part1_EventNum{-1};
    G4double part1_Time{-1};
    G4double primaryParticle{-1};
    G4int part1_particleSource{-1};
private:
    G4GeneralParticleSource *fpParticleGun;
    G4ParticleGun *fParticleGun;
    G4String fPS_data;

};