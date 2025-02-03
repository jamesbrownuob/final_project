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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh" // For vascular transport model
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"

#include "CommandLineParser.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

#include "RunAction.hh"

#include "G4RandomTools.hh" // Include the random tools header

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes invisGrey(false, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes visRed(true, G4Colour(0, 0, 1));

// Define constants
const int numTubes = 20;
const G4double tubeRadius = 5.0 * micrometer; // 10 micrometers in diameter
const G4double tubeHeight = 3.0 * mm;
const G4double seedRadius = 0.15 * mm;
const G4double halfWaterBoxLength = 15.0 * mm; // Half the size of the water box

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  R = {155 * micrometer, 175 * micrometer, 195 * micrometer, 215 * micrometer, 235 * micrometer, 255 * micrometer, 275 * micrometer, 295 * micrometer, 315 * micrometer, 335 * micrometer};

  fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *S_Steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

  G4double nucleusSize = 300*nm;
  G4double margin = 0*nm;

  G4Box *solidWorld = new G4Box("world", 100 * mm, 100 * mm, 100 * mm);

  G4Box *solidWater = new G4Box("water", halfWaterBoxLength, halfWaterBoxLength, halfWaterBoxLength); 

  G4Tubs *solidSeed = new G4Tubs("seed", 0., seedRadius, 3 * mm, 0, 360 * degree);

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    air,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicWater = new G4LogicalVolume(solidWater,
                                                    waterMaterial,
                                                    "water");
  G4PVPlacement *physiWater = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                logicWater,
                                                "water",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicSeed = new G4LogicalVolume(solidSeed,
                                                   S_Steel,
                                                   "seed");

  G4PVPlacement *physiSeed = new G4PVPlacement(0,
                                               G4ThreeVector(),
                                               logicSeed,
                                               "seed",
                                               logicWater,
                                               0,
                                               false,
                                               0);

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.01*nucleusSize;
  auto fStepLimit = new G4UserLimits(maxStep);

  SetCells(Rmin, Rmax, Nrings);
  G4int numberRadialDivisions = NperRing;

  for (G4int r = 0; r < R.size(); ++r)
  {
        G4Tubs *solidCell = new G4Tubs("cell", R[r]-(nucleusSize/2+ margin), R[r]+(nucleusSize/2+ margin), 3 * mm, 0, 360 * degree);

        G4LogicalVolume *logicCell = new G4LogicalVolume(solidCell,
                                                         waterMaterial,
                                                         "cell");
        logicCell->SetUserLimits(fStepLimit);
        logicCell->SetVisAttributes(&visRed);

        G4PVPlacement *physiCell = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     logicCell,
                                                     "cell",
                                                     logicWater,
                                                     0,
                                                     r,
                                                     0);
      }
    
  // Create 20 cylindrical tubes with random paths through the box of water
  G4Material *bloodMaterial = man->FindOrBuildMaterial("G4_WATER"); // Example material for blood

  for (int i = 0; i < numTubes; ++i)
  {
      // Generate random start and end points on the faces of the cube
      G4ThreeVector startPoint, endPoint;
      do {
          // Randomly select a face for the start point
          int startFace = static_cast<int>(G4UniformRand() * 6);
          switch (startFace) {
              case 0: startPoint = G4ThreeVector(-halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 1: startPoint = G4ThreeVector(halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 2: startPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, -halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 3: startPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 4: startPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, -halfWaterBoxLength); break;
              case 5: startPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, halfWaterBoxLength); break;
          }

          // Randomly select a face for the end point
          int endFace = static_cast<int>(G4UniformRand() * 6);
          switch (endFace) {
              case 0: endPoint = G4ThreeVector(-halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 1: endPoint = G4ThreeVector(halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 2: endPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, -halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 3: endPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength); break;
              case 4: endPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, -halfWaterBoxLength); break;
              case 5: endPoint = G4ThreeVector(G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, G4UniformRand() * 2 * halfWaterBoxLength - halfWaterBoxLength, halfWaterBoxLength); break;
          }
      } while (startPoint.perp() <= seedRadius || endPoint.perp() <= seedRadius);

      // Calculate the direction and length of the tube
      G4ThreeVector direction = endPoint - startPoint;
      G4double length = direction.mag();
      G4ThreeVector midPoint = 0.5 * (startPoint + endPoint);
      direction /= length; // Normalize the direction

      // Create the cylindrical tube
      G4Tubs *solidTube = new G4Tubs("Tube", 0, tubeRadius, length / 2, 0, 360 * degree);

      G4LogicalVolume *logicTube = new G4LogicalVolume(solidTube, bloodMaterial, "Tube");

      // Place the tube at the mid-point with the calculated rotation
      G4RotationMatrix* rotMatrix = new G4RotationMatrix();
      rotMatrix->rotateX(atan2(direction.y(), direction.z()));
      rotMatrix->rotateY(atan2(direction.z(), direction.x()));
      rotMatrix->rotateZ(atan2(direction.y(), direction.x()));

      new G4PVPlacement(rotMatrix, midPoint, logicTube, "Tube", logicWater, false, i, false);
  }

  logicWorld->SetVisAttributes(&invisGrey);
  logicWater->SetVisAttributes(&invisGrey);
  logicSeed->SetVisAttributes(&visGrey);

  return physiWorld;
}

void DetectorConstruction::SetMin(G4double min)
{
  Rmin = min;
  RunAction *myRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());
  if (myRunAction) { // Added null check
    myRunAction->setRmin(min);
  }
}

void DetectorConstruction::SetMax(G4double max)
{
  Rmax = max;
  RunAction *myRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());
  if (myRunAction) { // Added null check
    myRunAction->setRmax(max);
  }
}

void DetectorConstruction::SetNrings(G4int N)
{
  Nrings = N;
  RunAction *myRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());
  if (myRunAction) { // Added null check
    myRunAction->setNrings(Nrings);
  }
}

void DetectorConstruction::SetNperRing(G4int N)
{
  NperRing = N;
  RunAction *myRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());
  myRunAction->setNperRing(NperRing);
}

void DetectorConstruction::SetCells(G4double min, G4double max, G4int Nrings)
{
  R.resize(Nrings);

  for (G4int i = 0; i < Nrings; ++i)
  {
    R[i] = min + i * (max - min) / (Nrings-1);
    G4cout << "R[" << i << "] = " << R[i] / um << G4endl;
  }
}