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
/// \file PS3SteppingAction.cc
/// \brief Implementation of the PS3SteppingAction class

#include "PS3SteppingAction.hh"
#include "PS3EventAction.hh"
#include "PS3DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"  // ADD THIS LINE

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3SteppingAction::PS3SteppingAction(PS3EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3SteppingAction::~PS3SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const PS3DetectorConstruction* detectorConstruction
      = static_cast<const PS3DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);  

  // Get track info and position
  G4Track* track = step->GetTrack();
  G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
  G4double depth = position.getZ()/cm + 20.0;

  // TEMPORARY DEBUG: Check large energy deposits
  //if (edepStep > 10.0*MeV) {
    //G4cout << "Large deposit: depth=" << depth 
         //  << " cm, E=" << edepStep/MeV 
          // << " MeV, particle=" << track->GetDefinition()->GetParticleName()
         //  << " parentID=" << track->GetParentID()
          // << " step length=" << step->GetStepLength()/mm << " mm"
          // << G4endl;
  //}

  // Record ALL energy deposition (realistic dose)
  if (edepStep > 0.001*MeV) {
    fEventAction->FillHistograms(edepStep, position.getZ()/cm, position.getX()/cm, position.getY()/cm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

