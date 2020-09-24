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
/// \file PS3RunAction.cc
/// \brief Implementation of the PS3RunAction class

#include "PS3RunAction.hh"
#include "PS3PrimaryGeneratorAction.hh"
#include "PS3DetectorConstruction.hh"
// #include "PS3Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3RunAction::PS3RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);

  fEdepLong = new std::vector<TH1F*>(0);
  fEdepRad = new std::vector<TH1F*>(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3RunAction::~PS3RunAction()
{
  for (auto h : fEdepLong)
    if (h) delete h;
  for (auto h : fEdepRad)
    if (h) delete h;
  delete fEdepLong;
  delete fEdepRad;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute average energy deposition and RMS
  //
  G4double edep  = fEdep.GetValue() / nofEvents;
  G4double edep2 = fEdep2.GetValue() / nofEvents;
  
  //G4double rms = edep2 - edep*edep/nofEvents;
  //if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  
  G4double rms = edep2 - edep*edep;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const PS3DetectorConstruction* detectorConstruction
   = static_cast<const PS3DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const PS3PrimaryGeneratorAction* generatorAction
   = static_cast<const PS3PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Average energy deposition per particle : " 
     << G4BestUnit(edep,"Energy") << " +/- " << G4BestUnit(rms,"Energy")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  TH1F* edepLongAverage = fEdepLong->at(0)->Clone("edepLongAverage");
  edepLongAverage->Reset();
  for (auto *h : fEdepLong)
    edepLongAverage->Add(h,1/(G4double)nofEvents);
  TH1F* edepRadAverage = fEdepRad->at(0)->Clone("edepRadAverage");
  edepRadAverage->Reset();
  for (auto *h : fEdepRad)
    edepRadAverage->Add(h,1/(G4double)nofEvents);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::StoreEdep(TH1F* edepLong, TH1F* edepRad) 
{
  size_t nEventsSaved = fEdepLong->size()
  fEdepLong->push_back(edepLong->Clone((std::to_string(nEventsSaved+1)+"_l").c_str()));
  fEdepRad->push_back(edepRad->Clone((std::to_string(nEventsSaved+1)+"_r").c_str()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

