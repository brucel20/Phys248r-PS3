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
#include "PS3Analysis.hh"
// #include "PS3Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "globals.hh"
#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3RunAction::PS3RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Register(fEdep);
  accumulableManager->Register(fEdep2);

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in PS3Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  analysisManager->SetVerboseLevel(1);
  if ( G4Threading::IsMultithreadedApplication() ) analysisManager->SetNtupleMerging(true);

  analysisManager->CreateNtuple("edep", "Total Event Energy Deposition");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PS3RunAction::~PS3RunAction()
{
  //delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  const PS3DetectorConstruction* detectorConstruction
      = static_cast<const PS3DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4LogicalVolume* scoringVolume = detectorConstruction->GetScoringVolume();
  G4String material = scoringVolume->GetMaterial()->GetName();
  G4cout << "RADIATION LENGTH: " << scoringVolume->GetMaterial()->GetRadlen() << G4endl;

  G4String fileName = "output_default.root";
  //analysisManager->OpenFile(fileName);

  if (!IsMaster()) {
    const PS3PrimaryGeneratorAction* generatorAction =
        static_cast<const PS3PrimaryGeneratorAction*>(
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    if (generatorAction) {
      const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
      G4String name = particleGun->GetParticleDefinition()->GetParticleName();
      G4double energy = particleGun->GetParticleEnergy();
      fileName = "ntuple_" + name + "_" + std::to_string((int)energy) + ".root";
    }
    analysisManager->OpenFile(fileName);

    //// Creating histograms
    analysisManager->CreateNtuple("showerEDep", "Energy Deposition in Volume");
    analysisManager->CreateNtupleDColumn("E");
    analysisManager->CreateNtupleDColumn("z");
    analysisManager->CreateNtupleDColumn("r");
    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    m_segment = 0.050*m;
    m_segment = 0.025*m;
    G4float offset = m_segment / 2;
    G4int nSegments = 1.9*m / m_segment;
    G4int histo2= analysisManager->CreateH2("EdepKTeV", "",
          nSegments, -0.95*m-offset, 0.95*m-offset,
          nSegments, -0.95*m-offset, 0.95*m-offset); // This doesn't appear to be filled?
    analysisManager->SetH2Activation(histo2, true);
    analysisManager->FinishNtuple();
  }
  // Creating histograms
  //analysisManager->OpenFile(fileName);
  //analysisManager->CreateNtuple("showerEDep", "Energy Deposition in Volume");
  //analysisManager->CreateNtupleDColumn("E");
  //analysisManager->CreateNtupleDColumn("z");
  //analysisManager->CreateNtupleDColumn("r");
  //analysisManager->CreateNtupleDColumn("x");
  //analysisManager->CreateNtupleDColumn("y");
  //std::cout << "here1" << std::endl;
  //m_segment = 0.050*m;
  //m_segment = 0.025*m;
  //G4float offset = m_segment / 2;
  //G4int nSegments = 1.9*m / m_segment;
  //G4int histo2= analysisManager->CreateH2("EdepKTeV", "",
  //      nSegments, -0.95*m-offset, 0.95*m-offset,
  //      nSegments, -0.95*m-offset, 0.95*m-offset); // This doesn't appear to be filled?
  //analysisManager->SetH2Activation(histo2, true);
  //analysisManager->FinishNtuple();

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
  // save histograms & ntuple
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PS3RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......