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

#include "G4H10Messenger.hh"
#include "G4PSH103D.hh"
#include "G4ScoringManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcommand.hh"
#include "G4VScoringMesh.hh"

G4H10Messenger::G4H10Messenger(G4ScoringManager* SManager) : fSMan(SManager) {
    QuantityCommands();
}

G4H10Messenger::~G4H10Messenger() {
    delete qH10Cmd;
}

void G4H10Messenger::SetNewValue(G4UIcommand* command, G4String newVal) {
    //
    // Get Current Mesh
    //
    G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
    if (!mesh) {
        G4cerr << "ERROR : No mesh is currently open. Open/create a mesh "
                  "first. Command ignored."
               << G4endl;
        return;
    }
    // Tokens
    G4TokenVec token;
    FillTokenVec(newVal, token);
    //
    if (command == qH10Cmd) {
        G4PSH103D* ps = new G4PSH103D(token[0]);
        ps->SetUnit(token[1]);
        mesh->SetPrimitiveScorer(ps);
    }
}

G4String G4H10Messenger::GetCurrentValue(G4UIcommand* /*command*/) {
    G4String val;

    return val;
}

void G4H10Messenger::FillTokenVec(G4String newValues, G4TokenVec& token) {
    G4Tokenizer next(newValues);
    G4String val;
    while (!(val = next()).isNull()) {
        token.push_back(val);
    }
}

void G4H10Messenger::QuantityCommands() {
    G4UIparameter* param;
    // Primitive Scorers
    qH10Cmd = new G4UIcommand("/score/quantity/H10", this);
    param = new G4UIparameter("qname", 's', false);
    qH10Cmd->SetParameter(param);
    param = new G4UIparameter("unit", 's', true);
    param->SetDefaultValue("uSv");
    qH10Cmd->SetParameter(param);
}
