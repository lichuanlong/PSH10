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

#ifndef G4PSH10_h
#define G4PSH10_h 1

#include <map>
#include "G4THitsMap.hh"
#include "G4VPrimitiveScorer.hh"

class G4PSH10 : public G4VPrimitiveScorer {
   public:  // with description
    G4PSH10(G4String name, G4int depth = 0);
    virtual ~G4PSH10();

   protected:  // with description
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    virtual G4double ComputeVolume(G4Step*, G4int idx);
    // copy from G4PSCellFlux
    virtual void DefineUnitAndCategory();

   public:
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);
    virtual void clear();
    virtual void DrawAll();
    virtual void PrintAll();
    G4double CalcConvCoeff(G4double E, G4String& pname);
    virtual void SetUnit(const G4String& unit);

   private:
    G4double CubicLagrange(G4double E, std::map<G4double, G4double>& pmap,
                           G4bool elog);

   private:
    G4int HCID;
    G4THitsMap<G4double>* EvtMap;
    std::map<G4double, G4double> ConvCoeffP;
    std::map<G4double, G4double> ConvCoeffN;
    G4double LowerEdgeEnergyP;
    G4double UpperEdgeEnergyP;
    G4double LowerEdgeEnergyN;
    G4double UpperEdgeEnergyN;
};
#endif
