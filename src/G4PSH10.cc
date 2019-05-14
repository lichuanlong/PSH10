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

#include "G4PSH10.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

G4PSH10::G4PSH10(G4String name, G4int depth)
    : G4VPrimitiveScorer(name, depth), HCID(-1) {
    DefineUnitAndCategory();
    SetUnit("uSv");
    // conversion coefficients for photon H10_phi
    // datas from ICRP. 74. Table A.21.
    // units MeV, pSvcm2
    ConvCoeffP[0.010] = 0.061;
    ConvCoeffP[0.015] = 0.83;
    ConvCoeffP[0.020] = 1.05;
    ConvCoeffP[0.030] = 0.81;
    ConvCoeffP[0.040] = 0.64;
    ConvCoeffP[0.050] = 0.55;
    ConvCoeffP[0.060] = 0.51;
    ConvCoeffP[0.080] = 0.53;
    ConvCoeffP[0.100] = 0.61;
    ConvCoeffP[0.150] = 0.89;
    ConvCoeffP[0.200] = 1.20;
    ConvCoeffP[0.300] = 1.80;
    ConvCoeffP[0.400] = 2.38;
    ConvCoeffP[0.500] = 2.93;
    ConvCoeffP[0.600] = 3.44;
    ConvCoeffP[0.800] = 4.38;
    ConvCoeffP[1.] = 5.20;
    ConvCoeffP[1.5] = 6.90;
    ConvCoeffP[2.] = 8.60;
    ConvCoeffP[3.] = 11.1;
    ConvCoeffP[4.] = 13.4;
    ConvCoeffP[5.] = 15.5;
    ConvCoeffP[6.] = 17.6;
    ConvCoeffP[8.] = 21.6;
    ConvCoeffP[10.] = 25.6;
    // conversion coefficients for neutron H10_phi
    // datas from ICRP. 74. Table A.42.
    // units MeV, pSvcm2
    ConvCoeffN[1.00e-09] = 6.60;
    ConvCoeffN[1.00e-08] = 9.00;
    ConvCoeffN[2.53e-08] = 10.6;
    ConvCoeffN[1.00e-07] = 12.9;
    ConvCoeffN[2.00e-07] = 13.5;
    ConvCoeffN[5.00e-07] = 13.6;
    ConvCoeffN[1.00e-06] = 13.3;
    ConvCoeffN[2.00e-06] = 12.9;
    ConvCoeffN[5.00e-06] = 12.0;
    ConvCoeffN[1.00e-05] = 11.3;
    ConvCoeffN[2.00e-05] = 10.6;
    ConvCoeffN[5.00e-05] = 9.90;
    ConvCoeffN[1.00e-04] = 9.40;
    ConvCoeffN[2.00e-04] = 8.90;
    ConvCoeffN[5.00e-04] = 8.30;
    ConvCoeffN[1.00e-03] = 7.90;
    ConvCoeffN[2.00e-03] = 7.70;
    ConvCoeffN[5.00e-03] = 8.00;
    ConvCoeffN[1.00e-02] = 10.5;
    ConvCoeffN[2.00e-02] = 16.6;
    ConvCoeffN[3.00e-02] = 23.7;
    ConvCoeffN[5.00e-02] = 41.1;
    ConvCoeffN[7.00e-02] = 60.0;
    ConvCoeffN[1.00e-01] = 88.0;
    ConvCoeffN[1.50e-01] = 132.;
    ConvCoeffN[2.00e-01] = 170.;
    ConvCoeffN[3.00e-01] = 233.;
    ConvCoeffN[5.00e-01] = 322.;
    ConvCoeffN[7.00e-01] = 375.;
    ConvCoeffN[9.00e-01] = 400.;
    ConvCoeffN[1.00e+00] = 416.;
    ConvCoeffN[1.20e+00] = 425.;
    ConvCoeffN[2.00e+00] = 420.;
    ConvCoeffN[3.00e+00] = 412.;
    ConvCoeffN[4.00e+00] = 408.;
    ConvCoeffN[5.00e+00] = 405.;
    ConvCoeffN[6.00e+00] = 400.;
    ConvCoeffN[7.00e+00] = 405.;
    ConvCoeffN[8.00e+00] = 409.;
    ConvCoeffN[9.00e+00] = 420.;
    ConvCoeffN[1.00e+01] = 440.;
    ConvCoeffN[1.20e+01] = 480.;
    ConvCoeffN[1.40e+01] = 520.;
    ConvCoeffN[1.50e+01] = 540.;
    ConvCoeffN[1.60e+01] = 555.;
    ConvCoeffN[1.80e+01] = 570.;
    ConvCoeffN[2.00e+01] = 600.;
    ConvCoeffN[3.00e+01] = 515.;
    ConvCoeffN[5.00e+01] = 400.;
    ConvCoeffN[7.50e+01] = 330.;
    ConvCoeffN[1.00e+02] = 285.;
    ConvCoeffN[1.25e+02] = 260.;
    ConvCoeffN[1.50e+02] = 245.;
    ConvCoeffN[1.75e+02] = 250.;
    ConvCoeffN[2.01e+02] = 260.;

    LowerEdgeEnergyP = 0.01 * MeV;
    UpperEdgeEnergyP = 10. * MeV;
    LowerEdgeEnergyN = 1.00e-09 * MeV;
    UpperEdgeEnergyN = 2.01e+02 * MeV;
}

G4PSH10::~G4PSH10() { ; }

G4bool G4PSH10::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4double trklength = aStep->GetStepLength();
    G4String pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
    if (trklength == 0.) return FALSE;
    // G4cout << trklength << G4endl;
    G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
    if (pname == "gamma") {
        if (energy < LowerEdgeEnergyP || energy > UpperEdgeEnergyP) {
            G4ExceptionDescription ED;
            ED << "The energy of gamma " << energy / MeV << "MeV"
               << " is out of range [" << LowerEdgeEnergyP / MeV << ", "
               << UpperEdgeEnergyP / MeV << "]MeV." << G4endl;
            ED << "Do not calculate H*(10) for this step." << G4endl;
            G4Exception("G4PSH10::ProcessHits", "GammaEnergyOutOfRange",
                        JustWarning, ED);
            return FALSE;
        }
    } else if (pname == "neutron") {
        if (energy < LowerEdgeEnergyN || energy > UpperEdgeEnergyN) {
            G4ExceptionDescription ED;
            ED << "The energy of neutron " << energy / MeV << "MeV"
               << " is out of range [" << LowerEdgeEnergyN / MeV << ", "
               << UpperEdgeEnergyN / MeV << "]MeV." << G4endl;
            ED << "Do not calculate H*(10) for this step." << G4endl;
            G4Exception("G4PSH10::ProcessHits", "NeutronEnergyOutOfRange",
                        JustWarning, ED);
            return FALSE;
        }
    }
    G4int idx =
        ((G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()))
            ->GetReplicaNumber(indexDepth);
    G4double cubicVolume = ComputeVolume(aStep, idx);
    // G4cout << cubicVolume << G4endl;
    trklength *= aStep->GetPreStepPoint()->GetWeight();
    trklength /= cubicVolume;
    trklength *= CalcConvCoeff(energy, pname);
    G4int index = GetIndex(aStep);
    EvtMap->add(index, trklength);
    // G4cout << trklength << G4endl;
    return TRUE;
}

void G4PSH10::Initialize(G4HCofThisEvent* HCE) {
    EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if (HCID < 0) {
        HCID = GetCollectionID(0);
    }
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSH10::EndOfEvent(G4HCofThisEvent* /*HCE*/) { ; }

void G4PSH10::clear() { EvtMap->clear(); }

void G4PSH10::DrawAll() { ; }

void G4PSH10::PrintAll() { ; }

G4double G4PSH10::CalcConvCoeff(G4double E, G4String& pname) {
    if (pname != "gamma" && pname != "neutron") {
        return 0.;
    }
    G4double x = E / MeV;
    G4double H10_phi = 0.;

    // find the 4 points which should be used in lagrangian interpolation
    std::map<G4double, G4double>::iterator itrL, itrU, itr1, itr2;
    G4int dist = 0;
    if (pname == "gamma") {  // for photon
        itrL = ConvCoeffP.lower_bound(x);
        itrU = ConvCoeffP.upper_bound(x);
        dist = std::distance(itrU, ConvCoeffP.end());
        if (itrL == ConvCoeffP.begin()) {
            itr1 = itrL;
        } else if (dist == 0) {
            std::advance(itrU, -4);
            itr1 = itrU;
        } else if (dist == 1) {
            std::advance(itrU, -3);
            itr1 = itrU;
        } else if (dist == 2) {
            std::advance(itrU, -2);
            itr1 = itrU;
        } else {
            std::advance(itrL, -1);
            itr1 = itrL;
        }
        itr2 = itr1;
        std::advance(itr2, 4);
        std::map<G4double, G4double> fourPoints(itr1, itr2);

        // use the 4points map to interpolate
        // false : photon linear-log lagrangian interpolation
        H10_phi = CubicLagrange(x, fourPoints, false);

    } else {  // for neutron
        itrL = ConvCoeffN.lower_bound(x);
        itrU = ConvCoeffN.upper_bound(x);
        dist = std::distance(itrU, ConvCoeffN.end());
        if (itrL == ConvCoeffN.begin()) {
            itr1 = itrL;
        } else if (dist == 0) {
            std::advance(itrU, -4);
            itr1 = itrU;
        } else if (dist == 1) {
            std::advance(itrU, -3);
            itr1 = itrU;
        } else if (dist == 2) {
            std::advance(itrU, -2);
            itr1 = itrU;
        } else {
            std::advance(itrL, -1);
            itr1 = itrL;
        }
        itr2 = itr1;
        std::advance(itr2, 4);
        std::map<G4double, G4double> fourPoints(itr1, itr2);

        // use the 4points map to interpolate
        // true : neutron log-log lagrangian interpolation
        H10_phi = CubicLagrange(x, fourPoints, true);
    }

    // G4cout << pname << " " << x << " MeV, H10_phi " << H10_phi << G4endl;
    return H10_phi * gray / (1e+12) * cm * cm;
}

G4double G4PSH10::ComputeVolume(G4Step* aStep, G4int idx) {
    G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4VPVParameterisation* physParam = physVol->GetParameterisation();
    G4VSolid* solid = 0;
    if (physParam) {  // for parameterized volume
        if (idx < 0) {
            G4ExceptionDescription ED;
            ED << "Incorrect replica number --- GetReplicaNumber : " << idx
               << G4endl;
            G4Exception("G4PSCellFlux::ComputeVolume", "DetPS0001", JustWarning,
                        ED);
        }
        solid = physParam->ComputeSolid(idx, physVol);
        solid->ComputeDimensions(physParam, idx, physVol);
    } else {  // for ordinary volume
        solid = physVol->GetLogicalVolume()->GetSolid();
    }

    return solid->GetCubicVolume();
}

void G4PSH10::DefineUnitAndCategory() {
    // H10
    new G4UnitDefinition("Sievert", "Sv", "Dose", (gray));
    new G4UnitDefinition("milliSievert", "mSv", "Dose", (gray / 1000));
    new G4UnitDefinition("microSievert", "uSv", "Dose", (gray / 1000000));
}

void G4PSH10::SetUnit(const G4String& unit) {
    if (unit == "") {
        CheckAndSetUnit("uSv", "Dose");
    } else {
        CheckAndSetUnit(unit, "Dose");
    }
}

G4double G4PSH10::CubicLagrange(G4double E, std::map<G4double, G4double>& pmap,
                                G4bool elog) {
    if (pmap.size() != 4) {
        return 0.0;
    }
    G4double result = 0.;
    G4double up, down;
    std::map<G4double, G4double>::iterator itr1, itr2;
    for (itr1 = pmap.begin(); itr1 != pmap.end(); itr1++) {
        up = 1.;
        down = 1.;
        for (itr2 = pmap.begin(); itr2 != pmap.end(); itr2++) {
            if (itr2 != itr1) {
                if (elog) {
                    up *= (log10(E) - log10(itr2->first));
                    down *= (log10(itr1->first) - log10(itr2->first));
                } else {
                    up *= (E - itr2->first);
                    down *= (itr1->first - itr2->first);
                }
            }
        }
        result += ((up / down) * log10(itr1->second));
    }

    return pow(10., result);
}
