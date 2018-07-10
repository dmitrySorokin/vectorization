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


#ifndef DORMAND_PRINCE_745_HH
#define DORMAND_PRINCE_745_HH

#include <G4FieldTrack.hh>

template <typename EquationOfMotion>
class DormandPrince745 {
public:
    DormandPrince745(EquationOfMotion* equation):
        fEquation(equation)
    {}
   
    void Stepper(
        const G4double yInput[],
        const G4double dydxInput[],
        G4double hstep,
        G4double yOutput[],
        G4double yError[]);

    G4double DistChord() const;

private:
    void makeStep(
        const G4double yInput[],
        const G4double dydx[],
        G4double hstep,
        G4double yOutput[],
        G4double yError[]);

    EquationOfMotion* fEquation;
    G4double fyIn[G4FieldTrack::ncompSVEC], fyOut[G4FieldTrack::ncompSVEC], fdydx[G4FieldTrack::ncompSVEC];
    G4double ak2[G4FieldTrack::ncompSVEC];
    G4double ak3[G4FieldTrack::ncompSVEC];
    G4double ak4[G4FieldTrack::ncompSVEC];
    G4double ak5[G4FieldTrack::ncompSVEC];
    G4double ak6[G4FieldTrack::ncompSVEC];
    G4double ak7[G4FieldTrack::ncompSVEC];
    G4double fhstep;
};

#include "DormandPrince745.icc"

#endif 
