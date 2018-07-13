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


#ifndef VDORMAND_PRINCE_745_HH
#define VDORMAND_PRINCE_745_HH

#include "DormandPrince745.hh"


template <typename EquationOfMotion>
class DormandPrince745<EquationOfMotion, Double_8v> {
public:
    DormandPrince745(EquationOfMotion* equation):
        fEquation(equation)
    {}
   
    void Stepper(
        const Double_8v& yInput,
        const Double_8v& dydxInput,
        G4double hstep,
        Double_8v& yOutput,
        Double_8v& yError);

    G4double DistChord() const;
    G4double DistChord2() const;
    G4double DistChord3() const;

    void Interpolate(G4double tau, Double_8v& y)
    {
        Interpolate4th(y, tau);
    }

private:
    struct ACoeffs {
        Double_8v fa[6];
        Double_8v& operator [] (size_t i) { return fa[i - 2]; }
        const Double_8v& operator [] (size_t i) const { return fa[i - 2]; }
    };

    void makeStep(
        const Double_8v& yInput,
        const Double_8v& dydx,
        G4double hstep,
        Double_8v& yOutput,
        Double_8v& yError,
        ACoeffs& ak) const;

    void Interpolate4th(Double_8v& y, G4double tau) const;

    EquationOfMotion* fEquation;
    Double_8v fyIn, fyOut, fdydx;
    ACoeffs fak;
    G4double fhstep;
};

#include "VDormandPrince745.icc"

#endif 
