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

#include "aliases.hh"
#include <G4types.hh>

template <typename EquationOfMotion, typename State>
class DormandPrince745 {
public:
    DormandPrince745(EquationOfMotion* equation):
        fEquation(equation)
    {}
   
    void Stepper(
        const State& yInput,
        const State& dydxInput,
        G4double hstep,
        State& yOutput,
        State& yError);

    G4double DistChord() const;
    G4double DistChord2() const;
    G4double DistChord3() const;

    void Interpolate(G4double tau, State& y)
    {
        Interpolate4th(y, tau);
    }

private:
    struct ACoeffs {
        State fa[6];
        State& operator [] (size_t i) { return fa[i - 2]; }
        const State& operator [] (size_t i) const { return fa[i - 2]; }
    };

    void makeStep(
        const State& yInput,
        const State& dydx,
        G4double hstep,
        State& yOutput,
        State& yError,
        ACoeffs& ak) const;

    void Interpolate4th(State& y, G4double tau) const;

    EquationOfMotion* fEquation;
    State fyIn, fyOut, fdydx;
    ACoeffs fak;
    G4double fhstep;
};

#include "DormandPrince745.icc"

#endif 
