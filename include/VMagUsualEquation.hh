#ifndef VMAG_USUAL_EQUATION
#define VMAG_USUAL_EQUATION

#include "aliases.hh"
#include "MagUsualEquation.hh"

template <typename Field>
class MagUsualEquation<Field, Double_8v> {
public:
    MagUsualEquation(Field* field);

    void SetChargeMomentumMass(G4ChargeState particleCharge,
                               G4double momentumXc,
                               G4double particleMass);

    void RightHandSide(const Double_8v& y, Double_8v& dydx);

private:
    G4double fCof;
    Field* fBField;
};

#include "VMagUsualEquation.icc"

#endif