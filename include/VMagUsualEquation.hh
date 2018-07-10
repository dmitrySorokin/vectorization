#ifndef G4MAG_USUAL_EQUATION
#define G4MAG_USUAL_EQUATION

#include "aliases.hh"

#include <G4ChargeState.hh>

template <typename Field>
class VMagUsualEquation {
public:
    VMagUsualEquation(Field* field);

    void SetChargeMomentumMass(G4ChargeState particleCharge,
                               G4double momentumXc,
                               G4double particleMass);

    Double_8v operator() (const Double_8v& y);

private:
    G4double fCof;
    Field* fBField;
};

#include "VMagUsualEquation.icc"

#endif