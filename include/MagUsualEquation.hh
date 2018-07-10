#ifndef MAG_USUAL_EQUATION
#define MAG_USUAL_EQUATION

#include "aliases.hh"

#include <G4ChargeState.hh>

template <typename Field>
class MagUsualEquation {
public:
    MagUsualEquation(Field* field);

    void SetChargeMomentumMass(G4ChargeState particleCharge,
                               G4double momentumXc,
                               G4double particleMass);

    void RightHandSide(const G4double y[], G4double dydx[]);

private:
    G4double fCof;
    Field* fBField;
};

#include "MagUsualEquation.icc"

#endif