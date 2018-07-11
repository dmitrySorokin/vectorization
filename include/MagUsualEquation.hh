#ifndef MAG_USUAL_EQUATION
#define MAG_USUAL_EQUATION


#include <G4ChargeState.hh>

template <typename Field, typename State>
class MagUsualEquation {
public:
    MagUsualEquation(Field* field);

    void SetChargeMomentumMass(G4ChargeState particleCharge,
                               G4double momentumXc,
                               G4double particleMass);

    void RightHandSide(const State& y, State& dydx);

private:
    G4double fCof;
    Field* fBField;
};

#include "MagUsualEquation.icc"

#endif