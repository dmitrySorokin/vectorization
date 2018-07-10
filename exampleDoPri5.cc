#include "VDormandPrince745.hh"
#include "VMagUsualEquation.hh"

#include <G4DormandPrince745.hh>
#include <G4UniformMagField.hh>
#include <G4Mag_UsualEqRhs.hh>
#include <G4FieldTrack.hh>
#include <G4Proton.hh>
#include <G4DynamicParticle.hh>
#include <G4FieldUtils.hh>

#include <VecCore/Timer.h>

#include <memory>
#include <iostream>
#include <fstream>

const G4int INTEGRATED_COMPONENTS = 6;
const G4int NUMBER_OF_INTEGRATION_STEPS = 10000;
using State = G4double[G4FieldTrack::ncompSVEC];

template <typename Stepper, typename Equation>
void test(
    Stepper& method, 
    const Equation& equation, 
    const State& state,
    G4double stepLength,
    std::ofstream out)
{
    Timer<milliseconds> timer;

    State y, dydx, error;
    memcpy(y, state, sizeof(G4double) * G4FieldTrack::ncompSVEC);

    timer.Start();
    for (G4int i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        equation->RightHandSide(y, dydx);
        method.Stepper(y, dydx, stepLength, y, error);
        out << field_utils::makeVector(y, field_utils::Value3D::Position) << "\n";
        out << field_utils::makeVector(error, field_utils::Value3D::Position) << "\n";
        out << method.DistChord() << "\n";
    }
    G4double time = timer.Elapsed();


    printf("------------------------------------------\n");
    printf(">>  Time \n");
    printf(">>  %f\n", time);
    printf("------------------------------------------\n");
}

template <typename Stepper, typename Equation>
void vtest(
    Stepper& method, 
    const Equation& equation, 
    const State& state,
    G4double stepLength,
    std::ofstream out)
{
    Timer<milliseconds> timer;

    Double_8v y, dydx, error;
    memcpy(&y, state, sizeof(G4double) * G4FieldTrack::ncompSVEC);

    timer.Start();
    for (G4int i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        dydx = (*equation)(y);
        method.Stepper(y, dydx, stepLength, y, error);
        out << field_utils::makeVector(y, field_utils::Value3D::Position) << "\n";
        out << field_utils::makeVector(error, field_utils::Value3D::Position) << "\n";
        out << method.DistChord() << "\n";
    }
    G4double time = timer.Elapsed();


    printf("------------------------------------------\n");
    printf(">>  Time \n");
    printf(">>  %f\n", time);
    printf("------------------------------------------\n");
}

template<typename Equation, typename Field>
std::unique_ptr<Equation> makeEquation(Field* field, const G4DynamicParticle dynParticle)
{
    auto equation = std::make_unique<Equation>(field);
    equation->SetChargeMomentumMass(
        {
            dynParticle.GetCharge(),
            dynParticle.GetSpin(),
            dynParticle.GetMagneticMoment()
        },
        dynParticle.GetMomentum().mag(),
        dynParticle.GetMass()
    );

    return std::move(equation);
}

int main()
{
    G4DynamicParticle dynParticle(
            G4Proton::Definition(),
            G4ThreeVector(1, 0, 1).unit(),
            1 * CLHEP::MeV);

    auto track =
        std::make_shared<G4FieldTrack>(
            G4ThreeVector{0., 0., 0.}, // start position
            0,                         // LaboratoryTimeOfFlight
            dynParticle.GetMomentumDirection(),
            dynParticle.GetKineticEnergy(),
            dynParticle.GetMass(),
            dynParticle.GetCharge(),
            dynParticle.GetPolarization());

    auto field = std::make_unique<G4UniformMagField>(
            G4ThreeVector(0, 0, 1 * CLHEP::tesla));

    auto vequation = makeEquation<VMagUsualEquation>(field.get(), dynParticle);
    auto equation = makeEquation<G4Mag_UsualEqRhs>(field.get(), dynParticle);

    G4DormandPrince745 method(equation.get(), INTEGRATED_COMPONENTS);
    VDormandPrince745<VMagUsualEquation> vmethod(vequation.get());
    State y;
    track->DumpToArray(y);
    G4double stepLength = 2.5 * CLHEP::mm;

    vtest(vmethod, vequation, y, stepLength, std::ofstream("outv.txt"));
    test(method, equation, y, stepLength, std::ofstream("out.txt"));

    return 0;
}
