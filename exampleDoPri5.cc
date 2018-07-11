#include "aliases.hh"
#include "VDormandPrince745.hh"
#include "VMagUsualEquation.hh"
#include "DormandPrince745.hh"
#include "MagUsualEquation.hh"

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
const G4int NUMBER_OF_INTEGRATION_STEPS = 1000;
using G4State = G4double[G4FieldTrack::ncompSVEC];
using VEquation = MagUsualEquation<G4UniformMagField, Double_8v>;
using Equation = MagUsualEquation<G4UniformMagField, G4State>;


template <typename D, typename S>
void copy(D& dst, const S& src)
{
    for (G4int i = 0; i < G4FieldTrack::ncompSVEC; ++i) {
        dst[i] = src[i];
    }
}

template <>
void copy(Double_8v& dst, const G4State& src)
{
    vecCore::Load(dst, src);
}


template <typename Stepper, typename Equation, typename State>
void test(
    Stepper& method, 
    const Equation& equation, 
    const State& state,
    G4double stepLength,
    std::ofstream out)
{
    Timer<milliseconds> timer;

    State y, dydx, error;
    copy(y, state);

    timer.Start();
    for (G4int i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        equation->RightHandSide(y, dydx);
        method.Stepper(y, dydx, stepLength, y, error);
        //out << field_utils::makeVector(y, field_utils::Value3D::Position) << "\n";
        //out << field_utils::makeVector(error, field_utils::Value3D::Position) << "\n";
        //out << method.DistChord() << "\n";
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

    auto vequation = makeEquation<VEquation>(field.get(), dynParticle);
    auto iequation = makeEquation<Equation>(field.get(), dynParticle);
    auto equation = makeEquation<G4Mag_UsualEqRhs>(field.get(), dynParticle);

    G4DormandPrince745 method(equation.get(), INTEGRATED_COMPONENTS);
    DormandPrince745<VEquation, Double_8v> vmethod(vequation.get());
    DormandPrince745<Equation, G4State> imethod(iequation.get());
    
    G4State y;
    track->DumpToArray(y);

    Double_8v vy;
    copy(vy, y);
    
    G4double stepLength = 2.5 * CLHEP::mm;

    test(vmethod, vequation, vy, stepLength, std::ofstream("vout.txt"));
    test(imethod, iequation, y, stepLength, std::ofstream("iout.txt"));
    test(method, equation, y, stepLength, std::ofstream("out.txt"));

    return 0;
}
