#include "VDormandPrince745.hh"

#include <VecCore/VecCore>

using namespace vecCore;
using Double_v = typename backend::VcSimdArray<8>::Double_v;

namespace {

void rightHandSide(
    G4EquationOfMotion* equation, 
    const Double_v& yIn, 
    Double_v& dydxOut)
{
    G4double y[8], dydx[8];
    Store(yIn, y);

    equation->RightHandSide(y, dydx);

    Load(dydxOut, dydx);
}

} // namepsace

VDormandPrince745::VDormandPrince745(
    G4EquationOfMotion* equation,
    G4int noIntegrationVariables)
    : G4MagIntegratorStepper(equation, noIntegrationVariables)
{
}

void VDormandPrince745::Stepper(
    const G4double yInput[],
    const G4double dydx[],
    G4double step,
    G4double yOutput[],
    G4double yError[])
{
    const G4double
        b21 = 0.2,
    
        b31 = 3.0 / 40.0, b32 = 9.0 / 40.0,
    
        b41 = 44.0/45.0, b42 = -56.0 / 15.0, b43 = 32.0 / 9.0,
        
        b51 = 19372.0 / 6561.0, b52 = -25360.0 / 2187.0, b53 = 64448.0 / 6561.0,
        b54 = -212.0 / 729.0 ,
        
        b61 = 9017.0 / 3168.0 , b62 =   -355.0 / 33.0,
        b63 =  46732.0 / 5247.0, b64 = 49.0 / 176.0,
        b65 = -5103.0 / 18656.0,
        
        b71 = 35.0 / 384.0, b72 = 0.,
        b73 = 500.0 / 1113.0, b74 = 125.0 / 192.0,
        b75 = -2187.0 / 6784.0, b76 = 11.0 / 84.0,
     
    
        dc1 = -(b71 - 5179.0 / 57600.0),
        dc2 = -(b72 - .0),
        dc3 = -(b73 - 7571.0 / 16695.0),
        dc4 = -(b74 - 393.0 / 640.0),
        dc5 = -(b75 + 92097.0 / 339200.0),
        dc6 = -(b76 - 187.0 / 2100.0),
        dc7 = -(- 1.0 / 40.0);


    Double_v yIn, dydxIn, yOut, yTemp, yErr;
    Double_v ak2, ak3, ak4, ak5, ak6, ak7;
    Load(yIn, yInput);
    Load(dydxIn, dydx);

    
    auto equation = GetEquationOfMotion();
    
    // RightHandSide(yIn, dydxIn);
    yTemp = yIn + b21 * step * dydxIn;

    rightHandSide(equation, yTemp, ak2);             
    yTemp = yIn + step * (b31 * dydxIn + b32 * ak2);

    rightHandSide(equation, yTemp, ak3);
    yTemp = yIn + step * (b41 * dydxIn + b42 * ak2 + b43 * ak3);

    rightHandSide(equation, yTemp, ak4);
    yTemp = yIn + step * (b51 * dydxIn + b52 * ak2 + b53 * ak3 + b54 * ak4);

    rightHandSide(equation, yTemp, ak5); 
    yTemp = yIn + step * 
        (b61 * dydxIn + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);

    rightHandSide(equation, yTemp, ak6);
    yOut = yIn + step * 
        (b71 * dydxIn + b72 * ak2 + b73 * ak3 + 
            b74 * ak4 + b75 * ak5 + b76 * ak6);

    rightHandSide(equation, yOut, ak7);
    yErr = step * 
        (dc1 * dydxIn + dc2 * ak2 + dc3 * ak3 + 
            dc4 * ak4 + dc5 * ak5 + dc6 * ak6 + dc7 * ak7);

    Store(yErr, yError);
    Store(yOut, yOutput);
}

G4double VDormandPrince745::DistChord() const
{
    assert(false);
    return 0;
}

