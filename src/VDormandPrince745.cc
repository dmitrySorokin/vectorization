#include "VDormandPrince745.hh"
#include "G4LineSection.hh"
#include <cmath>

#include "VecCore/Timer.h"
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

}// namepsace

//Constructor
VDormandPrince745::VDormandPrince745(
    G4EquationOfMotion* EqRhs,
    G4int noIntegrationVariables)
    : G4MagIntegratorStepper(EqRhs, noIntegrationVariables)
{
}


//	  The coefficients and the algorithm have been adapted from
//    Table 2 : Coefficients of RK5(4)7M
//	  ---Ref---
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//		Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//		------------------

//    The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
//
//    0   |
//    1/5 | 1/5
//    3/10| 3/40        9/40
//    4/5 | 44/45      −56/15      32/9
//    8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
//    1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
//    1   | 35/384      0         500/1113    125/192 −2187/6784    11/84
//    ------------------------------------------------------------------------
//          35/384       0        500/1113    125/192  −2187/6784    11/84   0
//          5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void VDormandPrince745::Stepper(
    const G4double yInput[],
    const G4double dydx[],
    G4double step,
    G4double yOutput[],
    G4double yError[])
{
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
        b21 = 0.2 ,
    
        b31 = 3.0 / 40.0, b32 = 9.0 / 40.0 ,
    
        b41 = 44.0/45.0, b42 = -56.0 / 15.0, b43 = 32.0 / 9.0,
        
        b51 = 19372.0 / 6561.0, b52 = -25360.0 / 2187.0, b53 = 64448.0 / 6561.0,
        b54 = -212.0 / 729.0 ,
        
        b61 = 9017.0 / 3168.0 , b62 =   -355.0 / 33.0,
        b63 =  46732.0 / 5247.0    , b64 = 49.0 / 176.0 ,
        b65 = -5103.0 / 18656.0 ,
        
        b71 = 35.0 / 384.0, b72 = 0.,
        b73 = 500.0 / 1113.0, b74 = 125.0 / 192.0,
        b75 = -2187.0 / 6784.0, b76 = 11.0 / 84.0,
     
    //Sum of columns, sum(bij) = ei
//    e1 = 0. ,
//    e2 = 1.0/5.0 ,
//    e3 = 3.0/10.0 ,
//    e4 = 4.0/5.0 ,
//    e5 = 8.0/9.0 ,
//    e6 = 1.0 ,
//    e7 = 1.0 ,
    
// Difference between the higher and the lower order method coeff. :
    // b7j are the coefficients of higher order
    
    dc1 = -(b71 - 5179.0 / 57600.0),
    dc2 = -(b72 - .0),
    dc3 = -(b73 - 7571.0 / 16695.0),
    dc4 = -(b74 - 393.0 / 640.0),
    dc5 = -(b75 + 92097.0 / 339200.0),
    dc6 = -(b76 - 187.0 / 2100.0),
    dc7 = -(- 1.0/40.0); //end of declaration


    Double_v yIn, dydxIn, yOut, yTemp, yErr;
    Double_v ak2, ak3, ak4, ak5, ak6, ak7;
    Load(yIn, yInput);
    Load(dydxIn, dydx);

    
    const G4int numberOfVariables = GetNumberOfVariables();
    auto equation = GetEquationOfMotion();
    
    
    // RightHandSide(yIn, DyDx);
    yTemp = yIn + b21 * step * dydxIn;

    rightHandSide(equation, yTemp, ak2);              // 2nd stage
    yTemp = yIn + step * (b31 * dydxIn + b32 * ak2);

    rightHandSide(equation, yTemp, ak3);              // 3rd stage
    yTemp = yIn + step * (b41 * dydxIn + b42 * ak2 + b43 * ak3);

    rightHandSide(equation, yTemp, ak4);              // 4th stage
    yTemp = yIn + step * (b51 * dydxIn + b52 * ak2 + b53 * ak3 + b54 * ak4);

    rightHandSide(equation, yTemp, ak5);              // 5th stage
    yTemp = yIn + step * 
        (b61 * dydxIn + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);

    rightHandSide(equation, yTemp, ak6);              // 6th stage
    yOut = yIn + step * 
        (b71 * dydxIn + b72 * ak2 + b73 * ak3 + 
            b74 * ak4 + b75 * ak5 + b76 * ak6);

    rightHandSide(equation, yOut, ak7);				//7th and Final stage
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

