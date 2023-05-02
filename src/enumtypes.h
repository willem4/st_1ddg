#ifndef ENUMTYPES_H
#define ENUMTYPES_H

using namespace std;

//! Enum with Time Integration Methods.
enum timeintmethod {EulerForward, RungeKutta3, CrankNicolsonWithPredictor};    
//! Enum with Test Type Description. 
enum testtypedescr {Advection, Burgers, SWEBurgers, SWELinearWaveSolution, SWERiemannProblemLeftRarefactionRightShock, SWERiemannProblemLeftShockRightRarefaction,  SWERiemannProblemLeftShockRightShock,  SWERiemannProblemLeftRarefactionRightRarefaction, SWEContinuousTopography, SWEDiscontinuousTopography, SWEFlowOverIsolatedContinuousRidgeI, SWEFlowOverIsolatedDiscontinuousRidgeI, SWEFlowOverIsolatedContinuousRidgeIV, SWEFlowOverIsolatedDiscontinuousRidgeIV,BurgersDiffusive, SedimentTransport, HLLGrass, LFGrass, LFGrassMomentum, DecoupledGrassSedimentOnly, SuspendedSediment, LFGrassMomentumDiffusion, SWEFlowOverSinyBedI, SWERiemannProblemLeftRarefactionRightShockWidthChange};
//! Enum with Boundary Conditions.
enum boundarycondition {Periodic, Transmissive, Wall, Prescribed};

#endif
