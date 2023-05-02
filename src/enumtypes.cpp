using namespace std;

enum timeintmethod {EulerForward, RungeKutta3, CrankNicolsonWithPredictor};    
enum testtypedescr {Advection, Burgers, SWEBurgers, SWELinearWaveSolution, SWERiemannProblemLeftRarefactionRightShock, SWERiemannProblemLeftShockRightRarefaction,  SWERiemannProblemLeftShockRightShock,  SWERiemannProblemLeftRarefactionRightRarefaction, SWEContinuousTopography, SWEDiscontinuousTopography, SWEFlowOverIsolatedContinuousRidgeI, SWEFlowOverIsolatedDiscontinuousRidgeI, SWEFlowOverIsolatedContinuousRidgeIV, SWEFlowOverIsolatedDiscontinuousRidgeIV,BurgersDiffusive, SedimentTransport, HLLGrass, LFGrass, LFGrassMomentum, DecoupledGrassSedimentOnly, SuspendedSediment, LFGrassMomentumDiffusion, SWEFlowOverSinyBedI};
enum boundarycondition {Periodic, Transmissive, Wall, Prescribed};
