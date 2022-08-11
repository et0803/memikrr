#ifndef __muscleEffortEstimation_h__
#define __muscleEffortEstimation_h__

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Tools/AnalyzeTool.h>
#include <OpenSim/Analyses/StaticOptimization.h>
#include <OpenSim/Common/Array.h>
#include <OpenSim/Common/Storage.h>
#include <OpenSim/Common/Logger.h>
#include <OpenSim/Simulation/Model/ActivationFiberLengthMuscle.h>
#include <OpenSim/Common/Object.h>
#include <OpenSim/Actuators/RegisterTypes_osimActuators.h>
#include <string>
#include <Eigen/Geometry>

class muscleEffortEstimation: public OpenSim::AnalyzeTool{

public:
    ~muscleEffortEstimation();
    muscleEffortEstimation(const std::string &_staticOptimizationSetupXMLFile);
    double getMuscleEffortUsingOpensimAnalyzaTool(Eigen::VectorXd& dependentJCoordValues);
    void createStatesStorageFromCoordinatesInEigenVector(OpenSim::Model* workingModel, SimTK::State& s, Eigen::VectorXd& dependentJCoordValues);
    void runStaticOptimization();


private:
    OpenSim::Storage *statesStore;
    OpenSim::Storage *coordinatesStore;
    double *activiationLevels;
    int activationNum;  //number of muscle activations
};

#endif // __muscleEffortEstimation_h__
