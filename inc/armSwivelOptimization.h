#ifndef __armSwivelOptimization_h__
#define __armSwivelOptimization_h__

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>

#include"muscleEffortEstimation.h"
#include"hybridIK.h"

enum globalStopCriterionType
{
    stopCriterionOne=1,
    stopCriterionTwo=2
};


class armSwivelOptimization{

public:
    ~armSwivelOptimization();
    armSwivelOptimization();
    armSwivelOptimization(std::string& _staticOptimizationSetupXMLFile,   // static optimization setup file for muscle force and acitivation simulation.
                          std::string& _urdfFile,   // kinematic model file for feasible ik solution computation with hybrid IK.
                          std::string& _simplifiedSRSFile,  // kinematic model file for computing initials for the hybrid IK algorithm with a desired arm swivel angle.
                          double _optimizaitonInitialGuessInterval_InDegree,  // arm swivel angle sampling intervel to get geometrically equidistant initial guesses for the global optimization.
                          double _localStopCriterion_kappa,  
                          double _globalStopCriterion1_epsilon,
                          double _globalStopCriterion2_delta,
                          enum globalStopCriterionType _globalStopCriterionSelectionTag);
    void setDesiredWristPose(int poseIndex, Eigen::Isometry3d& _desiredWristPose);
    void appendLocalMinimusToFile(std::string& _outputFile);
    void runTwoPhaseOptimization();
    void localOptimization(Eigen::VectorXd _dependentJointCoordsAsIterationInitial, double modeledArmSwivelAngleForInitial);

private:
    Eigen::Isometry3d desiredWristPose;
    int currentPoseIndex;

    muscleEffortEstimation* armModelForMuscleEffortStaticOptimization; // with a musculoskeletal arm osim model
    hybridIK* IKForMusculoskeletalArmModel;

    double optimizaitonInitialGuessInterval_InDegree;
    double localStopCriterion_kappa;
    double globalStopCriterion1_epsiln;
    double globalStopCriterion2_delta;
    enum globalStopCriterionType globalStopCriterionSelectionTag;  //tag for global stop criterion selection

    //temp variable for local minimization
    Eigen::VectorXd delta_independentQ; // delta q for independent joints. All the coordinates in this function are in urdf definition.
    Eigen::VectorXd delta_dependentQNormalized; // normalized delta q for dependent joints
    Eigen::VectorXd iterationStep_dependentQ; // final iteration step;

    Eigen::MatrixXd jacobianMatrixForDependentQ;
    Eigen::MatrixXd jacobianMatrixForIndependentQ;
    Eigen::MatrixXd nullSpace_independentQ;

    Eigen::VectorXd forwardStep_dependentQ;  // q + rho*n in radian  in urdf definition
    Eigen::VectorXd backwardStep_dependentQ; // q - rho*n in radian  in urdf definition
    Eigen::VectorXd current_dependentQ;      // q in radian in urdf definition

    Eigen::VectorXd forwardStep_dependentQForOpensim;  // q + rho*n in radian in osim definition.
    Eigen::VectorXd backwardStep_dependentQForOpensim; // q - rho*n in radian in osim definition.
    Eigen::VectorXd current_dependentQForOpensim;      // q in radian in osim definition.

    std::vector<std::pair<Eigen::VectorXd,Eigen::Vector2d>> foundLocalMinimums; // each element of the vector contains a joint configuration and [minimum muscule activity cost, modeledArmSwivelAngleForInitial]
};

void transformIndependentCoordsFromURDF2OSIM(Eigen::VectorXd& input_independentCoordinatesInURDFDefinition, Eigen::VectorXd& output_independentCoordinatesInOSIMDefinition);  // transform the independent coordinates in urdf definition to indipendent coordinates in osim definition for muscle activiation simulation.

#endif // __armSwivelOptimization_h__
