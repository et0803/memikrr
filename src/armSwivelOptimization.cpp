#include"armSwivelOptimization.h"

armSwivelOptimization::~armSwivelOptimization()
{
    delete armModelForMuscleEffortStaticOptimization;
    delete IKForMusculoskeletalArmModel;
}
armSwivelOptimization::armSwivelOptimization()
{
    std::cout<<"can not establish the arm swivel optimization with proper parameters"<<std::endl;
}

armSwivelOptimization::armSwivelOptimization(std::string& _staticOptimizationSetupXMLFile,
                                             std::string& _urdfFile,
                                             std::string& _simplifiedSRSFile,
                                             double _optimizaitonInitialGuessInterval_InDegree,
                                             double _localStopCriterion_kappa,
                                             double _globalStopCriterion1_epsilon,
                                             double _globalStopCriterion2_delta,
                                             enum globalStopCriterionType _globalStopCriterionSelectionTag)
{
    armModelForMuscleEffortStaticOptimization = new muscleEffortEstimation(_staticOptimizationSetupXMLFile);
    IKForMusculoskeletalArmModel=new hybridIK(_urdfFile,_simplifiedSRSFile);

    optimizaitonInitialGuessInterval_InDegree = _optimizaitonInitialGuessInterval_InDegree;
    localStopCriterion_kappa = _localStopCriterion_kappa;
    globalStopCriterion1_epsiln = _globalStopCriterion1_epsilon;
    globalStopCriterion2_delta = _globalStopCriterion2_delta;
    globalStopCriterionSelectionTag = _globalStopCriterionSelectionTag;
}

void armSwivelOptimization::runTwoPhaseOptimization()
{
    // temporary variable declaration
    foundLocalMinimums.clear();
    unsigned int armSwivelIntervalDivisionScale=1;  //just for indicating the round of two-phase optimization
    bool globalSearchStopCriterionReached=false;

    Eigen::MatrixXd armSwivelIterationSeeds;                 //1 row, only store the arm swivel angle value
    Eigen::MatrixXd armSwivelIterationSeedsInLimitsTag;
    Eigen::MatrixXd dependentJointCoordValueOfArmSwivelIterationSeeds;  // the dependent joint coordinate value in the arm swivel configuration of the iteration seeds
    Eigen::MatrixXd armSwivelIterationSeedsLocalOptimizaitionFinishTag;
    Eigen::MatrixXd armSwivelIntervalDivisionKronMatrix; // a matrix for dividing the arm swivel angle sweep step by half with kron operation.

    armSwivelIntervalDivisionKronMatrix.resize(1,2);
    armSwivelIntervalDivisionKronMatrix.matrix()(0,0)=1;
    armSwivelIntervalDivisionKronMatrix.matrix()(0,1)=0;

    armSwivelIterationSeeds.resize(1,int(360.0/optimizaitonInitialGuessInterval_InDegree));
    for(unsigned int i=0; i< armSwivelIterationSeeds.cols(); i++)
        armSwivelIterationSeeds.matrix()(0,i) = i * optimizaitonInitialGuessInterval_InDegree/180.0*M_PI;
    armSwivelIterationSeedsInLimitsTag.resize(1,armSwivelIterationSeeds.cols());
    armSwivelIterationSeedsInLimitsTag.setZero();
    dependentJointCoordValueOfArmSwivelIterationSeeds.resize(IKForMusculoskeletalArmModel->armModel_urdf->robotChain.getNrOfJoints(),armSwivelIterationSeeds.cols());
    dependentJointCoordValueOfArmSwivelIterationSeeds.setZero();
    armSwivelIterationSeedsLocalOptimizaitionFinishTag.resize(1,armSwivelIterationSeeds.cols());
    armSwivelIterationSeedsLocalOptimizaitionFinishTag.setZero();

    // perform hybrid IK for arm swivel iteration seeds
    bool IKSolutionInLimitsTag;
    Eigen::MatrixXd feasibleIKSolutionsWithDesiredWristPoseAndArmSwivelAngle;

    for(unsigned int i=0; i< armSwivelIterationSeeds.cols(); i++){
        IKSolutionInLimitsTag = IKForMusculoskeletalArmModel->getIKSolutionForMusculoskeletalArmModel(armSwivelIterationSeeds.matrix()(0,i),
                                                                                                      desiredWristPose,
                                                                                                      feasibleIKSolutionsWithDesiredWristPoseAndArmSwivelAngle,
                                                                                                      currentPoseIndex); //return the feasible IK solution in joint limits
        if(IKSolutionInLimitsTag){
            armSwivelIterationSeedsInLimitsTag.matrix()(0,i)=1;
            dependentJointCoordValueOfArmSwivelIterationSeeds.col(i)=feasibleIKSolutionsWithDesiredWristPoseAndArmSwivelAngle.col(0); //just chose a feasible one if there are multiple ones.
            //std::cout<<currentPoseIndex<<" "<<armSwivelIterationSeeds.matrix()(0,i)*180.0/M_PI<<" "<<dependentJointCoordValueOfArmSwivelIterationSeeds.col(i).transpose()<<std::endl;

            //check forward kinematic error (wrist pose error)
            /*
            Eigen::VectorXd ikSolution= dependentJointCoordValueOfArmSwivelIterationSeeds.col(i);
            Eigen::Isometry3d wristPose = IKForMusculoskeletalArmModel->armModel_urdf->FK_getEndEffectorPoseInBaseFrame(ikSolution);
            Eigen::Matrix<ScalarType,6,1> twistForPoseError = poseDiff_fromFrame1ToFrame2InBase(desiredWristPose,wristPose);
            std::cout<< "pose error is (3 elements in m, 3 element in axis angle_radian): "<<twistForPoseError.transpose()<<" when arm swivel angle of "<<armSwivelIterationSeeds.matrix()(0,i)*180.0/M_PI<<" degrees for "<<currentPoseIndex<<" wrist pose"<<std::endl;
            */
        }
        else{
            //std::cout<<"no feasible IK solution for arm model with arm swivel angle of "<<armSwivelIterationSeeds.matrix()(0,i)*180.0/M_PI<<" degrees for "<<currentPoseIndex<<" wrist pose"<<std::endl;
        }
    }

    //global phase optimization
    while (!globalSearchStopCriterionReached) {
        for(unsigned int i=0; i < armSwivelIterationSeeds.cols(); i++)
            if(armSwivelIterationSeedsInLimitsTag.matrix()(0,i)>0 && armSwivelIterationSeedsLocalOptimizaitionFinishTag.matrix()(0,i)<1){
                localOptimization(dependentJointCoordValueOfArmSwivelIterationSeeds.col(i),armSwivelIterationSeeds.matrix()(0,i)/M_PI*180.0);  // local phase optimization
                armSwivelIterationSeedsLocalOptimizaitionFinishTag.matrix()(0,i)=1;  // set the tag indicating finishing local optimization
            }

        double stopCriterionValue;
        switch (globalStopCriterionSelectionTag)
        {
        case stopCriterionOne:
            stopCriterionValue = (double)foundLocalMinimums.size()*(armSwivelIterationSeedsInLimitsTag.sum()-1.0)/(double)(armSwivelIterationSeedsInLimitsTag.sum() - foundLocalMinimums.size()-2.0);
            if(stopCriterionValue<(foundLocalMinimums.size()+globalStopCriterion1_epsiln))
                globalSearchStopCriterionReached = true;
            break;
        case stopCriterionTwo:
            stopCriterionValue = (double)(foundLocalMinimums.size()*(foundLocalMinimums.size()+1.0))/(double)(armSwivelIterationSeedsInLimitsTag.sum()*(armSwivelIterationSeedsInLimitsTag.sum()-1.0));
            if(stopCriterionValue<globalStopCriterion2_delta)
                globalSearchStopCriterionReached = true;
            break;
        default:
            globalSearchStopCriterionReached = false;
            break;
        }
        if(globalSearchStopCriterionReached){
            break;
        }
        else{
            armSwivelIntervalDivisionScale = armSwivelIntervalDivisionScale * 2;

            armSwivelIterationSeeds = kron(armSwivelIterationSeeds,armSwivelIntervalDivisionKronMatrix);
            armSwivelIterationSeedsInLimitsTag = kron(armSwivelIterationSeedsInLimitsTag, armSwivelIntervalDivisionKronMatrix);
            armSwivelIterationSeedsLocalOptimizaitionFinishTag = kron(armSwivelIterationSeedsLocalOptimizaitionFinishTag,armSwivelIntervalDivisionKronMatrix);
            dependentJointCoordValueOfArmSwivelIterationSeeds = kron(dependentJointCoordValueOfArmSwivelIterationSeeds,armSwivelIntervalDivisionKronMatrix);

            for(unsigned int i=0; i<armSwivelIterationSeeds.size(); i++){
                if(i%2){
                    unsigned int lastSeedIndex=(i-1+armSwivelIterationSeeds.size())%armSwivelIterationSeeds.size();
                    unsigned int nextSeedIndex=(i+1+armSwivelIterationSeeds.size())%armSwivelIterationSeeds.size();
                    if(armSwivelIterationSeedsInLimitsTag.matrix()(0,lastSeedIndex)>0 && armSwivelIterationSeedsInLimitsTag.matrix()(0,nextSeedIndex)>0){
                        armSwivelIterationSeedsInLimitsTag.matrix()(0,i)=1;

                        if(abs(armSwivelIterationSeeds.matrix()(0,lastSeedIndex)-armSwivelIterationSeeds.matrix()(0,nextSeedIndex))>M_PI)
                            armSwivelIterationSeeds.matrix()(0,i)=(armSwivelIterationSeeds.matrix()(0,lastSeedIndex)+armSwivelIterationSeeds.matrix()(0,nextSeedIndex))/2.0 + M_PI;
                        IKSolutionInLimitsTag = IKForMusculoskeletalArmModel->getIKSolutionForMusculoskeletalArmModel(armSwivelIterationSeeds.matrix()(0,i),
                                                                                                                      desiredWristPose,
                                                                                                                      feasibleIKSolutionsWithDesiredWristPoseAndArmSwivelAngle,
                                                                                                                      currentPoseIndex);
                        if(IKSolutionInLimitsTag){
                            armSwivelIterationSeedsInLimitsTag.matrix()(0,i)=1.0;
                            dependentJointCoordValueOfArmSwivelIterationSeeds.col(i)=feasibleIKSolutionsWithDesiredWristPoseAndArmSwivelAngle.col(0); //just chose a feasible one if there are multiple ones.
                        }
                    }
                }
            }
        }
    }
}

void armSwivelOptimization::setDesiredWristPose(int _poseIndex, Eigen::Isometry3d& _desiredWristPose){
    desiredWristPose = _desiredWristPose;
    currentPoseIndex = _poseIndex;
}

void armSwivelOptimization::localOptimization(Eigen::VectorXd _dependentJointCoordsAsIterationInitial,       //21 joint coordinates in urdf definition
                                              double modeledArmSwivelAngleForInitial)
{
    double iterationStepSize_rho=0.05;  //about 1 degree in arm swivel angle
    double accumulateRho=0; // if accumulated rho is much more larger than the value cooresponding to the distance between two initial guesses, then this initial guess will not converge a local minimum more probably than the initial guess closer to it.

    bool inJointLimitsTag=true;

    double forwardMuscleEffort;   // for q + rho*n
    double backwardMuscleEffort;  // for q - rho*n
    double currentMuscleEffort;   // for q

    current_dependentQ = _dependentJointCoordsAsIterationInitial;
    transformIndependentCoordsFromURDF2OSIM(current_dependentQ,current_dependentQForOpensim);
    std::cout<<"coordinates for osim for a new local optimization starts."<<std::endl;
    std::cout<<currentPoseIndex<<" "<<modeledArmSwivelAngleForInitial<<" "<< current_dependentQForOpensim.transpose()<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;

    currentMuscleEffort = armModelForMuscleEffortStaticOptimization->getMuscleEffortUsingOpensimAnalyzaTool(current_dependentQForOpensim);

    while(iterationStepSize_rho>localStopCriterion_kappa  && accumulateRho < 0.15){  //less than 3 degrees far from the iteration seed when optimizaitonInitialGuessInterval_InDegree is 2.0 degrees.
        IKForMusculoskeletalArmModel->armModel_urdf->getJacobian_endEffectorInBaseFrame(current_dependentQ,jacobianMatrixForDependentQ);
        jacobianMatrixForIndependentQ = jacobianMatrixForDependentQ.matrix().block(0,jacobianMatrixForDependentQ.cols() - dependentCoordNum,jacobianMatrixForDependentQ.rows(),dependentCoordNum) *  IKForMusculoskeletalArmModel->W_fromIndependentCoordToDependentCoord;

        nullSpace_independentQ = jacobianMatrixForIndependentQ.fullPivLu().kernel();  //assume that jacobianMatrixForIndependentQ is full row rank
        delta_independentQ = nullSpace_independentQ.col(0);
        delta_dependentQNormalized = (IKForMusculoskeletalArmModel->W_fromIndependentCoordToDependentCoord * delta_independentQ).normalized();
        iterationStep_dependentQ = iterationStepSize_rho * delta_dependentQNormalized;

        forwardStep_dependentQ =  current_dependentQ;
        forwardStep_dependentQ.matrix().block(current_dependentQ.rows() - dependentCoordNum,0,dependentCoordNum,1) = forwardStep_dependentQ.matrix().block(current_dependentQ.rows() - dependentCoordNum,0,dependentCoordNum,1)+iterationStep_dependentQ;
        backwardStep_dependentQ = current_dependentQ;
        backwardStep_dependentQ.matrix().block(current_dependentQ.rows() - dependentCoordNum,0,dependentCoordNum,1) = backwardStep_dependentQ.matrix().block(current_dependentQ.rows() - dependentCoordNum,0,dependentCoordNum,1)-iterationStep_dependentQ;

        transformIndependentCoordsFromURDF2OSIM(forwardStep_dependentQ,forwardStep_dependentQForOpensim);
        transformIndependentCoordsFromURDF2OSIM(backwardStep_dependentQ,backwardStep_dependentQForOpensim);

        // call opensim simulation to determine update direction for joints value update

        // Caution! If the inDegrees=yes is checked in mot file,
        // in mot file for opensim 4.1 model, the unit of 2 coordinates deviation and flexion still need radian while the rest can be in degree. However, in mot file for opensim 3.3 gui application, all coordinate are in degrees.
        // This can be due to something wrong with opensim 4.1 gui application software. It is that all independent coordinate slider controls are in degree except the deviation and flexion are in radian.
        // Opensim 4.+ gui application can identify and tranform into gui slider controls the coordinates in mot file acoording to the inDegree flag expcet the two coordinates deviation and flexion which are always taken in radian.
        // Therefore, for C++ api to set the dependent coordinates, all coordinates can be in radian as implied in the osim file.
        // So, in mot file for opensim 3.3 and opensim 4.+, we can set inDegrees=no for simplification and consistency.

        forwardMuscleEffort = armModelForMuscleEffortStaticOptimization->getMuscleEffortUsingOpensimAnalyzaTool(forwardStep_dependentQForOpensim);
        backwardMuscleEffort = armModelForMuscleEffortStaticOptimization->getMuscleEffortUsingOpensimAnalyzaTool(backwardStep_dependentQForOpensim);

        accumulateRho=accumulateRho+iterationStepSize_rho;

        std::cout<<currentPoseIndex<<"(pose index): Current modeled arm swivel angle is "<<modeledArmSwivelAngleForInitial<<". Current rho is "<<iterationStepSize_rho<<". Accumulated rho is "<<accumulateRho<<". And squared total muscle-acitivation for three arm configuration are (back, current, forward): ("<<backwardMuscleEffort<<", "<<currentMuscleEffort<<", "<<forwardMuscleEffort<<")"<<std::endl;
        if((forwardMuscleEffort-currentMuscleEffort)*(backwardMuscleEffort-currentMuscleEffort)<0)  //E(q) is between E(q + d) and  E(q - d)
        {
            if(forwardMuscleEffort<currentMuscleEffort)
            {
                if(IKSolutionFeasibilityCheck(forwardStep_dependentQ,IKForMusculoskeletalArmModel->W_fromDependentCoordToInependentCoord,IKForMusculoskeletalArmModel->W_fromIndependentCoordToDependentCoord)){
                    current_dependentQ = forwardStep_dependentQ;
                    currentMuscleEffort = forwardMuscleEffort;
                }
                else
                    inJointLimitsTag=false;
            }
            else
            {
                if(IKSolutionFeasibilityCheck(backwardStep_dependentQ,IKForMusculoskeletalArmModel->W_fromDependentCoordToInependentCoord,IKForMusculoskeletalArmModel->W_fromIndependentCoordToDependentCoord)){
                    current_dependentQ = backwardStep_dependentQ;
                    currentMuscleEffort = backwardMuscleEffort;
                }
                else
                    inJointLimitsTag=false;
            }

            if(inJointLimitsTag){
                std::cout<<"Joint coordinates are in limits, and squared total muscle-acitivation for next current is: "<<currentMuscleEffort<<std::endl;
                std::cout<<std::endl;
                std::cout<<std::endl;
            }
            else{
                std::cout<<"current arm configuration is out of joint limits"<<std::endl;
                std::cout<<std::endl;
                std::cout<<std::endl;
                break;
            }
        }
        else{
            iterationStepSize_rho = iterationStepSize_rho/2;
            std::cout<<"The rho for next iteration is "<<iterationStepSize_rho<<", and squared total muscle-acitivation for next current is: "<<currentMuscleEffort<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl;
        }
    }
    std::cout<<"local iteration ended, the rho is "<<iterationStepSize_rho<<", and squared total muscle-acitivation for current is: "<<currentMuscleEffort<<". Number of found local minimums is "<<foundLocalMinimums.size()<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;

    // check whether the iteration result is a newly found local mininmum. If it is, then append it to the foundLocalMinimums;
    if(!inJointLimitsTag || ((forwardMuscleEffort>currentMuscleEffort)&&(backwardMuscleEffort>currentMuscleEffort))){
        bool addToNewFoundLocalMinimum=true;
        if (foundLocalMinimums.size()>0){
            for(unsigned int i=0; i<foundLocalMinimums.size();i++){
                if((foundLocalMinimums.at(i).first - current_dependentQ).norm()<1e-1){
                    addToNewFoundLocalMinimum=false;
                    break;
                }
            }
        }
        else
            addToNewFoundLocalMinimum=true;

        if(addToNewFoundLocalMinimum){
            std::pair<Eigen::VectorXd,Eigen::Vector2d> newlyFoundLocalMinimum;
            newlyFoundLocalMinimum.first = current_dependentQ;
            newlyFoundLocalMinimum.second[0] = currentMuscleEffort;
            newlyFoundLocalMinimum.second[1] = modeledArmSwivelAngleForInitial;
            foundLocalMinimums.push_back(newlyFoundLocalMinimum);
            std::cout<<currentPoseIndex<<"(pose index): squared total muscle-acitivation for current "<<currentMuscleEffort<<" is added to found local minimums"<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl;
        }
    }
}

void armSwivelOptimization::appendLocalMinimusToFile(std::string& _outputFile){
    std::ofstream outputFile;
    outputFile.open(_outputFile.c_str(),std::ios::app);  //append new data to the exist file or non-exist file

    if(outputFile.fail())
    {
        std::cout<<"write (append) result to file error!"<<std::endl;
        return;
    }

    std::ostringstream ss;
    ss.precision(std::numeric_limits<double>::digits10);
    for(unsigned int i=0;i<foundLocalMinimums.size();i++){  //
        //wrist pose index, modeled arm swivel angle for initial, muscle effort, arm joint coordinates
        ss<< currentPoseIndex<<" " << foundLocalMinimums.at(i).second[1]<<" "<< foundLocalMinimums.at(i).second[0]<<" "<<foundLocalMinimums.at(i).first.transpose() <<"\n";
    }
    outputFile << ss.str();
    outputFile.close();
}

void transformIndependentCoordsFromURDF2OSIM(Eigen::VectorXd& input_independentCoordinatesInURDFDefinition, Eigen::VectorXd& output_independentCoordinatesInOSIMDefinition){
    output_independentCoordinatesInOSIMDefinition = input_independentCoordinatesInURDFDefinition;

    int flexionCoordinateIndex=input_independentCoordinatesInURDFDefinition.rows() - dependentCoordNum + 17;

    output_independentCoordinatesInOSIMDefinition(flexionCoordinateIndex)=output_independentCoordinatesInOSIMDefinition(flexionCoordinateIndex)*2;
}
