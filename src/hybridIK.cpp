#include"hybridIK.h"

hybridIK::~hybridIK(){
    delete armModel_urdf;
    delete simplifiedArmModel_urdf;
}

hybridIK::hybridIK(){
    std::cout<<"can not build the hybridIK algorithm for the musculoskeletal arm model"<<std::endl;
}


hybridIK::hybridIK(const std::string& _urdfFile,
                   const std::string& _simplifiedSRSFile){
    //rootLink is base. endLink is wrist.
    armModel_urdf = new urdfModel(_urdfFile, "base", "wrist_hand");
    simplifiedArmModel_urdf = new urdfModel(_simplifiedSRSFile,"base","wrist");

    // in urdf, axis in <joint> section is defined in a joint frame according to the parameters in <origin>, not in the parent link frame.
    // http://wiki.ros.org/urdf/XML/joint   (including definition of mimic joints)

    //constant coordinates transformation matrix
    K_fromSRSInitialToArmIndependentCoord.resize(independentCoordNum,independentCoordNum);
    K_fromSRSInitialToArmIndependentCoord.setIdentity();
    K_fromSRSInitialToArmIndependentCoord.matrix()(2,0)=1;

    //post-multiply joint coordinates column vector to get the transformed one
    W_fromIndependentCoordToDependentCoord.resize(dependentCoordNum,independentCoordNum);
    W_fromIndependentCoordToDependentCoord.setZero();
    W_fromIndependentCoordToDependentCoord.matrix()(0,1)= -0.633555/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(1,1)= 0.322013/3.14159;
    W_fromIndependentCoordToDependentCoord.matrix()(2,1)=-0.322013/3.14159;
    W_fromIndependentCoordToDependentCoord.matrix()(3,1)=0.633555/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(4,1)=-0.128282/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(5,1)=1.03673/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(6,1)=0.466003/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(7,1)=-0.466003/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(8,1)=-1.03673/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(9,1)=0.128282/2.61799;
    W_fromIndependentCoordToDependentCoord.matrix()(10,0)=1;
    W_fromIndependentCoordToDependentCoord.matrix()(11,1)=1;
    W_fromIndependentCoordToDependentCoord.matrix()(12,0)=-1;
    W_fromIndependentCoordToDependentCoord.matrix()(13,2)=1;
    W_fromIndependentCoordToDependentCoord.matrix()(14,3)=1;
    W_fromIndependentCoordToDependentCoord.matrix()(15,4)=1;
    W_fromIndependentCoordToDependentCoord.matrix()(16,5)=0.5;               // this is for urdf definition, which is same for the following three dependent coordinates.
    W_fromIndependentCoordToDependentCoord.matrix()(17,6)=0.785398/1.5708;   //
    W_fromIndependentCoordToDependentCoord.matrix()(18,5)=0.5;               //
    W_fromIndependentCoordToDependentCoord.matrix()(19,6)=0.610865/1.22173;  //

    W_fromDependentCoordToInependentCoord.resize(independentCoordNum,dependentCoordNum); // here the independent coordinates are cooresponding to the IK solution for urdf definition.
    W_fromDependentCoordToInependentCoord.setZero();
    W_fromDependentCoordToInependentCoord.matrix()(0,10)=1;
    W_fromDependentCoordToInependentCoord.matrix()(1,11)=1;
    W_fromDependentCoordToInependentCoord.matrix()(2,13)=1;
    W_fromDependentCoordToInependentCoord.matrix()(3,14)=1;
    W_fromDependentCoordToInependentCoord.matrix()(4,15)=1;
    W_fromDependentCoordToInependentCoord.matrix()(5,16)=2;
    W_fromDependentCoordToInependentCoord.matrix()(6,17)=2;

    // SRS analytical IK solution parameters initialization

    /* sub link serials not work with KDL tree
    Eigen::VectorXd jointPosition_1, jointPosition_2;
    jointPosition_1.resize(2);
    jointPosition_1 <<  thoraxRotation_rzInDegree/180*M_PI, 0;
    jointPosition_2.resize(1);
    jointPosition_2 << 0;
    SRSShoulderFrameInBaseFrame = simplifiedArmModel_urdf->getPoseOfSpecificChildFrameInParentFrame("base","upperarm1",jointPosition_1) ;
    SRSUpperarmLength= simplifiedArmModel_urdf->getPoseOfSpecificChildFrameInParentFrame("upperarm3","ulna",jointPosition_2).translation().norm();
    SRSForearmLength= simplifiedArmModel_urdf->getPoseOfSpecificChildFrameInParentFrame("radius","proximal_row",jointPosition_2).translation().norm();
    */
    //r_z=thoraxRotation_rzInDegree
    SRSGravityDirectionInBaseFrame << 0,-1,0;
    SRSShoulderFrameInBaseFrame.matrix()(0, 0)=0.9928;
    SRSShoulderFrameInBaseFrame.matrix()(0, 1)=-0.1194;
    SRSShoulderFrameInBaseFrame.matrix()(0, 2)=0;
    SRSShoulderFrameInBaseFrame.matrix()(0, 3)=-0.0075;
    SRSShoulderFrameInBaseFrame.matrix()(1, 0)=0.1194;
    SRSShoulderFrameInBaseFrame.matrix()(1, 1)=0.9928;
    SRSShoulderFrameInBaseFrame.matrix()(1, 2)=0;
    SRSShoulderFrameInBaseFrame.matrix()(1, 3)=-0.0041;
    SRSShoulderFrameInBaseFrame.matrix()(2, 0)=0;
    SRSShoulderFrameInBaseFrame.matrix()(2, 1)=0;
    SRSShoulderFrameInBaseFrame.matrix()(2, 2)=1;
    SRSShoulderFrameInBaseFrame.matrix()(2, 3)=0.2164;
    SRSShoulderFrameInBaseFrame.matrix()(3, 0)=0;
    SRSShoulderFrameInBaseFrame.matrix()(3, 1)=0;
    SRSShoulderFrameInBaseFrame.matrix()(3, 2)=0;
    SRSShoulderFrameInBaseFrame.matrix()(3, 3)=1;
    SRSUpperarmLength=0.303365;
    SRSForearmLength=0.247625;
}

bool hybridIK::getIKSolutionForMusculoskeletalArmModel(double armSwivelAngle,
                                                       Eigen::Isometry3d& _desiredWristPoseInBase,
                                                       Eigen::MatrixXd& finalIKSolutions,
                                                       int ind){
    // the base frame of the high-fidelity arm urdf model, the simplified SRS arm urdf model and opensim model have been aligned.
    Eigen::MatrixXd SRSIKSolutionsForInitials;
    analyticalIKSolverForSRSArm(SRSShoulderFrameInBaseFrame,
                                _desiredWristPoseInBase,
                                SRSGravityDirectionInBaseFrame,
                                SRSUpperarmLength,
                                SRSForearmLength,
                                armSwivelAngle,
                                SRSIKSolutionsForInitials);  //each column in SRSIKSolutionsForInitials is a SRS IK solution
    //std::cout<<SRSIKSolutionsForInitials.transpose()<<std::endl;

    Eigen::MatrixXd IKSolutions;
    IKSolutions.resize(armModel_urdf->robotChain.getNrOfJoints(),SRSIKSolutionsForInitials.cols());  //each column in finalIKSolutions is a IK solution
    IKSolutions.setZero();
    unsigned int feasibleConfigurationNum=0;
    Eigen::VectorXd refinedIKSolution;
    Eigen::MatrixXd initialsForIndependentCoordinates;

    initialsForIndependentCoordinates = K_fromSRSInitialToArmIndependentCoord * SRSIKSolutionsForInitials; //each column in initialsForIndependentCoordinates is a IK solution
    enum t_iterationStatus KDLIKStatus;
    for(unsigned int i=0;i<initialsForIndependentCoordinates.cols();i++){
        KDLIKStatus = InitialStatusFlag;
        KDLIKStatus = KDLIKWithIndependentCoordInitial(armModel_urdf,
                                                       initialsForIndependentCoordinates.col(i),
                                                       _desiredWristPoseInBase,
                                                       W_fromIndependentCoordToDependentCoord,
                                                       refinedIKSolution);

        //std::cout<<"arm swivel angle before rectifying to [-pi, pi]: "<<armSwivelAngle<<"; "<<i <<"th solution: ("<<refinedIKSolution.transpose() << ") with KDL IK status code: "<< KDLIKStatus<<std::endl;
        int independentJointCoordinatesInLimitTag=0;
        if(IKSolutionFeasibilityCheck(refinedIKSolution,W_fromDependentCoordToInependentCoord, W_fromIndependentCoordToDependentCoord) && KDLIKStatus==E_NOERROR){
            IKSolutions.col(feasibleConfigurationNum)=refinedIKSolution;
            feasibleConfigurationNum++;
            independentJointCoordinatesInLimitTag=1;
        }

        //std::cout<<"arm swivel angle after  rectifying to [-pi, pi]: "<<armSwivelAngle/M_PI*180<<"; "<<i <<"th solution: ("<<refinedIKSolution.transpose() << ") with KDL IK status code: "<< KDLIKStatus<<" in limit:"<<independentJointCoordinatesInLimitTag<<std::endl;
        //std::cout<<"arm swivel angle after  rectifying to [-pi, pi], (arm swivel angle of "<<armSwivelAngle/M_PI*180<<" degrees), "<<i <<"th solution: ("<<(W_fromDependentCoordToInependentCoord*refinedIKSolution.block(1,0,20,1)).transpose() << ") with KDL IK status code: "<< KDLIKStatus<<" in limit:"<<independentJointCoordinatesInLimitTag<<std::endl;
    }
    if(feasibleConfigurationNum){
        //std::cout<<ind<<" "<<armSwivelAngle/M_PI*180<<" "<<refinedIKSolution.transpose() <<std::endl;
        //std::cout<<"final selected IK solution (arm swivel angle of "<<armSwivelAngle/M_PI*180<<" degrees): "<<(W_fromDependentCoordToInependentCoord*IKSolutions.matrix().block(1,0,20,1)).transpose() <<std::endl;
        finalIKSolutions=IKSolutions.matrix().block(0,0,armModel_urdf->robotChain.getNrOfJoints(),1);  //21 rows, each column in finalIKSolutions is a IK solution
        //std::cout<<"dependent coordinate: "<<finalIKSolutions.transpose()<<std::endl;
        return true;
    }
    else{
        finalIKSolutions.resize(0,0);
        return false;
    }
}

/*
model and algorithm definition:
1, this algorithm is designed for right arms. In the shoulder frame, x is to the front, y is to the up, z is from the left shoulder to the right shoulder (defined according to the opensim model)
2, in zero position, the frames at each of 7 joint-related links are defined same as the shoulder frame.
3, in zero position, the upper arm and forearm are towards -y direction (down along the gravity direction)
4, when the wrist is in front of the body and the elbow is at the right side of the vertical swivel reference plane, the arm swivel angle is positive.
5, all angle are in radian
6, there are three sigularities (shoulder singularity, wrist singularity, elbow zero position singularity)
7, rotation axis: 1(y), 2(-x), 3(y), 4(z), 5(y), 6(-z), 7(x)
8, when elbow singularity occurs only, joint5 is assume to be at zero position and joint3 serves the arm in-line rotation at first. And the further partition between joint3 and joint5 can be balanced according to jacobian requirements.
9, when elbow and wrist singularity occurs only, joint5 and joint7 assume to be at zero position and joint3 serves the arm in-line rotation first. And the further partition between joint3, joint5 and joint7 can be balanced according to jacobian requirements.
10, when elbow and shoulder singularity occurs only, joint1 and joint5 assume to be at zero position and joint3 serves the arm in-line rotation first. And the furthur partition between joint1, joint3 and joint5 can be balanced according to jacobian requirements.
11, when elbow, wrist and shoulder singularity occurs simultaneously, joint1, joint5, joint7 assume to be at zero position and joint3 serves the arm in-line rotation first. And the further partition between joint1, joint3, joint5 and joint7 can be balanced according to jacobian requirements.
12, when elbow singularity is absent and shoulder singularity occurs, joint1 is assume to be at zero position and only joint3 serves the  upper arm in-line rotation. And the further partition between joint1 and joint3 can be balanced according to jacobian requirements.
13, when elbow singularity is absent and wrist singularity occurs, joint7 is assume to be at zero position and only joint5 serves the forearm in-line rotation. And the further partition between joint5 and joint7 can be balanced according to jacobian requirements.

additional comments:
1, S: shoulder point, E: elbow point, W:wrist point,
2, shoulderBaseFrame and wristEndFrame are represented in a world reference frame.
3, gravity direction is down towards ground
4, gravity direciton starting from the shoulder point and the wrist point form the reference plane of the arm swivel angle.
5, each row of IKSolutionSet contains 7 elements as a IK solution.
6, all frame and direction data are transformed into the shoulder frame for further computation
7, any rotation matrix can be resolved into any type of euler angle (https://www.wikiwand.com/en/Euler_angles)
*/
bool analyticalIKSolverForSRSArm(Eigen::Isometry3d& shoulderBaseFrame,
                                 Eigen::Isometry3d& wristEndFrame,
                                 Eigen::Vector3d& gravityDirectionInWorldReferenceFrame,
                                 double upperArmLength,
                                 double foreArmLength,
                                 double armSwivelAngle,
                                 Eigen::MatrixXd& IKSolutionSet)
{
    //calculate the angle of the revolution joint at elbow in the range of [0,M_PI], otherwise [-M_PI, M_PI)
    double cosValue_Angle_ES_EW;
    Eigen::Isometry3d wristFrameInShoulderBaseFrame = shoulderBaseFrame.inverse()*wristEndFrame;  //shoulderBaseFrame and wristEndFrame are represented in a world reference frame.
    cosValue_Angle_ES_EW = (upperArmLength*upperArmLength + foreArmLength*foreArmLength - pow(wristFrameInShoulderBaseFrame.translation().norm(),2))/(2*upperArmLength*foreArmLength);
    double theta4;
    if(cosValue_Angle_ES_EW>1)
        theta4= M_PI;
    else
        if(cosValue_Angle_ES_EW<-1)
            theta4=0;
        else
            theta4=M_PI-acos(cosValue_Angle_ES_EW);

    //calculate the position of elbow
    double cosValue_Angle_SW_SE,angle_SW_SE;
    cosValue_Angle_SW_SE = (upperArmLength*upperArmLength + pow(wristFrameInShoulderBaseFrame.translation().norm(),2) - foreArmLength*foreArmLength)/(2*upperArmLength*wristFrameInShoulderBaseFrame.translation().norm());
    if(cosValue_Angle_SW_SE>1)
        angle_SW_SE=0;
    else
        if(cosValue_Angle_SW_SE<-1)
            angle_SW_SE=M_PI;
        else
            angle_SW_SE = acos(cosValue_Angle_SW_SE);

    Eigen::Vector3d gravityDirectionInShoulderFrame = shoulderBaseFrame.inverse().rotation()*gravityDirectionInWorldReferenceFrame;
    //Eigen::Vector3d gravityDirectionInShoulderFrame = gravityDirectionInWorldReferenceFrame;
    Eigen::Vector3d negativeArmSwivelDirectionFromReferencePlaneInShoulderFrame = (wristFrameInShoulderBaseFrame.translation().cross(gravityDirectionInShoulderFrame)).normalized();
    Eigen::Vector3d armSwivelRotationAxisInShoulderFrame = -wristFrameInShoulderBaseFrame.translation().normalized();

    Eigen::AngleAxisd rotationForAngle_SW_SE(angle_SW_SE,negativeArmSwivelDirectionFromReferencePlaneInShoulderFrame);
    Eigen::AngleAxisd rotationForArmSwivel(armSwivelAngle,armSwivelRotationAxisInShoulderFrame);

    Eigen::Vector3d elbowPositionWithDesiredArmSwivelAngleInShoulderFrame;
    elbowPositionWithDesiredArmSwivelAngleInShoulderFrame = rotationForArmSwivel*rotationForAngle_SW_SE*wristFrameInShoulderBaseFrame.translation().normalized()*upperArmLength;

    //compute the frames of the upper arm and forearm in consideration of elbow and wrist singularity
    Eigen::Matrix3d rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame;
    rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(1) = - elbowPositionWithDesiredArmSwivelAngleInShoulderFrame.normalized(); //y
    Eigen::Vector3d unNormalizedZ = (wristFrameInShoulderBaseFrame.translation()-elbowPositionWithDesiredArmSwivelAngleInShoulderFrame).cross(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(1));  //possible invalid due to elbow singularity
    rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(2) = unNormalizedZ.normalized(); //z, possible invalid due to elbow singularity
    rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(0) = ((rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(1)).cross(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(2))).normalized(); //x

    Eigen::Matrix3d rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame;
    rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(1) = (elbowPositionWithDesiredArmSwivelAngleInShoulderFrame-wristFrameInShoulderBaseFrame.translation()).normalized(); //y
    rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(2) = rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(2);  //z, possible invalid due to elbow singularity
    rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(0) = ((rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(1)).cross(rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(2))).normalized(); //x

    Eigen::Matrix3d rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4;

    if(unNormalizedZ.norm() < 1e-10){  //joint4 is at zero. elbow singularity is due to colinearlity in y axis of upper arm and forearm

        Eigen::Vector3d yOfForearmFrameAtElbowAfterJoint4InWristEndFrame = wristFrameInShoulderBaseFrame.rotation().transpose() * rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(1); //forearm and upper arm are colinear.

        // forearm frame at elbow after joint4 represented in wrist end frame can be formulated by the euler rotation series of xzy [from wrist end frame to forearm frame at elbow after joint4]
        double XZYEulerAngle_beta = asin(-yOfForearmFrameAtElbowAfterJoint4InWristEndFrame(0));  // alpha, beta and gama represent rotation angles corresponding to the first, second and third rotation, namely, around x, z and y in this euler angle type.
        Eigen::Matrix3d rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame;

        if(fabs(fabs(yOfForearmFrameAtElbowAfterJoint4InWristEndFrame(0))-1) < 1e-10){  //wrist singularity and the constraint is: joint3, joint5 and joint7 sum to a specific total rotation angle.
            // in this case, assume that joint5 and joint7 are at zero position. And the colinear rotation is only caused by joint3 (also joint1 when shoulder singularity occurs)
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(0)=Eigen::Vector3d(cos(XZYEulerAngle_beta),sin(XZYEulerAngle_beta),0);
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(1)=yOfForearmFrameAtElbowAfterJoint4InWristEndFrame;
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(2)=Eigen::Vector3d(0,0,1);
        }
        else{ //wrist is not in singular configuration, joint6 is not at +-90 degree.
            double XZYEulerAngle_alpha = aCosSin(yOfForearmFrameAtElbowAfterJoint4InWristEndFrame(1)/cos(XZYEulerAngle_beta),
                                                 yOfForearmFrameAtElbowAfterJoint4InWristEndFrame(2)/cos(XZYEulerAngle_beta));
            // in this case, assume that joint5 is at zero position. And the colinear rotation is only caused by joint3 (also joint1 when shoulder singularity occurs)
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(0) = Eigen::Vector3d(cos(XZYEulerAngle_beta),cos(XZYEulerAngle_alpha)*sin(XZYEulerAngle_beta), sin(XZYEulerAngle_alpha)*sin(XZYEulerAngle_beta));
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(1) = yOfForearmFrameAtElbowAfterJoint4InWristEndFrame;
            rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.col(2) = Eigen::Vector3d(0,-sin(XZYEulerAngle_alpha), cos(XZYEulerAngle_alpha));
        }

        rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame = wristFrameInShoulderBaseFrame.rotation()*rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame;
        rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(2) = rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.col(2);
        rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(0) = (rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(1).cross(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.col(2))).normalized();

        rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4 = rotationPart_forearmFrameAtElbowAfterJoint4InWristFrame.transpose();
    }
    else{ //nonsingular case, projection for unitary matrix accounting for computation errors
        rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4 = rotationPart_forearmFrameAtElbowAfterJoint4InShoulderFrame.transpose()*wristFrameInShoulderBaseFrame.rotation();

        Eigen::JacobiSVD<Eigen::Matrix3d> svd1(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame, Eigen::ComputeFullU|Eigen::ComputeFullV);
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4, Eigen::ComputeFullU|Eigen::ComputeFullV);

        rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame = svd1.matrixU() * svd1.matrixV().transpose(); // \Sigma = I, J'=U*V^T
        rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4 = svd2.matrixU() * svd2.matrixV().transpose();
    }

    // compute the angles of joint1, joint2, joint3 into column vectors
    Eigen::MatrixXd shoulderJointSet;
    if(fabs(1-fabs(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,1)))<1e-10){ // shoulder singularity occurs
        shoulderJointSet.resize(3,1);

        shoulderJointSet.matrix()(1,0) = acos(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,1));  // 0 or M_PI

        if(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,1)>0){ // theta3 + theta1 = specific value, and assume theta1 is zero.
            shoulderJointSet.matrix()(0,0) = 0;
            shoulderJointSet.matrix()(2,0) = aCosSin(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,0),
                                                     rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,2));
        }
        else{   // theta3 -theta1 = specific value
            shoulderJointSet.matrix()(0,0) = 0;
            shoulderJointSet.matrix()(2,0) = aCosSin(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,0),
                                                     rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,2));
        }
    }
    else{ //shoulder singularity is absent
        shoulderJointSet.resize(3,2);

        shoulderJointSet.matrix()(1,0) = acos(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,1));
        shoulderJointSet.matrix()(1,1)= -acos(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,1));

        shoulderJointSet.matrix()(0,0) = aCosSin(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(2,1)/sin(shoulderJointSet.matrix()(1,0)),
                                                 rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,1)/sin(shoulderJointSet.matrix()(1,0)));
        shoulderJointSet.matrix()(0,1) = aCosSin(rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(2,1)/sin(shoulderJointSet.matrix()(1,1)),
                                                 rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(0,1)/sin(shoulderJointSet.matrix()(1,1)));
        shoulderJointSet.matrix()(2,0) = aCosSin(-rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,2)/sin(shoulderJointSet.matrix()(1,0)),
                                                 rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,0)/sin(shoulderJointSet.matrix()(1,0)));
        shoulderJointSet.matrix()(2,1)= aCosSin(-rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,2)/sin(shoulderJointSet.matrix()(1,1)),
                                                rotationPart_upperArmFrameAtElbowBeforeJoint4InShoulderFrame.matrix()(1,0)/sin(shoulderJointSet.matrix()(1,1)));
    }

    //compute the angles of joint5, joint6, joint7 into column vectors
    Eigen::MatrixXd wristJointSet;
    if(fabs(1-fabs(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0)))<1e-10){ // wrist singularity occurs
        wristJointSet.resize(3,1);

        wristJointSet.matrix()(1,0)= asin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0));
        if(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0)>0){  //theta5 + theta7 = specific value, and assume theta7 is zero
            wristJointSet.matrix()(0,0)= aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,2),
                                                 rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,1));
            wristJointSet.matrix()(2,0) = 0;
        }
        else{ //theta7 - theta5 = specific value, and assume theta7 is zero
            wristJointSet.matrix()(0,0) = - aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,2),
                                                    rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,1));
            wristJointSet.matrix()(2,0) = 0;
        }
    }
    else{  //wrist singularity is absent
        wristJointSet.resize(3,2);

        if(fabs(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0))<1e-10){ // when theta6 equals zero, asin can only offer one solution of 0
            wristJointSet.matrix()(1,0) = 0;
            wristJointSet.matrix()(1,1) = M_PI;
        }
        else{ // expand the solution space of [-M_PI/2, M_PI/2] to [-M_PI, M_PI]
            wristJointSet.matrix()(1,0)= asin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0)); // [-M_PI/2, M_PI/2]
            double flag = (rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0)>0)?(1.0):(-1.0);
            wristJointSet.matrix()(1,1) = flag*M_PI - asin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,0)); //[-M_PI, -M_PI/2) U (M_PI/2, M_PI]
        }

        wristJointSet.matrix()(0,0)  = aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(0,0)/cos(wristJointSet.matrix()(1,0)),
                                               -rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,0)/cos(wristJointSet.matrix()(1,0)));
        wristJointSet.matrix()(0,1)  = aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(0,0)/cos(wristJointSet.matrix()(1,1)),
                                               -rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(2,0)/cos(wristJointSet.matrix()(1,1)));
        wristJointSet.matrix()(2,0)  = aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,1)/cos(wristJointSet.matrix()(1,0)),
                                               -rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,2)/cos(wristJointSet.matrix()(1,0)));
        wristJointSet.matrix()(2,1)  = aCosSin(rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,1)/cos(wristJointSet.matrix()(1,1)),
                                               -rotationPart_wristFrameInForearmFrameAtElbowAfterJoint4.matrix()(1,2)/cos(wristJointSet.matrix()(1,1)));
    }
    shoulderJointSet.row(1) =  - shoulderJointSet.row(1); // the srs kinematic for the human arm have two opposite rotation direction as referred to the euler angle.
    wristJointSet.row(1) = - wristJointSet.row(1);        // 1(y), 2(-x), 3(y), 4(z), 5(y), 6(-z), 7(x)

    unsigned int totalCols = shoulderJointSet.cols()*wristJointSet.cols();  //each column is a SRS IK solutions
    unsigned int totalRows = 7;
    IKSolutionSet.resize(totalRows,totalCols);
    IKSolutionSet.block(0,0,shoulderJointSet.rows(),totalCols) = kron(shoulderJointSet, Eigen::MatrixXd::Ones(1,wristJointSet.cols()));
    IKSolutionSet.block(shoulderJointSet.rows(),0,1,totalCols) = Eigen::MatrixXd::Ones(1,shoulderJointSet.cols()*wristJointSet.cols())*theta4;
    IKSolutionSet.block(shoulderJointSet.rows()+1,0,wristJointSet.rows(),totalCols) = kron(Eigen::MatrixXd::Ones(1,shoulderJointSet.cols()), wristJointSet);

    return true;
}

enum t_iterationStatus KDLIKWithIndependentCoordInitial(urdfModel* highFidelityArmModel,
                                                        const Eigen::VectorXd& initialForIndependentCoordinates,
                                                        Eigen::Isometry3d& _desiredWristPoseInBase,
                                                        Eigen::MatrixXd& W_fromIndependentCoordToDependentCoord,
                                                        Eigen::VectorXd& refinedIKSolution){
    // modified from kdl ik [int ChainIkSolverPos_LMA::CartToJnt(const KDL::JntArray& q_init, const KDL::Frame& T_base_goal, KDL::JntArray& q_out) in chainiksolverpos_lma.cpp]
    // https://github.com/orocos/orocos_kinematics_dynamics
    // benefiting from the initials close to the desired ik solution 
    // jacobian-based iterative IK method can be regarded as a twist-diff based motion plan. As a result, its good iteration step may suffer from the singularity of the jacobian matrix at a singular configuration.
    // Therefore the increment of the joint may be too small even the pose diff is large, where the desired pose can not be reached with a proper iteration step.
    // The only solution is to offer a closer initial joint configuration to reach the desired pose without singularity configuration at the iteration path.

    bool IKInformationPrint=false;

    double v      = 2;
    double tau    = 10;
    double rho;
    double lambda;
    double delta_pos_norm;
    Eigen::Matrix<ScalarType,6,1> delta_pos;
    Eigen::Matrix<ScalarType,6,1> delta_pos_new;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    Eigen::MatrixXd jac;
    Eigen::VectorXd original_Aii;
    Eigen::VectorXd lastSV;
    Eigen::VectorXd q_out;
    Eigen::VectorXd tmp;
    Eigen::VectorXd diffq, diffqOfIndependentCoord;
    Eigen::VectorXd grad;
    Eigen::VectorXd q_new;
    enum t_iterationStatus returnCode;
    returnCode=InitialStatusFlag;

    double eps= 1E-5;
    unsigned int maxiter=500;
    double eps_joints = 1E-15;

    int lastNrOfIter=0;
    double lastDifference;
    double lastTransDiff;
    double lastRotDiff;

    Eigen::VectorXd q;
    Eigen::Isometry3d _endPoseInBase;
    Eigen::Matrix<ScalarType,6,1> L;

    //fix the shoulder elevation angle or the elbow_flexion so that the jacobian can be square(6x6)
    int fixedCoordinateIndInLMIterations;
    int fixedCoordinateIndCandidates[2]={1,3};
    int fixedCoordinateIndCandidateNum=2;
    int CurrentFixedCoordinateCandidateInd=0;

    //Loss = K*TransDiff(m) + RotDiff(rad);   K = 1 degree/1mm = 17.45 rad/m.   1 mm and 1 degree error contribute to loss equally
    L.matrix()(0) = 17.45;  // first three elements indicate loss weight of translational part with unit of m
    L.matrix()(1) = 17.45;
    L.matrix()(2) = 17.45;
    L.matrix()(3) = 1; // second three elements indicate loss weight of rotational part with unit of rad
    L.matrix()(4) = 1;
    L.matrix()(5) = 1;

    //initializing joint configuration
    unsigned int musculoskeletalArmUrdfModelCoordNum_wristToBase = highFidelityArmModel->robotChain.getNrOfJoints(); //21 joints
    q.resize(musculoskeletalArmUrdfModelCoordNum_wristToBase);
    q.setZero();
    q.matrix()(0,0)=thoraxRotation_rzInDegree/180*M_PI;
    q.matrix().block(musculoskeletalArmUrdfModelCoordNum_wristToBase-dependentCoordNum,0,dependentCoordNum,1)= W_fromIndependentCoordToDependentCoord * initialForIndependentCoordinates; //column vector

    //compute end pose with initial joint configuration and the diff of end pose
    _endPoseInBase=highFidelityArmModel->FK_getEndEffectorPoseInBaseFrame(q);
    delta_pos = poseDiff_fromFrame1ToFrame2InBase(_endPoseInBase,_desiredWristPoseInBase); // [dx, dy, dz, drx, dry, drz]'

    //weight the tranlationl and rotational components of the pose diff into a scalar loss
    delta_pos=L.asDiagonal()*delta_pos;
    delta_pos_norm = delta_pos.norm();

    //check the initial loss of pose diff. return if it is small enough
    if (delta_pos_norm<eps) {
        q_out        = q;
        returnCode   = E_NOERROR_Inital;
    }
    else
    {
        //get the weighted jacobian matrix with the initial joint configuration for all joints (fixed joints, dependpend joints and independent joints)
        highFidelityArmModel->getJacobian_endEffectorInBaseFrame(q,jac); // void Jacobian::setColumn, [vx, vy, vz, wx, wy, wz]'
        jac = L.asDiagonal()*jac;

        //transform the raw weighted jacobian matrix into the jacobian matrix for independent joints only (remove coefficient for fixed joint)
        jac = jac.matrix().block(0,musculoskeletalArmUrdfModelCoordNum_wristToBase-dependentCoordNum,jac.rows(),dependentCoordNum) * W_fromIndependentCoordToDependentCoord;

        //remove the 2nd column or the 4th column, the jacobian column corresponding to fix the shoulder elevation angle or the elbow_flexion so that the jacobian can be square(6x6) and invertible.
        fixedCoordinateIndInLMIterations = fixedCoordinateIndCandidates[CurrentFixedCoordinateCandidateInd];
        removeColumn(jac,fixedCoordinateIndInLMIterations);

        //start kdl IK iteration: Levenberg-Marquardt algorithm
        lambda = tau;
        double dnorm = 1;
        diffqOfIndependentCoord.resize(independentCoordNum);
        for (unsigned int i=0;i<maxiter;++i)
        {
            lastNrOfIter   = i;

            svd.compute(jac,Eigen::ComputeFullU|Eigen::ComputeFullV);
            original_Aii = svd.singularValues();

            for (unsigned int j=0;j<original_Aii.rows();++j)
                original_Aii(j) = original_Aii(j)/( original_Aii(j)*original_Aii(j)+lambda);  // LM weight between gradient decent and Newton method

            tmp = svd.matrixU().transpose()*delta_pos;
            tmp = original_Aii.cwiseProduct(tmp);
            diffq = svd.matrixV()*tmp;
            grad = jac.transpose()*delta_pos;
            if (IKInformationPrint) {
                std::cout << "------- iteration " << i << " ----------------\n"
                          << "  q              = " << q.transpose() << "\n"
                          << "  weighted jac   = \n" << jac << "\n"
                          << "  lambda         = " << lambda << "\n"
                          << "  eigenvalues    = " << svd.singularValues().transpose() << "\n"
                          << "  difference     = "   << delta_pos.transpose() << "\n"
                          << "  difference norm= "   << delta_pos_norm << "\n"
                          << "  proj. on grad. = "   << grad << "\n";
                std::cout << std::endl;
            }
            dnorm = diffq.lpNorm<Eigen::Infinity>();
            if (dnorm < eps_joints) {
                returnCode = S_INCREMENT_JOINTS_TOO_SMALL;
                break;
            }

            if (grad.transpose()*grad < eps_joints*eps_joints ) {
                returnCode = S_GRADIENT_JOINTS_TOO_SMALL ;
                break;
            }

            // update the dependent coordinates related to the 6 selected independent joints.
            diffqOfIndependentCoord.matrix().block(0,0,fixedCoordinateIndInLMIterations,1) = diffq.matrix().block(0,0,fixedCoordinateIndInLMIterations,1);
            diffqOfIndependentCoord(fixedCoordinateIndInLMIterations)=0;        //the shoulder elevation angle or the elbow_flexion is fixed in the kdl ik iteration
            diffqOfIndependentCoord.matrix().block(fixedCoordinateIndInLMIterations+1,0,independentCoordNum-1-fixedCoordinateIndInLMIterations,1) = diffq.matrix().block(fixedCoordinateIndInLMIterations,0,diffq.size()-fixedCoordinateIndInLMIterations,1);

            // switch the choice of the fixed coordinate
            CurrentFixedCoordinateCandidateInd++;
            CurrentFixedCoordinateCandidateInd=CurrentFixedCoordinateCandidateInd%fixedCoordinateIndCandidateNum;

            q_new                     = q;
            q_new.matrix().block(musculoskeletalArmUrdfModelCoordNum_wristToBase-dependentCoordNum,0,dependentCoordNum,1) = q_new.matrix().block(musculoskeletalArmUrdfModelCoordNum_wristToBase-dependentCoordNum,0,dependentCoordNum,1) + W_fromIndependentCoordToDependentCoord * diffqOfIndependentCoord;

            _endPoseInBase            = highFidelityArmModel->FK_getEndEffectorPoseInBaseFrame(q_new);
            delta_pos_new             = poseDiff_fromFrame1ToFrame2InBase(_endPoseInBase,_desiredWristPoseInBase);
            delta_pos_new             = L.asDiagonal()*delta_pos_new;
            double delta_pos_new_norm = delta_pos_new.norm();
            rho                       = delta_pos_norm*delta_pos_norm - delta_pos_new_norm*delta_pos_new_norm;
            rho                      /= diffq.transpose()*(lambda*diffq + grad);
            if (rho > 0) {
                q               = q_new;
                delta_pos       = delta_pos_new;
                delta_pos_norm  = delta_pos_new_norm;

                if (delta_pos_norm<eps) {
                    returnCode = E_NOERROR;
                    break;
                }
                highFidelityArmModel->getJacobian_endEffectorInBaseFrame(q_new,jac);;
                jac = L.asDiagonal()*jac;
                //transform the raw weighted jacobian matrix into the jacobian matrix for independent joints only (remove coefficient for fixed joint)
                jac = jac.matrix().block(0,musculoskeletalArmUrdfModelCoordNum_wristToBase-dependentCoordNum,jac.rows(),dependentCoordNum) * W_fromIndependentCoordToDependentCoord;

                //remove the 2nd column or the 4th column, the jacobian column corresponding to fix the shoulder elevation angle or the elbow_flexion so that the jacobian can be square(6x6) and invertible.
                fixedCoordinateIndInLMIterations = fixedCoordinateIndCandidates[CurrentFixedCoordinateCandidateInd];
                removeColumn(jac,fixedCoordinateIndInLMIterations);

                double tmp=2*rho-1;
                lambda = lambda*max(1/3.0, 1-tmp*tmp*tmp);
                v = 2;
            }
            else {
                lambda = lambda*v;
                v      = 2*v;
            }
        }
    }

    if(returnCode == InitialStatusFlag){
        lastNrOfIter = maxiter;
        returnCode   = E_MAX_ITERATIONS_EXCEEDED;
    }
    else
        lastSV         = svd.singularValues();

    lastDifference = delta_pos_norm;
    q_out          = q;  // q and q_out contain all 21 joints in high-fidelity arm urdf model (fixed joints, dependpend joints and independent joints)
    refinedIKSolution = q_out;

    _endPoseInBase = highFidelityArmModel->FK_getEndEffectorPoseInBaseFrame(q);
    delta_pos      = poseDiff_fromFrame1ToFrame2InBase(_endPoseInBase,_desiredWristPoseInBase);
    lastTransDiff  = delta_pos.topRows(3).norm();
    lastRotDiff    = delta_pos.bottomRows(3).norm();

    // the best iteration-based IK return E_NOERROR and E_NOERROR_Inital,
    // and other return codes imply a trade-off ik solution.
    return returnCode;
}

bool IKSolutionFeasibilityCheck(Eigen::VectorXd& refinedIKSolution, Eigen::MatrixXd& W_fromDependentCoordToInependentCoord, Eigen::MatrixXd& W_fromIndependentCoordToDependentCoord){
    Eigen::VectorXd IKSolution_independentCoord;
    Eigen::VectorXd IKSolution_dependentCoord;
    bool outOfLimitFlag = false;

    IKSolution_dependentCoord=refinedIKSolution.matrix().block(1,0,refinedIKSolution.size()-1,1); //Eigen::VectorXd is column vector.
    IKSolution_independentCoord = W_fromDependentCoordToInependentCoord * IKSolution_dependentCoord;

    //std::cout<<"Before retifying joint value into [-pi, pi], the independent joint value is: "<<IKSolution_independentCoord.transpose()<<std::endl;
    for(unsigned int i=0; i< IKSolution_independentCoord.size(); i++){
        // rectify joint values to [-pi, pi]
        IKSolution_independentCoord(i) = IKSolution_independentCoord(i)- 2*M_PI*(int)(IKSolution_independentCoord(i)/(2*M_PI));

        if(IKSolution_independentCoord(i) > M_PI)
            IKSolution_independentCoord(i) = IKSolution_independentCoord(i) - 2*M_PI;
        if(IKSolution_independentCoord(i) < -M_PI)
            IKSolution_independentCoord(i) = IKSolution_independentCoord(i) + 2*M_PI;

        if(IKSolution_independentCoord(i)>(jointUpperLimitsInDegreeForURDF[i]*M_PI/180.0) || IKSolution_independentCoord(i)<(jointLowerLimitsInDegreeForURDF[i]*M_PI/180.0))
            outOfLimitFlag =true;
    }
    //std::cout<<"After  retifying joint value into [-pi, pi], the independent joint value is: "<<IKSolution_independentCoord.transpose()<<std::endl;
    refinedIKSolution.matrix().block(1,0,refinedIKSolution.size()-1,1) = W_fromIndependentCoordToDependentCoord*IKSolution_independentCoord;

    return !outOfLimitFlag;
}


// given the cos and sin function value of a angle to calculate the angle value in [-pi, pi)
double aCosSin(double cosValue, double sinValue){
    double angleValue;
    if(sinValue>0){
        if(cosValue>0)
            angleValue = asin(sinValue);
        else
            angleValue = M_PI - asin(sinValue);
    }
    else{
        if(cosValue>0)
            angleValue = asin(sinValue);
        else
            angleValue = -M_PI - asin(sinValue);
    }
    return  angleValue;
}

//matrix Kronecker tensor product
Eigen::MatrixXd kron( Eigen::MatrixXd m1, Eigen::MatrixXd m2 ){
    int m1R,m1C,m2R,m2C;
    m1R = m1.rows();
    m1C = m1.cols();

    m2R = m2.rows();
    m2C = m2.cols();

    Eigen::MatrixXd m3(m1R*m2R,m1C*m2C);

    for (int i = 0; i < m1R; i++) {
        for (int j = 0; j <  m1C; j++) {
            m3.block(i*m2R, j*m2C, m2R, m2C ) = m1(i,j)*m2;
        }
    }
    return m3;
}

/**
 * determines the rotation axis necessary to rotate the frame b1 to the same orientation as frame b2 and the vector
 * necessary to translate the origin of b1 to the origin of b2, and stores the result in a Twist datastructure.
 */
Eigen::Matrix<ScalarType,6,1> poseDiff_fromFrame1ToFrame2InBase(Eigen::Isometry3d& _pose1InBase,
                                       Eigen::Isometry3d& _pose2InBase,
                                       double deltaTime){
    Eigen::Matrix<ScalarType,6,1> motionTwist;
    motionTwist.matrix().block(0,0,3,1)= (_pose2InBase.translation()-_pose1InBase.translation())/deltaTime;

    Eigen::AngleAxisd rotationFromPose1ToPose2InFrame1;
    rotationFromPose1ToPose2InFrame1.fromRotationMatrix(_pose1InBase.rotation().transpose()*_pose2InBase.rotation());
    motionTwist.matrix().block(3,0,3,1)=_pose1InBase.rotation()*rotationFromPose1ToPose2InFrame1.axis()*rotationFromPose1ToPose2InFrame1.angle()/deltaTime;

    return motionTwist;
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}
