#ifndef __hybridIK_h__
#define __hybridIK_h__

#include<kdl/chainiksolver.hpp>
#include<math.h>
#include"urdfModel.h"

// caution: eigen matrix operator overload of "=" is achieved via std::fill_n, and cast for element data type is essential.

enum t_iterationStatus{  // tags for KDL IK with LM and dependent (mimic) joints.
/// Converged but degraded solution (e.g. WDLS with psuedo-inverse singular)
    E_DEGRADED         = +1,
//! No error
    E_NOERROR          =  0,
//! Failed to converge
    E_NO_CONVERGE      = -1,
//! Undefined value (e.g. computed a NAN, or tan(90 degrees) )
    E_UNDEFINED        = -2,
//! Chain size changed
    E_NOT_UP_TO_DATE = -3,
//! Input size does not match internal state
    E_SIZE_MISMATCH = -4,
//! Maximum number of iterations exceeded
    E_MAX_ITERATIONS_EXCEEDED = -5,
//! Requested index out of range
    E_OUT_OF_RANGE = -6,
//! Not yet implemented
    E_NOT_IMPLEMENTED = -7,
//! Internal svd calculation failed
    E_SVD_FAILED = -8,
//! No error with the initials
    E_NOERROR_Inital = -9,
//! joint gradient too small with success
    S_GRADIENT_JOINTS_TOO_SMALL = -10,
//! jointn increment too small with success
    S_INCREMENT_JOINTS_TOO_SMALL = -11,
//! Iteration initial status flag
    InitialStatusFlag = -100
};

typedef double ScalarType;

const unsigned int jointNum = 7;
const unsigned int independentCoordNum=7;
const unsigned int dependentCoordNum=20;
const double thoraxRotation_rzInDegree = 6.8557;

//{elevation plane angle, shoulder elevation, shoulder rotation, elbow flexion, forearm rotation, wrist deviation, wrist flexion}
// The adopted arm musculoskeletal model is from https://simtk.org/projects/upexdyn

// For consistency between SRS IK solution, dependent coordinates in urdf file, independent coordinates in urdf file (W_fromDependentCoordToInependentCoord) are different from the special definition of the independent coordinates in osim file.
// The range of motion of independent coordinte deviation in urdf file is twice the RoM of the independent coordinate deviation in osim, from [-10,25] to [-20,50]
const double jointUpperLimitsInDegreeForURDF[jointNum]={130,180,120,130,90,50,70}; // this RoM is for the independent coordinates in urdf file.
const double jointLowerLimitsInDegreeForURDF[jointNum]={-95,0,-90,0,-90,-20,-70};

const double jointUpperLimitsInDegreeForOSIM[jointNum]={130,180,120,130,90,25,70}; // this RoM is for the independent coordinates in osim file.
const double jointLowerLimitsInDegreeForOSIM[jointNum]={-95,0,-90,0,-90,-10,-70};

/*
* diffrence of deviation and flexion independent coordinates between urdf definition and osim definition.
* 1, urdf denifinition is used in hybrid IK, while osim definition is used for set the coordinate values as referred to the mot file for muscle effort estimation with opensim.
* 2, in urdf definition, the identical rotation angles of the two dependent joints, proximal_row2radius_deviation and  hand2proximal_row_wrist_hand_r1, sum up to the rotation angle of the independent coordinate deviation.
* 3, in urdf definition, the identical rotation angles of the two dependent joints, proximal_row2radius_flexion and hand2proximal_row_wrist_hand_r3, sum up to the rotation angle of the independent coordinate flexion.
*
* 4, in osim definition, the value of independent coordinate deviation is equal to the actual rotation angle of the 17th and 19th dependent joints, proximal_row2radius_deviation and  hand2proximal_row_wrist_hand_r1. (No r_z configuration)
* 5, in osim definition, the value of independent coordinate flexion is twice the actual rotation angle of the 18th and 20th dependent joints, proximal_row2radius_flexion and hand2proximal_row_wrist_hand_r3. (No r_z configuration)
* 6, in mot file for opensim simulation, values of the two independent coordinates deviation and flexion according to the osim definition are needed to be filled in.
*/


class hybridIK{
public:
    ~hybridIK();
    hybridIK();
    hybridIK(const std::string& _urdfFile, const std::string& _simplifiedSRSFile);

    bool getIKSolutionForMusculoskeletalArmModel(double armSwivelAngle, Eigen::Isometry3d& _desiredWristPoseInBase, Eigen::MatrixXd& jointsValue, int ind);

    Eigen::Isometry3d SRSShoulderFrameInBaseFrame;
    Eigen::Vector3d SRSGravityDirectionInBaseFrame;
    double  SRSUpperarmLength, SRSForearmLength;

    Eigen::MatrixXd K_fromSRSInitialToArmIndependentCoord;  //these three linear transformation matrix are only used for urdf definition.
    Eigen::MatrixXd W_fromIndependentCoordToDependentCoord;
    Eigen::MatrixXd W_fromDependentCoordToInependentCoord;

    urdfModel* armModel_urdf;  //arm model with equivalent kinematics of the musculoskeletal arm osim model for jacobian computation using robotic library
    urdfModel* simplifiedArmModel_urdf;  //simplified SRS arm model for compute the initial guess of the proposed hybridIK algorithm
};

bool analyticalIKSolverForSRSArm(Eigen::Isometry3d& shoulderBaseFrame,
                                 Eigen::Isometry3d& wristEndFrame,
                                 Eigen::Vector3d& gravityDirectionInWorldReferenceFrame,
                                 double upperArmLength,
                                 double foreArmLength,
                                 double armSwivelAngle,
                                 Eigen::MatrixXd& IKSolutionSet);

enum t_iterationStatus KDLIKWithIndependentCoordInitial(urdfModel* highFidelityArmModel,   //independent coordinates are in urdf definition.
                                                        const Eigen::VectorXd& initialForIndependentCoordinates,
                                                        Eigen::Isometry3d& _desiredWristPoseInBase,
                                                        Eigen::MatrixXd& W_fromIndependentCoordToDependentCoord,
                                                        Eigen::VectorXd& refinedIKSolution);

bool IKSolutionFeasibilityCheck(Eigen::VectorXd& refinedIKSolution, Eigen::MatrixXd& W_fromDependentCoordToInependentCoord, Eigen::MatrixXd& W_fromIndependentCoordToDependentCoord);  //independent coordinates are in urdf definition.

double aCosSin(double cosValue, double sinValue);
Eigen::MatrixXd kron( Eigen::MatrixXd m1, Eigen::MatrixXd m2);

Eigen::Matrix<ScalarType,6,1> poseDiff_fromFrame1ToFrame2InBase(Eigen::Isometry3d& _pose1InBase,
                                       Eigen::Isometry3d& _pose2InBase,
                                       double deltaTime=1.0);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

inline double max(double a,double b) {
    if (b<a)
        return a;
    else
        return b;
}
#endif // __hybridIK_h__
