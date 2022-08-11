#ifndef __urdfModel_h__
#define __urdfModel_h__

#include <kdl/kdl.hpp>
#include <kdl/chain.hpp>
#include <kdl/tree.hpp>
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainjnttojacsolver.hpp>

#include <stdio.h>
#include <iostream>
#include <string>
#include <Eigen/Geometry>

class urdfModel{

public:
    ~urdfModel();
    urdfModel();
    urdfModel(const std::string& urdfFile, const std::string& baseLink, const std::string& endLink);

    Eigen::Isometry3d FK_getEndEffectorPoseInBaseFrame(Eigen::VectorXd& _jointsValue);
    void getJacobian_endEffectorInBaseFrame(Eigen::VectorXd& _jointsValue,
                                            Eigen::MatrixXd& _jacobianMatrixOutput);
    Eigen::Isometry3d getPoseOfSpecificChildFrameInParentFrame(const std::string& _parentLinkName,
                                                               const std::string& _childLinkName,
                                                               Eigen::VectorXd& _subChainJointsValue);
    KDL::Tree robotTree;
    KDL::Chain robotChain;
    KDL::ChainFkSolverPos_recursive* FKsolver;
    KDL::ChainJntToJacSolver* jntToJacSovler;
};

#endif // __urdfModel_h__
