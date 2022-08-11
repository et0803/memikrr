#include"urdfModel.h"

urdfModel::~urdfModel(){
    delete FKsolver;
    delete jntToJacSovler;
}

urdfModel::urdfModel(){
    std::cout<<"can not build the robotic model with empty urdf"<<std::endl;
}

urdfModel::urdfModel(const std::string& urdfFile, const std::string& baseLink, const std::string& endLink){
    // kdl robot model notation: https://medium.com/@sarvagya.vaish/forward-kinematics-using-orocos-kdl-da7035f9c8e
    // frame0 : the default reference frame(world frame, base frame, default frame), dont need to be defined explicitly
    // frame1 : defined related to frame0 with (Rotation& R,const Vector& V)) and fixed with the world base
    // joint1 : joint motion axis defined in frame0 and passing through frame0's origin is used to move the attached frame1 to a new pose. the aligned identical frame1' is attached two base link.
    // segment1 : combination of <<joint1 axis in frame0>> and <<frame1>>, used to put this combination into a chain to form a serially articulated robot.

    if(kdl_parser::treeFromFile(urdfFile,robotTree)){
        std::cout<<"building robot success with urdf file: "<< urdfFile <<std::endl;
        robotTree.getChain(baseLink, endLink, robotChain);
        FKsolver = new KDL::ChainFkSolverPos_recursive(robotChain);
        jntToJacSovler = new KDL::ChainJntToJacSolver(robotChain);
    }
    else
        std::cout<<"building robot failure with urdf file: "<< urdfFile <<std::endl;

}

Eigen::Isometry3d urdfModel::FK_getEndEffectorPoseInBaseFrame(Eigen::VectorXd& _jointsValue){
    Eigen::Isometry3d _endLinkPoseOutput;
    _endLinkPoseOutput.setIdentity();

    if (robotChain.getNrOfJoints()==_jointsValue.size()){
        KDL::JntArray q(_jointsValue.size());
        for(unsigned int i=0;i<_jointsValue.size();i++)
            q(i)=_jointsValue(i);

        KDL::Frame eeFrame;
        FKsolver->JntToCart(q, eeFrame);

        //transform KDL::frame to Eigen::Isometry3d
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                _endLinkPoseOutput.matrix()(i, j) = eeFrame(i,j);
            }
        }
    }
    else{
        std::cout<<"input number of joint value not equal to the robot chain"<<std::endl;
    }
    return _endLinkPoseOutput;
}

void urdfModel::getJacobian_endEffectorInBaseFrame(Eigen::VectorXd& _jointsValue,
                                                   Eigen::MatrixXd& _jacobianMatrixOutput){
    if (robotChain.getNrOfJoints()==_jointsValue.size()){
        KDL::JntArray q(_jointsValue.size());
        for(unsigned int i=0;i<_jointsValue.size();i++)
            q(i)=_jointsValue(i);

        //get the jacbian at _jointsValue configuration
        KDL::Jacobian jac;  //eigen matrix in jac need to be resized before calculate jacobian matrix.

        if(0){   //segement and joint names check.
            std::cout<<"Number of segements is "<<robotChain.getNrOfSegments()<<". and number of joints is "<<robotChain.getNrOfJoints()<<std::endl;
            for(unsigned int i=0;i<robotChain.getNrOfSegments();i++){
                std::cout<<i<<"th segement: "<< robotChain.getSegment(i).getName()<<std::endl;
                std::cout<<i<<"th joint: "<< robotChain.getSegment(i).getJoint().getName()<<std::endl;
            }
        }

        jac.resize(robotChain.getNrOfSegments());
        jntToJacSovler->JntToJac(q,jac);
        _jacobianMatrixOutput=jac.data;
    }
    else{
        std::cout<<"input number of joint value not equal to the robot chain"<<std::endl;
    }
}

Eigen::Isometry3d urdfModel::getPoseOfSpecificChildFrameInParentFrame(const std::string& _parentLinkName,
                                                                      const std::string& _childLinkName,
                                                                      Eigen::VectorXd& _subChainJointsValue){
    KDL::Chain specificLinkChain;
    robotTree.getChain(_parentLinkName, _childLinkName, specificLinkChain);
    KDL::ChainFkSolverPos_recursive specificFK(robotChain);

    KDL::JntArray q(_subChainJointsValue.size());
    for(unsigned int i=0;i<_subChainJointsValue.size();i++)
        q(i)=_subChainJointsValue(i);

    KDL::Frame eeFrame;
    specificFK.JntToCart(q, eeFrame);

    //transform KDL::frame to Eigen::Isometry3d
    Eigen::Isometry3d _endLinkPoseOutput;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            _endLinkPoseOutput.matrix()(i, j) = eeFrame(i,j);
        }
    }
    return _endLinkPoseOutput;
}
