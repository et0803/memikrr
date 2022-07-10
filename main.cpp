#include<armSwivelOptimization.h>
#include <ctime>  // clock(), clock_t, CLOCKS_PER_SEC
#include <math.h>

// inputed arguments
std::string staticOptimizationSetupXMLFile;
std::string urdfFile;
std::string simplifiedSRSFile;
double optimizaitonInitialGuessInterval_InDegree;
double localStopCriterion_kappa;
double globalStopCriterion1_epsilon;
double globalStopCriterion2_delta;
enum globalStopCriterionType globalStopCriterionSelectionTag;
std::string wristPosesFile;
std::string optimizationResultOutputFile;

//parse command line argument
bool parseArguments(int argc, char** argv){
    char str[PATH_MAX];
    if(argc==11){
        realpath(argv[1],str);
        staticOptimizationSetupXMLFile=str;

        realpath(argv[2],str);
        urdfFile=str;

        realpath(argv[3],str);
        simplifiedSRSFile=str;

        optimizaitonInitialGuessInterval_InDegree=atof(argv[4]);
        localStopCriterion_kappa=atof(argv[5]);
        globalStopCriterion1_epsilon=atof(argv[6]);
        globalStopCriterion2_delta=atof(argv[7]);
        globalStopCriterionSelectionTag=atoi(argv[8])<2?stopCriterionOne:stopCriterionTwo;

        realpath(argv[9],str);
        wristPosesFile=str;

        realpath(argv[10],str);
        optimizationResultOutputFile=str;

        std::cout<<"intput arguments success"<<std::endl;
        std::cout<<"1: static optimization setup for arm osim model file relative path:"<< staticOptimizationSetupXMLFile<<std::endl;
        std::cout<<"2: musculoskeletal arm urdf model file relative path:"<<urdfFile<<std::endl;
        std::cout<<"3: simplified SRS arm urdf model file relative path:"<<simplifiedSRSFile<<std::endl;
        std::cout<<"4: optimizaiton initial guess interval of arm swivel angle in degree:"<<optimizaitonInitialGuessInterval_InDegree<<std::endl;
        std::cout<<"5: local optimization stop criterion kappa:"<<localStopCriterion_kappa<<std::endl;
        std::cout<<"6: global optimization stop criterion1 epsilon:"<<globalStopCriterion1_epsilon<<std::endl;
        std::cout<<"7: global optimization stop criterion2 delta:"<<globalStopCriterion2_delta<<std::endl;
        std::cout<<"8: the index of the selected global stop criterion:"<<globalStopCriterionSelectionTag<<std::endl;
        std::cout<<"9: wrist poses file relative path:"<<wristPosesFile<<std::endl;
        std::cout<<"10: optimization result output file relative path:"<<optimizationResultOutputFile<<std::endl;

        return true;
    }
    else{
        std::cout<<"intput arguments error"<<std::endl;
        std::cout<<"1: static optimization setup for arm osim model file relative path (string)"<<std::endl;
        std::cout<<"2: musculoskeletal arm urdf model file relative path (string)"<<std::endl;
        std::cout<<"3: simplified SRS arm urdf model file relative path (string)"<<std::endl;
        std::cout<<"4: optimizaiton initial guess interval of arm swivel angle in degree (double)"<<std::endl;
        std::cout<<"5: local optimization stop criterion kappa (double)"<<std::endl;
        std::cout<<"6: global optimization stop criterion1 epsilon (double)"<<std::endl;
        std::cout<<"7: global optimization stop criterion2 delta (double)"<<std::endl;
        std::cout<<"8: the index of the selected global stop criterion (int)"<<std::endl;
        std::cout<<"9: wrist poses file relative path (string)"<<std::endl;
        std::cout<<"10: optimization result output file relative path (string)"<<std::endl;
        return false;
    }
}

// load wrist poses for arm posture optimization from txt file: row format [id, qx,qy,qz,qw,x,y,z]
bool loadWristPosesFromFile(std::vector<Eigen::Isometry3d>& _wristPosesForArmPostureOptimization)
{
    std::ifstream file;
    file.open(wristPosesFile.c_str());
    std::string str;
    int poseID;

    Eigen::Quaterniond q;
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Isometry3d Trans=Eigen::Isometry3d::Identity();

    while(getline(file, str) && str!=""){     //[id, qx,qy,qz,qw,x,y,z]
        std::istringstream strStream(str);
        double extractedData;

        strStream >> extractedData;
        poseID=extractedData;

        strStream >> extractedData;
        q.x()=extractedData;
        strStream >> extractedData;
        q.y()=extractedData;
        strStream >> extractedData;
        q.z()=extractedData;
        strStream >> extractedData;
        q.w()=extractedData;

        strStream >> extractedData;
        t.x()=extractedData;
        strStream >> extractedData;
        t.y()=extractedData;
        strStream >> extractedData;
        t.z()=extractedData;

        R=q.normalized().toRotationMatrix();
        Trans.setIdentity();
        Trans.rotate(R);
        Trans.pretranslate(t);
        _wristPosesForArmPostureOptimization.push_back(Trans);

        std::cout<<"pose id: "<<poseID<<std::endl<<Trans.matrix() <<std::endl;
    }
    if(_wristPosesForArmPostureOptimization.size()){
        std::cout<<"wrist pose file load success"<<std::endl;
        return true;
    }
    else{
        std::cout<<"wrist pose file load failure"<<std::endl;
        return true;
    }
}

int main(int argc, char** argv)
{
    // parse command line arguments into global variable
    if(!parseArguments(argc, argv))
        return -1;
    //load the pose file to vector wristPosesForArmPostureOptimization
    std::vector<Eigen::Isometry3d> wristPosesForArmPostureOptimization;
    loadWristPosesFromFile(wristPosesForArmPostureOptimization);
    Eigen::VectorXd timeCostForWristPoses(wristPosesForArmPostureOptimization.size());

    //registered type.
    RegisterTypes_osimActuators();
//    OpenSim::Array<std::string> typeNames;
//    OpenSim::Object::getRegisteredTypenames(typeNames);
//    for(int i=0;i<typeNames.size();i++)
//        std::cout<<typeNames[i]<<std::endl;

    // create optimizaiton for muscle-effort-minimization kinematic redundancy resolution of the musculoskeletal arm model
    armSwivelOptimization optimizationObj(staticOptimizationSetupXMLFile,   //in the xml setup file for the model_file path can be relative to the xml setup file, and the coordinate file should be relative to the executable
                                          urdfFile,
                                          simplifiedSRSFile,
                                          optimizaitonInitialGuessInterval_InDegree,
                                          localStopCriterion_kappa,
                                          globalStopCriterion1_epsilon,
                                          globalStopCriterion2_delta,
                                          globalStopCriterionSelectionTag);

    //start optimization for each hand pose, record the time cost and save the optimization path to a result file
    for(unsigned long i=0; i<wristPosesForArmPostureOptimization.size();i++){
        optimizationObj.setDesiredWristPose(i,wristPosesForArmPostureOptimization.at(i));
        std::clock_t startTime = std::clock();
        optimizationObj.runTwoPhaseOptimization();
        std::clock_t endTime = std::clock();
        timeCostForWristPoses(i)=(endTime-startTime)*1.0/CLOCKS_PER_SEC;
        std::cout<<i<<"th wrist pose costs  "<<timeCostForWristPoses(i)<< " seconds!"<<std::endl;
        optimizationObj.appendLocalMinimusToFile(optimizationResultOutputFile);
    }

    // overall time cost statistics
    std::cout<<"total time cost for "<<wristPosesForArmPostureOptimization.size()<< "wrist poses:"<<std::endl<<timeCostForWristPoses.transpose()<<std::endl;
    double meanTimeCost = timeCostForWristPoses.array().mean();
    double varTimeCost = sqrt((timeCostForWristPoses.array()-meanTimeCost).square().sum()/(timeCostForWristPoses.size()-1));
    std::cout<<"performance  of the optimization algorithm (time cost) is: mean("<< meanTimeCost <<") and std variance("<<varTimeCost<<")"<<std::endl;

    return 0;
}
