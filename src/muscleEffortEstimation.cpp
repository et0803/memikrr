#include"muscleEffortEstimation.h"

muscleEffortEstimation::~muscleEffortEstimation(){
    delete statesStore;
    delete coordinatesStore;
    free(activiationLevels);
}

muscleEffortEstimation::muscleEffortEstimation(const std::string &_staticOptimizationSetupXMLFile):
    OpenSim::AnalyzeTool(_staticOptimizationSetupXMLFile)// aFileName is the xml setting file for AnalysisSet
{
    setPrintResultFiles(false);   //dont print static optimization result to file
    setLoadModelAndInput(false);   //load model and motion data from files
    statesStore = new OpenSim::Storage(512,"states");
    coordinatesStore = new OpenSim::Storage(getCoordinatesFileName());
    activiationLevels=NULL;
    _model->initSystem();
}

double muscleEffortEstimation::getMuscleEffortUsingOpensimAnalyzaTool(Eigen::VectorXd& dependentJCoordValues){

    double totalSquaredActiviation=0;

    OpenSim::Model* _modelWorkingCopy;
    _modelWorkingCopy = _model->clone();

    // adapt from bool AnalyzeTool::run(bool plotting) to perform static optimization
    SimTK::State& s = _modelWorkingCopy->initSystem();
    createStatesStorageFromCoordinatesInEigenVector(_modelWorkingCopy ,s, dependentJCoordValues);
    _modelWorkingCopy->getMultibodySystem().realize(s, SimTK::Stage::Position );

    bool completed = true;
    try {

    // VERIFY THE CONTROL SET, STATES, AND PSEUDO STATES ARE TENABLE
    verifyControlsStates();

    // SET OUTPUT PRECISION
    OpenSim::IO::SetPrecision(_outputPrecision);

    // ANALYSIS SET
    OpenSim::AnalysisSet& analysisSet = _modelWorkingCopy->updAnalysisSet();
    if(analysisSet.getSize()<=0) {
        std::string msg = "AnalysisTool.run: ERROR- no analyses have been set.";
        throw OpenSim::Exception(msg,__FILE__,__LINE__);
    }

    // Call helper function to process analysis
    /*Array<double> bounds;
    bounds.append(_ti);
    bounds.append(_tf);
    const_cast<Storage &>(aStatesStore).interpolateAt(bounds);*/
    double ti,tf;
    int iInitial = statesStore->findIndex(_ti);
    int iFinal = statesStore->findIndex(_tf);
    statesStore->getTime(iInitial,ti);
    statesStore->getTime(iFinal,tf);

    // It is ridiculous to start before the specified time! So check we aren't doing something stupid.
    //while(ti < _ti){
    //  _statesStore->getTime(++iInitial,ti);
    //}

    OpenSim::log_info("Executing the analyses from {} to {}...", ti, tf);
    run(s, *_modelWorkingCopy, iInitial, iFinal, *statesStore, _solveForEquilibriumForAuxiliaryStates);
    _modelWorkingCopy->getMultibodySystem().realize(s, SimTK::Stage::Position );
    } catch (const OpenSim::Exception& x) {
        x.print(std::cout);
        completed = false;
        throw OpenSim::Exception(x.what(),__FILE__,__LINE__);
    }

    // PRINT RESULTS
    // TODO: give option to write partial results if not completed
    if (completed)
    {
        int index=0;
        OpenSim::Analysis& analysis = _modelWorkingCopy->updAnalysisSet().get(index);  //in printResults of opensim-core/OpenSim/Simulation/Model/AbstractTool.cpp  and opensim-core/OpenSim/Simulation/Model/AnalysisSet.cpp
        if(getName()=="StaticOptimization"){
            OpenSim::StaticOptimization* analysisSO;
            analysisSO = (OpenSim::StaticOptimization*)(&analysis);
            OpenSim::Storage* activiations = analysisSO->getActivationStorage();

            OpenSim::Array<std::string> muslceLabels = activiations->getColumnLabels();
            double ti=activiations->getFirstTime();       //in print of opensim-core/OpenSim/Common/Storage.cpp
            if(activiationLevels==NULL)
            {
                activationNum = activiations->getStateVector(1)->getSize();
                activiationLevels= new double[activationNum];
            }
            activationNum = activiations->getDataAtTime(ti,activationNum, &activiationLevels);   //in print of opensim-core/OpenSim/Common/Storage.cpp

            Eigen::VectorXd activationEigenVector(activationNum);
            for(int i=0;i<activationNum;i++){
                //std::cout<<"the "<<i<<"th muscle: "<<muslceLabels[i+1]<<", the activiation level is "<<activiationLevels[i]<<std::endl;
                activationEigenVector[i] = activiationLevels[i];
                totalSquaredActiviation += activiationLevels[i]*activiationLevels[i];
            }
            std::cout<<std::endl;
            std::cout<<std::endl;
            std::cout<<"Arm coordinates are: "<<dependentJCoordValues.transpose()<<std::endl;
            std::cout<<"Arm muscle acitivations are: "<<activationEigenVector.transpose()<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl;
        }
    }

    delete _modelWorkingCopy;
    return totalSquaredActiviation;
}

void muscleEffortEstimation::createStatesStorageFromCoordinatesInEigenVector(OpenSim::Model* workingModel, SimTK::State& s, Eigen::VectorXd& dependentJCoordValues)
{
    // only configure the _statesStore, a private member variable in class AnalyzeTool;
    statesStore->reset(0);
    OpenSim::Storage *qStore=NULL;
    OpenSim::Storage *uStore=NULL;
    workingModel->getSimbodyEngine().formCompleteStorages(s, *coordinatesStore, qStore, uStore);  //in void AnalyzeTool::loadStatesFromFile(SimTK::State& s), return qStore and uStore in radian


    // assign new values to the qStore and uStore
    int nq = workingModel->getNumCoordinates();
    int nu = workingModel->getNumSpeeds();
    int ny = workingModel->getNumStateVariables();

    OpenSim::Array<std::string> stateNames("", ny);
    OpenSim::Array<std::string> qLabels = qStore->getColumnLabels();
    OpenSim::Array<std::string> uLabels = uStore->getColumnLabels();
    stateNames = workingModel->getStateVariableNames();
    stateNames.insert(0, "time");
//    for(int i=0;i<=ny;i++)
//        std::cout<<"Before checking, the "<<i<<"th column name is "<<stateNames[i]<<std::endl;
    // Preserve the labels from the data file which are typically abbreviated
    // label[0] = time
    for(int i=1; i<=nq; ++i){
        stateNames[i] = qLabels[i];
    }
    for(int i=1; i<=nu; ++i){
        stateNames[i+nq] = uLabels[i];
    }
//    for(int i=0;i<=ny;i++)
//        std::cout<<"After checking, the "<<i<<"th column name is "<<stateNames[i]<<std::endl;

    statesStore->setColumnLabels(stateNames);
    OpenSim::Array<double> y(0.0,ny);
    for(int i=0;i<nq;i++){
        y[i] = dependentJCoordValues[i]; // order of coordinates in dependentJCoordValues should be same to those in coordinates_file from getCoordinatesFileName().
        //std::cout<<"the "<<i+1<<"th column name is "<<stateNames[i+1]<<", and the value is "<<y[i]<<std::endl;

    }

    //set multiple same row for a static motion
    for(int index=0;index<qStore->getSize();index++){
        double t;
        qStore->getTime(index,t);
        statesStore->append(t,ny,&y[0]);
    }
    setStatesStorage(*statesStore);
    delete qStore;
    delete uStore;
}
