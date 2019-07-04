#include <cstdlib>
#include <iostream>
#include "anyoption.h"
#include "BayesianNetwork.h"
#include "factorTree.h"
#include "BayesianNetwork.h"
#include "io_utils.h"
#include "utils.h"
#include "varElimOrder.h"

using namespace qtm;
AnyOption* readOptions(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    
    /* process input params */
    AnyOption *opt = readOptions(argc,argv);
    string filenameBN = opt->getValue("bn_file");
    string elimOrderOpt = opt->getValue("elimination_order");
    int indexSize = strToInt(opt->getValue("k"));
    
    /* read Bayesian network */
    qtm::BayesianNetwork *myNet = new qtm::BayesianNetwork();
    std::cout << "reading input Bayesian Network " << std::endl;
    read_bif_format(myNet, filenameBN);
    
    int nrVars = myNet->nrNodes;
    std::cout << "total nr of variables in BN: " << nrVars << std::endl;
    std::cout << "total nr of edges in BN: " << myNet->nrEdges << std::endl;
    
    /* process query arguments */
    string queryByOpt = opt->getValue("queryBy");
    Query q;
    int queryId = 0;
    if (queryByOpt.compare("variable_label") == 0) {
        process_queryByLabel_arguments(myNet, opt->getValue("query"), q.queryVars, q.evidenceStates);
    }
    else if (queryByOpt.compare("variable_id") == 0) {
        process_queryById_arguments(myNet, opt->getValue("query"), q.queryVars, q.evidenceStates);
    }
    else if (queryByOpt.compare("query_file") == 0)  {
        queryId = process_queryFromFile_arguments(myNet, opt->getValue("query"), q.queryVars, q.evidenceStates);
    }
    else {
        std::cerr << "ERROR: unrecognized queryBy option, exiting " << std::endl;
        exit(1);
    }
    
    std::cout << "*** starting query processing for queryId: " << queryId << ", k: " << indexSize <<  ", and order: " << elimOrderOpt << " ***" << std::endl;
    
    /* compute elimination ordering */
    std::vector<int> elimOrder;
    elimOrder.reserve(myNet->nrNodes);
    computeEliminationOrder(myNet, elimOrder, elimOrderOpt);
    
    /* create elimination tree and perform inference */
    long double tableCreationTimes = 0.0;
    
    FactorTree *factorTree = new FactorTree(myNet, elimOrder, q, elimOrderOpt, indexSize);
    
    //    // control ////////
    //    for (int i = 0; i < myNet->nodesVector.size(); i++) {
    //        std::cout << i << " : " << myNet->nodesVector[i]->nrStates << std::endl;
    //    }
    //    for (int i = 0; i < factorTree->postOrder.size(); i++) {
    //        Factor *f = factorTree->postOrder[i];
    //        if(f->isLeaf)
    //            continue;
    //        std::cout << f->factorId << ": (" ;
    //        for (int j = 0; j < f->scope.size(); j++) {
    //            std::cout << f->scope[j] << ", ";
    //        }
    //        std::cout << ") & " ;
    //        std::cout << f->factorTableLength << " - " << f->cost << std::endl;
    //    }
    //    factorTree->printPrettySecondary(factorTree->rootFactor, "", true);
    
    bool withinTimeLimit = true;
    
    if(indexSize == 0) {
        std::cout << "executing plain VE" << std::endl;
        withinTimeLimit = factorTree->plainVariableElimination(myNet, tableCreationTimes);
    }
    else {
        std::cout << "executing indexed VE with indexSize " << indexSize << std::endl;
        withinTimeLimit = factorTree->indexedVariableElimination(myNet, tableCreationTimes);
    }
    
    std::cout << "total time taken by VE (in sec) " << tableCreationTimes << std::endl;
    
    size_t peakMemoryUsage = getPeakRSS() / 1048576 ; // in MB
    std::cout << "peak memory (in MB) " << peakMemoryUsage << std::endl;
    
    /* writing query results to screen */
    std::cout << " *** printing query results *** " << std::endl;
    std::cout << "scope " << factorTree->rootFactor->scope << std::endl;
    for (int i = 0; i < factorTree->rootFactor->factorTableProbs.size(); i++) {
        std::cout << factorTree->rootFactor->factorTableProbs[i] << std::endl;
    }
    
    /* writing query processing statistics to file */
    string command1 = string("mkdir -p qpResults");
    system(command1.c_str());
    
    string qPFilename = std::string("qpResults") + OS_SEP + std::string("qp_") + elimOrderOpt + std::string("_k") + intToStr(indexSize) + std::string("_q") + intToStr(queryId) + std::string(".txt");
        
    ofstream outMasterStream;
    
    if(outMasterStream.is_open())
        outMasterStream.close();
    
    outMasterStream.open(qPFilename, ios::out | ios::app);
    
    if (outMasterStream.is_open() == false) {
        cerr << "ERROR: cannot open " << qPFilename << " for writing " << endl;
        exit(1);
    }
    
    if (queryId == 0) {
        outMasterStream << "queryId k r nrUseful nrSkipped nrCreated benefitVsCost runtime peakMemory\n";
    }
    
    int querySize = q.queryVars.size();
    
    if(withinTimeLimit){
        double benefitVsCost = 0;
        if (indexSize > 0) {
            benefitVsCost= 100.0 * (factorTree->totalObtainedBenefit / (double) factorTree->totalComputationCost);
        }
        outMasterStream << queryId << " " << indexSize << " " << querySize << " " << factorTree->nrUseful << " " << factorTree->nrSkipped << " " << factorTree->nrCreated << " " << benefitVsCost << " " << tableCreationTimes << " " << peakMemoryUsage << "\n";
    }
    
    else{
        outMasterStream << queryId << " " << indexSize << " " << querySize << " " << "-1" << " " << "-1" << " " << "-1" << " " << "-1" << " " << "-1" << " " << peakMemoryUsage << "\n";
    }
    
    
    outMasterStream.close();
}


AnyOption* readOptions(int argc, char* argv[]) {
    
    AnyOption *opt = new AnyOption();
    
    // ignore POSIX style options
    opt->noPOSIX();
    opt->setVerbose(); /* print warnings about unknown options */
    opt->autoUsagePrint(true); /* print usage for bad options */
    
    opt->addUsage("");
    opt->addUsage("Usage: ");
    opt->addUsage("");
    opt->addUsage("-help Prints this help ");
    opt->addUsage(" -c <config_file> Specify config file ");
    opt->addUsage("");
    
    opt->setOption("bn_file"); // path to bif file
    opt->setOption("elimination_order");
    opt->setOption("queryBy");
    opt->setOption("query");
    opt->setOption("k"); // budget for precomputation
    opt->setCommandFlag("help");
    opt->setCommandOption("c");
    opt->processCommandArgs(argc, argv);
    
    if(opt->getFlag( "help" )) {
        opt->printUsage();
        delete opt;
        exit(0);
    }
    
    const char* configFile = opt->getValue("c");
    if (configFile == NULL) {
        cout << "Config file not mentioned" << endl;
        opt->printUsage();
        delete opt;
        exit(0);
    }
    
    cout << "Config file : " << configFile << endl;
    opt->processFile(configFile);
    opt->processCommandArgs( argc, argv );
    cout << endl << endl;
    return opt;
    
}
