#include <cstdlib>
#include <iostream>
#include "anyoption.h"
#include "BayesianNetwork.h"
#include "factorTree.h"
#include "BayesianNetwork.h"
#include "io_utils.h"
#include "utils.h"
#include "varElimOrder.h"
#include "dp.h"

using namespace qtm;
AnyOption* readOptions(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    
    /* process input params */
    AnyOption *opt = readOptions(argc,argv);
    string filenameBN = opt->getValue("bn_file");
    string elimOrderOpt = opt->getValue("elimination_order");
    string workloadType = opt->getValue("workload_type");
    int indexSize = strToInt(opt->getValue("k"));
    string indexFolderName = std::string("index_files_") + elimOrderOpt + std::string("_") + intToStr(indexSize);
    
    if (indexSize == 0) {
        std::cerr << "nr of factors to index cannot be 0, exiting " << std::endl;
        exit(1);
    }

    /* read Bayesian network */
    qtm::BayesianNetwork *myNet = new qtm::BayesianNetwork();
    std::cout << "reading input Bayesian Network " << std::endl;
    read_bif_format(myNet, filenameBN);
    myNet->setNrParameters(); // compute nr parameters
    
    std::cout << "total nr of variables in BN: " << myNet->nrNodes << std::endl;
    std::cout << "total nr of edges in BN: " << myNet->nrEdges << std::endl;
    std::cout << "total nr of parameters of BN: " << myNet->getNrParameters() << "\n";
    
//    // control CPTs
//    CPT::iterator it = myNet->nodesVector[0]->condTable.begin();
//    while(it != myNet->nodesVector[0]->condTable.end()) {
//        std::cout << "probs: " << it->second << std::endl;
//        it++;
//    }
//    exit(1);
    
    // check connectivity
    bool isConnected = true;
    for (int i = 0; i < myNet->nodesVector.size(); i++) {
        if (myNet->nodesVector[i]->children.size() == 0 && myNet->nodesVector[i]->parents.size() == 0) {
            std::cout << "this Bayesian node is not connected: " << myNet->nodesVector[i]->nodeLabel << std::endl;
            isConnected = false;
        }
    }
    if (!isConnected) {
        std::cout << "input Bayesian Network should be connected, so exiting" << std::endl;
        exit(1);
    }

    /* compute elimination ordering */
    std::vector<int> elimOrder;
    elimOrder.reserve(myNet->nrNodes);
    computeEliminationOrder(myNet, elimOrder, elimOrderOpt);
    std::cout << "starting with index selection for order " << elimOrderOpt << std::endl;
    /* create elimination tree and perform precomputation */
    std::chrono::steady_clock::time_point begin_all = std::chrono::steady_clock::now();
    FactorTree *factorTree = new FactorTree(myNet, elimOrder, elimOrderOpt);
    
    // check connectivity again for disconnected components
    if (factorTree->nrSecondaryFactors != myNet->nrNodes) {
        std::cout << "input Bayesian Network should be connected, so exiting" << std::endl;
        exit(1);
    }
    
    // assigns usefulness probabilities
    factorTree->assignUsefulnessProbs(myNet, workloadType);
//    factorTree->printPrettySecondary(factorTree->rootFactor, "", true);
    
    std::chrono::steady_clock::time_point begin_dp = std::chrono::steady_clock::now();
    createIndexPotentialsDP(factorTree, indexSize);
    std::chrono::steady_clock::time_point end_dp = std::chrono::steady_clock::now();
    factorTree->fillAndWriteIndexTables(indexFolderName, indexSize);
    std::chrono::steady_clock::time_point end_all = std::chrono::steady_clock::now();
    long double runningTime_dp = getElapsedWallTime(begin_dp, end_dp);
    long double runningTime_all = getElapsedWallTime(begin_all, end_all);
    std::cout << "total running time of DP (in sec) " << runningTime_dp << std::endl;
    std::cout << "total time taken (in sec) " << runningTime_all << std::endl;
    size_t peakMemoryUsage = getPeakRSS() / 1048576 ; // in MB
    std::cout << "peak memory (in MB) " << peakMemoryUsage << std::endl;
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
    
    // input:
    opt->setOption("k"); // budget for precomputation
    opt->setOption("bn_file"); // path to input bif files
    opt->setOption("elimination_order");
    opt->setOption("workload_type");
    
    // keep these options even not used to work with a single config file for both online processing and preprocessing
//    opt->setOption("queryBy");
//    opt->setOption("query");
//    opt->setOption("use_precomputation");
    
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
