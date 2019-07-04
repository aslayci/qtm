#ifndef _FACTOR_TREE_H
#define _FACTOR_TREE_H
#include "BayesianNetwork.h"
#include "factor.h"
#include <iomanip>

namespace qtm {
    
//        typedef std::vector<Factor*> FactorsList;
    
    class FactorTree {
        
    public:
        
        // bu ilk constructora daha sonra evidence da gelecek, simdilik boyle birak
        // constructor for query processing
        FactorTree(BayesianNetwork *myNet, std::vector<int> &eliminationOrder, Query &query, string elimOrderType, int k); // constructor for preprocessing
        FactorTree(BayesianNetwork *myNet,std::vector<int> &elimOrder, string elimOrderType); // constructor for preprocessing
        ~FactorTree();
        
        // factor containers
        Factor *rootFactor;
        std::vector<std::vector<Factor*> > elimOrderMap; // indices inline with eliminationOrder vector
        std::vector<Factor*> postOrder;
        std::vector<Factor*> binaryPostOrder; // required only for DP
        std::vector<Factor*> levelOrder; 
        
        std::vector<int> queryVariables;
        std::vector<int> evidenceStates;
        
        bool deleteThings;
        
        // workload node in query probs assignment
        void assignUsefulnessProbs(BayesianNetwork *myNet, string workloadType);
        void assignUniformProbNodeInQuery(BayesianNetwork *myNet);
        void assignRootBiasedProbNodeInQuery(BayesianNetwork *myNet);
        void writeProbNodeInQueryToFile(BayesianNetwork *myNet, string workloadType);

        
        // keep a flag based on which constructor is called to differentiate DP from inference -- not to perform unnecessary operations
        bool forDP;
        string elimOrderOption;
        
        int nrFactors;
        int nrSecondaryFactors;
        
        void fillInitialFactorTables(); // bu evidence ile de calismali!
        void createInitialFactorsFromCPTs(BayesianNetwork *myNet, std::vector<int> &eliminationOrder);
        void buildFactorEliminationTree(BayesianNetwork *myNet, std::vector<int> &elimOrder);
       
        // print to screen operations
        void printPretty(Factor *f, std::string indent, bool last);
        void printPrettySecondary(Factor *f, std::string indent, bool last);
        void printPrettySecWithLenthRatios(Factor *f, std::string indent, bool last);
        void printGivenLevel(Factor *f, int level);
        void printLevelOrder(Factor *f);
        void printPostorder(Factor *f);
        
        // recursive tree operations
        void getGivenLevel(Factor *f, int level, std::vector<Factor*> &levelOrder);
        void getLevelOrder(Factor *f, std::vector< std::vector<Factor*>> &levels);
        void getPostorder(Factor *f, std::vector<Factor*> &postOrder);
        bool getAncestors(Factor *root, Factor *f, std::vector<Factor*> &ancestors);
        int height(Factor *f);
        int secondaryHeight(Factor *f);
        
        void getBinaryLevelOrder(Factor *f, std::vector< std::vector<Factor*>> &levels);
        void getBinaryGivenLevel(Factor *f, int level, std::vector<Factor*> &levelOrder);
        
        // precomputation savings
        int nrUseful;
        int nrSkipped;
        int nrCreated;
        int64 totalObtainedBenefit; // total benefit coming from useful factors for the query
        int64 totalComputationCost; // total computation cost of computing tables
        
        // inference operations
        bool plainVariableElimination(BayesianNetwork *myNet, long double &tableCreationTimes);
        bool indexedVariableElimination(BayesianNetwork *myNet, long double &tableCreationTimes);
        void findUsefulIndexFactors(Factor *f);
        void readIndexFactorIds(string folderName, std::vector<int> &indexFactorIds);
        std::vector<int> indexFactorIds;
        string folderName;
        int indexSize;
        
        
        // DP operations
        void fillAndWriteIndexTables(string folderName, int indexSize);
        std::unordered_map<int, int> factor_to_poPosition; // secondary factor nodes' ids postorder
        
        void outputCostsForCorrelationAnalysis(BayesianNetwork *myNet);

        // binary conversion operations -- for DP
        void splitToBinary(BayesianNetwork *myNet, Factor *root, vector<Factor*> remainingSecChildren);
        void getBinaryPostorder(Factor *f, std::vector<Factor*> &postOrder);
        bool getBinaryAncestors(Factor *root, Factor *f, std::vector<Factor*> &ancestors);
        void printPrettyBinary(Factor *f, std::string indent, bool last);
        int binaryHeight(Factor *f);
        
        
//        void fillFactorTables(BayesianNetwork *myNet); // this was just for testing
        
        
//        void printPostorder(Factor *f);
//        int treeWidth(Factor *f);
//        void recFindDepths(Factor* f, std::unordered_map<int, int>& levelFreqs, int depth);
//        int height(Factor *f);
//        void printLevelOrder(Factor *f);
//        void printGivenLevel(Factor *f, int level);
        
    };
    
}

#endif
