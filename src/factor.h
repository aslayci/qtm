#ifndef _FACTOR_H
#define _FACTOR_H

#include "BayesianNetwork.h"


namespace qtm {
    
    class Factor {
        
        friend std::ostream& operator<<(std::ostream& os, const Vi& vi);
        friend std::ostream& operator<<(std::ostream& os, const Vvi& vvi);
        
    public:
        
        Factor(BayesianNetwork *myNet);
        Factor(BayesianNetwork *myNet, Factor *f1, Factor *f2);
        Factor(BayesianNetwork *myNet, double value);
        ~Factor();
        
        // sum product ops
        void assignDomainInfo(); //offset and size
        Factor* product(Factor *f1, Factor *f2);
        int position_consistent_valuation(vector<int> valuation, Factor *f);
        int position_consistent_valuation(vector<int> valuation, Factor *f, int var, int value);
        
        void next_valuation(vector<int> &valuation);
        
        void sum_out(Factor *f);
        void createFactorTableFromSumProduct();
        void createOrganicFactorTable();
        
        std::string factorLabel;
        int var_to_sum_out;
        int factorId;
        
        BayesianNetwork *myNet;
        std::vector<int> scope;
        int width, size;
        std::unordered_map<int, int> var_to_index;
        std::vector<int> offset;
        bool inScope(int);
        
//        int sec_height; // height of the factor in the tree with only secondary factors
//        int sec_level;
        
        bool isRoot;
        bool isLeaf;
        bool isSecondaryLeaf; // true if it contains children that are organic factors from CPDs
        bool inQuery; // whether it is a free / bound var in query or not
        bool isIndex;
        bool isUseful;
        bool canSkip; // whether we can skip table computations for this node (so that we can jump to precomputed info)
        bool isDummy; // for DP binary conversion
        double tableCreationTime; // to measure the correlation of table creation runtime with costs
        
        std::vector<Factor*> childFactors; // pointers to child factors -- points to "all" children that are organic or secondary
        int getElimParentOrderId(std::vector<int> &eliminationOrder, int elimIndex);
        std::vector<Factor*> ancestors; // needed only for DP -- keeps ancestors in binary-converted-tree
        std::vector<Factor*> secChildFactors; // needed only for DP
        std::vector<Factor*> binaryChildFactors; // needed only for DP 

        
        // ********** for inference ********** //
//        FactorTable factorTable;
//        std::vector<std::vector<int> > factorTableStates;
        std::vector<double> factorTableProbs;
        
        uint64 factorTableLength; // length (nr lines) of the potential table at this factor node (var summed out)
        std::unordered_set<int> factorsInSubtree; // contains secondary and also organic factors which has negative factorIds -- to be able to skip organic leaves when they have useful ancestor
        
        
//        bool createFactorTableFromSumProduct(BayesianNetwork *myNet);
        void readFactorTablesFromIndex(string folderName);
        
        
        // ********** for DP ********** //
        uint64 cost;
        uint64 utility; // b(u) from total cost of subtree
        double usefulAloneProb; //singleton usefulness probability
        double probSubtreeInQuery; // probability that at least one node in subtree is in query
        int subtreeSize; // 0 for the leaves that are organic potentials so not counting those ones! for DP
        
        // outer vector is aligned with the order of u's ancestor vector
        // each inside vector per ancestor is of size min(k,|T_u|), keeping the values of F[u,Kappa,v] for all possible Kappa in [1,min(k,|T_u|)]
        std::vector<std::vector<double> > FVal;
        
        // Finclude[ancestor][kappa] keeps whether node u is included or not, i.e., outcome of F+(u, kappa, ancestor) > F-(u, kappa, ancestor)
        std::vector<std::vector<bool> > Finclude;
        
        void writeIndexFactorTables(string folderName);
        
    };
    
    
}

#endif
