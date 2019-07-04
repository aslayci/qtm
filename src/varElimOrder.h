#ifndef _ELIM_ORDER_H
#define _ELIM_ORDER_H


#include "BayesianNetwork.h"

namespace qtm {
    
    void computeEliminationOrder(BayesianNetwork *myNet, std::vector<int> &elimOrderVec, string elimOrderOpt);
    int minimum_neighbors(std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars);
    int minimum_weight(BayesianNetwork *myNet, std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars, int maxNrStates);
    int min_fill(std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars);
    int weighted_min_fill(BayesianNetwork *myNet, std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars, int maxNrStates);
    void my_min_fill(BayesianNetwork *myNet, std::vector<int> &elimOrderVec);
    bool isConnected(std::vector<std::vector<int> > &neighbors, int node1, int node2);
    
}

#endif //
