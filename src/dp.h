#ifndef _DP_H
#define _DP_H

//#include <map>
//#include <vector>
//#include <list>
//#include <iostream>
#include <algorithm>
#include <iterator>
#include "factor.h"
#include "factorTree.h"
#include "utils.h"

typedef std::numeric_limits< double > dbl;

namespace qtm {
    
    void createIndexPotentialsDP(FactorTree *factorTree, int indexSize);
    
    void memoizeDPTable(FactorTree *factorTree, int indexSize);

    void computeFValVec(Factor *f, int indexSize);

    void computeFplus(Factor *f, std::vector<std::vector<int> > &possiblePartitions, std::vector<double> &FplusChildSolutions);

    void computeFminus(Factor *f, std::vector<std::vector<int> > &possiblePartitions, std::vector<double> &FminusChildSolutions, int ancestorID);

    double probUseful(Factor *f, Factor *ancestor);

    void getAllPossiblePartitions(std::vector<std::vector<int> > &possiblePartitions, std::vector<int> &childrenST, int maxTotal);
    
    void getMaximizerPartition(Factor *f, int ancestorID, std::vector<int> &targetPartition, int targetTotal);
    
    void cart_product(Vvi& out, Vvi& in, int targetTotal, bool forMaximizer);

    int getAncestorId(Factor *f, Factor *ancestor);

    void constructDPSolution(Factor *f, Factor *ancestor, int solutionSize);
    
}

#endif //
