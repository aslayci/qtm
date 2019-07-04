#include "dp.h"


namespace qtm {
    
    // linkte tum nodelarla run edince sorun cikiyor, o yuzden simdilik asagidaki gibi yapiyorum
    
    void createIndexPotentialsDP(FactorTree *factorTree, int indexSize) {
        std::cout << "starting with DP memoization " << std::endl;
        
        memoizeDPTable(factorTree, indexSize);
        std::cout << "starting with DP solution construction " << std::endl;
        constructDPSolution(factorTree->rootFactor, factorTree->rootFactor->ancestors[0], indexSize);
        
    }
    
    void memoizeDPTable(FactorTree *factorTree, int indexSize) {
        
        for (int i = 0; i < factorTree->binaryPostOrder.size(); i++) {
            
            Factor *f = factorTree->binaryPostOrder[i];
            
            if (f->isLeaf) 
                continue;
           
            computeFValVec(f, indexSize);

//            /* kontrol start */
//            std::cout << "******************************************************************" << std::endl;
//            std::cout << "***** kontrolling FVal table for " << f->factorLabel << " with id " << f->factorId  << std::endl;
//            for (int j = 0; j < f->ancestors.size(); j++) {
//                std::cout << "** printing entries wrt ancestor " << f->ancestors.at(j)->factorLabel << std::endl;
//                for (int z = 0; z < f->FVal[j].size(); z++) {
//                    int snsz = z + 1;
//                    cout.precision(dbl::max_digits10);
//                    std::cout  << snsz << ", " << f->FVal[j][z] << ", " << f->Finclude[j][z] << std::endl;
//                }
//                std::cout << std::endl;
//            }
//            std::cout << std::endl;
//            std::cout << "******************************************************************" << std::endl;
//           /* kontrol end */
            
        }
    
//         exit(1);
    }
    
    void computeFValVec(Factor *f, int indexSize) {
        
        if (f->isSecondaryLeaf) { // direk F+ cagirip oyle halletmek lazim for each ancestor combinationlar ile ugrasmadan
            // compute F+[u,1,v] = B_u(u,v) for each ancestor v
            for (int v = 0; v < f->ancestors.size(); v++) {
                Factor* ancestor = f->ancestors[v];
                f->FVal.push_back(std::vector<double>());
                f->Finclude.push_back(std::vector<bool>());
                double partialSingletonBenefit = probUseful(f, ancestor) * f->utility; // B_u(u,v)
                f->FVal[v].push_back(partialSingletonBenefit);
                f->Finclude[v].push_back(true);
            }
        } // end of if isSecondaryLeaf
        
        else { // if this node contains at least one secondary child
            
            int max_kappa_at_node, maxSolutionSizeFromST;
            if (f->isDummy) {
                max_kappa_at_node = std::min(f->subtreeSize, indexSize);
                maxSolutionSizeFromST = max_kappa_at_node;
            }
            else {
                maxSolutionSizeFromST = std::min(f->subtreeSize - 1, indexSize); // max nr of nodes selectable from this factor's children's subtrees -- excluding organic factors
                max_kappa_at_node = std::min(f->subtreeSize, indexSize);
            }
            
            // get all possible feasible nrSecChildren-partitioning of the number maxSolutionSizeFromST
            std::vector<int> childrenSTSizes;
            for (int i = 0; i < f->binaryChildFactors.size(); i++) { // to work only with secondary factors for DP
                childrenSTSizes.push_back(f->binaryChildFactors[i]->subtreeSize);
            }
            // get all valid partitioning of maxSolutionSizeFromST into children's solution sizes
            std::vector<std::vector<int> > possiblePartitions;
            getAllPossiblePartitions(possiblePartitions, childrenSTSizes, maxSolutionSizeFromST);
            
             // kontrol start
//            std::cout << "******************************************************************" << std::endl;
//            std::cout << " **** checking all possible partitions for factor **** " << f->factorLabel << std::endl;
//            for (int p = 0; p < possiblePartitions.size(); p++) {
//                for (int z = 0; z < possiblePartitions[p].size(); z++) {
//                    std::cout << possiblePartitions[p][z] << " ";
//                }
//                std::cout << std::endl;
//            }
             // kontrol end
            
            
            // compute F+(T_u w/o u, Kappa, u) for each Kappa -- this is agnostic to u's ancestors so computing once
            std::vector<double> FplusChildSolutions(maxSolutionSizeFromST, 0.0);
            if (!f->isDummy) {
                computeFplus(f, possiblePartitions, FplusChildSolutions);
            }
            
            // compute  F- for each ancestor and possible solution size
            for (int v = 0; v < f->ancestors.size(); v++) {
                
                Factor* ancestor = f->ancestors[v];
                
                f->FVal.push_back(std::vector<double>());
                f->Finclude.push_back(std::vector<bool>());
                
                std::vector<double> FminusChildSolutions(maxSolutionSizeFromST, 0.0);
                computeFminus(f, possiblePartitions, FminusChildSolutions, v); // v gets plus 1 for reaching FVal of children of course
                
                
                double partialSingletonBenefit = probUseful(f, ancestor) * f->utility; // B_u(u,v)
                
                if (f->isDummy) {
                    partialSingletonBenefit = 0;
                }
                
//                if (f->factorId == 16 || f->factorId == 622 || f->factorId == 629 ) {
//                    cout.precision(dbl::max_digits10);
//                    std::cout << f->factorId << " " << ancestor->factorId << " " << partialSingletonBenefit << std::endl;
//
//                }
                
                // kontrol start
//                std::cout << "B(" << f->factorLabel << " | " << ancestor->factorLabel << ") = " << partialSingletonBenefit << std::endl;
//                std::cout << probUseful(f, ancestor) << " lan " << f->utility << std::endl;
//                std::cout << "prob alones " << f->usefulAloneProb << " lan " << ancestor->usefulAloneProb << std::endl;
                // kontrol end
                
                // now compare F+ and F- for kappa in [1, maxSolutionSizeFromST]
                int kappaIndex;
                for (kappaIndex = 0; kappaIndex < maxSolutionSizeFromST; kappaIndex++) { // kappaIndex corresponds to solution size kappaIndex + 1
                    // for dummy
                    if (kappaIndex == 0) {
                        if (f->isDummy) {
                            f->FVal[v].push_back(FminusChildSolutions[0]);
                            f->Finclude[v].push_back(false);
                        }
                        else {
                            if (partialSingletonBenefit > FminusChildSolutions[0]) {// if F+(u,1,v) >= F-(u,1,v)
                                f->FVal[v].push_back(partialSingletonBenefit);
                                f->Finclude[v].push_back(true);
                            }
                            else {
                                f->FVal[v].push_back(FminusChildSolutions[0]);
                                f->Finclude[v].push_back(false);
                            }
                        }
                        continue;
                    } // if kappaIndex = 0
                    if (f->isDummy) {
                        f->FVal[v].push_back(FminusChildSolutions[kappaIndex]);
                        f->Finclude[v].push_back(false);
                    }
                    
                    else {
                        if (FplusChildSolutions[kappaIndex - 1] + partialSingletonBenefit > FminusChildSolutions[kappaIndex]) {
                            f->FVal[v].push_back(FplusChildSolutions[kappaIndex - 1] + partialSingletonBenefit);
                            f->Finclude[v].push_back(true);
                        }
                        else {
                            f->FVal[v].push_back(FminusChildSolutions[kappaIndex]);
                            f->Finclude[v].push_back(false);
                        }
                    }
                    
                } // end of for loop
                // cide: binary conversion'dan etkilenmii galiba bu for fummy ones
                if (max_kappa_at_node == maxSolutionSizeFromST + 1) { // indexSize > |T_u| - 1 then compute kappa = |T_u| only from F+
                    f->FVal[v].push_back(FplusChildSolutions[maxSolutionSizeFromST - 1] + partialSingletonBenefit);
                    f->Finclude[v].push_back(true);
                }
            }
            
        } // end of else if this node is not a secondary leaf
    }
    
    // for a given factor compute F+(child_u,K_i,u) for each child for all possible children solution sizes w/o B_u(u,v) part, below code is ancestor agnostic called once hence no redundant calls for each ancestor as it is the same
    void computeFplus(Factor *f, std::vector<std::vector<int> > &possiblePartitions, std::vector<double> &FplusChildSolutions) {
        
        double totValFChildren;
        int total;
        
        for (int i = 0; i < possiblePartitions.size(); i++) {
            total = 0;
            totValFChildren = 0;
            for (int j = 0; j < possiblePartitions[i].size(); j++) {
                total += possiblePartitions[i][j];
                if(possiblePartitions[i][j] != 0) { //if zero then we add nothing since FVal would be 0
                    totValFChildren += f->binaryChildFactors[j]->FVal[0][possiblePartitions[i][j]-1];
                }
            }
            
            
//            if (f->factorId == 7) {
//                std::cout << "f8 pp " << possiblePartitions[i] << std::endl;
//                std::cout << "t0tal " << total<< std::endl;
//                std::cout << "totovalf " << totValFChildren << std::endl;
//                std::cout << "plusChildSolutions[total-1] " << FplusChildSolutions[total-1] << std::endl;
//            }
            
            if (FplusChildSolutions[total-1] < totValFChildren) {
                FplusChildSolutions[total-1] = totValFChildren;
            }
            
        }
        
    }
    
    void computeFminus(Factor *f, std::vector<std::vector<int> > &possiblePartitions, std::vector<double> &FminusChildSolutions, int ancestorID) {
        
        double totValFChildren;
        int total;
        
        for (int i = 0; i < possiblePartitions.size(); i++) {
            total = 0;
            totValFChildren = 0;
            for (int j = 0; j < possiblePartitions[i].size(); j++) {
                total += possiblePartitions[i][j];
                if(possiblePartitions[i][j] != 0) { //if zero then we add nothing since FVal would be 0
                    totValFChildren += f->binaryChildFactors[j]->FVal[ancestorID + 1][possiblePartitions[i][j]-1];
                }
            }
            if (FminusChildSolutions[total-1] < totValFChildren) {
                FminusChildSolutions[total-1] = totValFChildren;
            }
        }
    }
    
    double probUseful(Factor *f, Factor *ancestor) {
        
        if (ancestor->isDummy) { // if dummy node
            return 0;
        }
        // for greedy, we need to check whether this is really the ancestor since greedy performs arbitrary usefulness computations
        return (f->usefulAloneProb - ancestor->usefulAloneProb);
    }
    
    void getAllPossiblePartitions(std::vector<std::vector<int> > &possiblePartitions, std::vector<int> &childrenST, int targetTotal) {
        
        std::vector<std::vector<int> > dene;
        
        for (int i = 0; i < childrenST.size(); i++) {
            dene.push_back(Vi());
            for (int j = 0; j <= childrenST[i]; j++) {
                dene[i].push_back(j);
            }
        }
        cart_product(possiblePartitions, dene, targetTotal, false);
    }
    
    void cart_product(Vvi& out, Vvi& in, int targetTotal, bool forMaximizer) {
        Vd vd;
        for(Vvi::const_iterator it = in.begin();
            it != in.end();
            ++it) {
            Digits d = {(*it).begin(), (*it).end(), (*it).begin()};
            vd.push_back(d);
        }
        
        while(1) {
            Vi result;
            int total = 0;
            for(Vd::const_iterator it = vd.begin(); it != vd.end(); it++) {
                total += (*(it->me));
                result.push_back(*(it->me));
            }
            if (!forMaximizer) {
                if (total != 0 && total <= targetTotal) {
                    out.push_back(result);
                }
            }
            else {
                if (total == targetTotal) {
                    out.push_back(result);
                }
            }
            
            for(Vd::iterator it = vd.begin(); ; ) {
                ++(it->me);
                if(it->me == it->end) {
                    if(it+1 == vd.end()) {
                        return;
                    } else {
                        it->me = it->begin;
                        ++it;
                    }
                } else {
                    break;
                }
            }
        }
    }
    
    void getMaximizerPartition(Factor *f, int ancestorID, std::vector<int> &targetPartition, int targetTotal) {
        
        std::vector<std::vector<int> > possiblePartitions;
        Vvi dene;
        std::vector<int> childrenST;
        for (int i = 0; i < f->binaryChildFactors.size(); i++) { // to work only with secondary factors for DP
            childrenST.push_back(f->binaryChildFactors[i]->subtreeSize);
        }
        
        for (int i = 0; i < childrenST.size(); i++) {
            dene.push_back(Vi());
            for (int j = 0; j <= childrenST[i]; j++) {
                dene[i].push_back(j);
            }
        }
        
        cart_product(possiblePartitions, dene, targetTotal, true);
        
        
//        if (f->factorId == 20) {
//            std::cout << "printing pps for 20 " << std::endl;
//            std::cout << possiblePartitions << std::endl;
//        }

        double totValFChildren = 0.0;
        double maxTotValChildren = 0.0;
        int bestID = -1;

        for (int i = 0; i < possiblePartitions.size(); i++) {
            totValFChildren = 0.0;
            for (int j = 0; j < possiblePartitions[i].size(); j++) {
                if (possiblePartitions[i][j] != 0) {
                    totValFChildren += f->binaryChildFactors[j]->FVal[ancestorID][possiblePartitions[i][j]-1];
//                    if (f->factorId == 20){
//                        std::cout << "totvalf for f 20: " << f->binaryChildFactors[j]->FVal[ancestorID][possiblePartitions[i][j]-1] << std::endl;
//                    }
                }
            }
            if (totValFChildren > maxTotValChildren) {
                maxTotValChildren = totValFChildren;
                bestID = i;
            }
        }
        
        targetPartition = possiblePartitions[bestID];
        
//        if (f->factorId == 20) {
//            std::cout << "bestid for 20 " << bestID << std::endl;
//        }
//
//        if (f->factorId == 20) {
//            std::cout << "printing target for 20 " << std::endl;
//            std::cout << targetPartition << std::endl;
//        }
        
    }
    
    int getAncestorId(Factor *f, Factor *ancestor) {
        
        int ancestorID = -1;
        
        std::vector<Factor*>::iterator it = std::find(f->ancestors.begin(), f->ancestors.end(), ancestor);
        if (it != f->ancestors.end()) {
            ancestorID = std::distance(f->ancestors.begin(), it);
        }
        
        return ancestorID;
    }
    
    
    void  constructDPSolution(Factor *f, Factor *ancestor, int solutionSize) {
       
//        std::cout << "******** start **********" << std::endl;
//        std::cout << f->factorLabel << std::endl;
//        std::cout << "f->binaryChildFactors.size() : " << f->binaryChildFactors.size() << std::endl;
//        std::cout << "solutionSize " <<  solutionSize << std::endl;
        
        
        if (f->isLeaf) // don't consider organic factors
            return;
        
        if(solutionSize == 0)
            return;
        
        // get the id of ancestor for operating on the node's FVal vector aligned with ancestors
        int ancestorID =  getAncestorId(f, ancestor);
//        std::cout << "ancestor id and label " << ancestorID << " - " << f->ancestors[ancestorID]->factorLabel << std::endl;
//        std::cout << "include " << f->Finclude[ancestorID][solutionSize-1] << std::endl;
//        cout.precision(dbl::max_digits10);
//        std::cout << "fval " << f->FVal[ancestorID][solutionSize-1] << std::endl;
        
        // if this node is included in the solution
        if (f->Finclude[ancestorID][solutionSize-1]) {
//            indexFactors.push_back(f);
            f->isIndex = true;
            std::cout << "index: " <<  f->factorId << std::endl;
//            std::cout << "index: " <<  f->factorId << " is dummy?? "<< f->isDummy << std::endl;
            solutionSize--;
//            std::cout << "solutionSize after icnld" <<  solutionSize << std::endl;
            if (solutionSize == 0)
                return;
            
            if (f->binaryChildFactors.size() == 1) {
//                std::cout << "dp 1 " << std::endl;
//                std::cout << "mp1: " << solutionSize << std::endl;
                constructDPSolution(f->binaryChildFactors[0], f, solutionSize);
            }
            else {
//                std::cout << "dp 2 " << std::endl;
                // get children's solution sizes for this node being their solution ancestor
                std::vector<int> maximizerPartition;
                maximizerPartition.reserve(f->binaryChildFactors.size());
                getMaximizerPartition(f, 0, maximizerPartition, solutionSize);
//                std::cout << "mp1: " << maximizerPartition << std::endl;
//                std::cout << "maximizerPartition" << maximizerPartition.size() << std::endl;
                // recurse on children with this node being their solution ancestor
                for (int i = 0; i < f->binaryChildFactors.size(); i++) {
                    constructDPSolution(f->binaryChildFactors[i], f, maximizerPartition[i]);
                }
            }

        }
        
        // if this node is not included in the solution
        else {
            // get children's solution sizes for this node's solution ancestor being their solution ancestor
//            std::cout << "else " << std::endl;
            if (f->binaryChildFactors.size() == 1) {
//                std::cout << "dp 1 " << std::endl;
//                std::cout << "mp2: " << solutionSize << std::endl;
                constructDPSolution(f->binaryChildFactors[0], ancestor, solutionSize);
            }
            else if (f->binaryChildFactors.size() > 1) {
//            else {
                std::vector<int> maximizerPartition;
                ++ancestorID;
                getMaximizerPartition(f, ancestorID, maximizerPartition, solutionSize);
//                std::cout << "mp2: " << maximizerPartition << std::endl;
//                std::cout << "solutionSize and ancestor name " <<  solutionSize << " - " << f->ancestors[ancestorID]->factorLabel << std::endl;
                // recurse on children with this node's solution ancestor being their solution ancestor
                for (int i = 0; i < f->binaryChildFactors.size(); i++) {
                    constructDPSolution(f->binaryChildFactors[i], ancestor, maximizerPartition[i]); // recursion using this node's ancestor
                }
                
            }

        }
    }
    
    
    
    
    
} //




