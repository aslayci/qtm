#include "factorTree.h"

namespace qtm {

    
    // constructor for queryprocessing
    FactorTree::FactorTree(BayesianNetwork *myNet, std::vector<int> &eliminationOrder, Query &query, string elimOrderType, int k){
        
        deleteThings = true;
        
        forDP = false;
        elimOrderOption = elimOrderType;
        
        this->queryVariables = query.queryVars;
        this->evidenceStates = query.evidenceStates;
        
        for (int i = 0; i < eliminationOrder.size(); i++) {
            elimOrderMap.push_back(std::vector<Factor*>());
        }
        
        nrFactors = 0;
        nrSecondaryFactors = 0;
        
        indexSize = k;
        if (indexSize > 0) {
            indexFactorIds.reserve(indexSize);
            folderName = std::string("index_files_") + elimOrderType + std::string("_") + intToStr(indexSize);
            readIndexFactorIds(folderName, indexFactorIds);
        }
        
        nrUseful = 0;
        nrSkipped = 0;
        nrCreated = 0;
        totalObtainedBenefit = 0;
        totalComputationCost = 0;
        
        createInitialFactorsFromCPTs(myNet, eliminationOrder);
        buildFactorEliminationTree(myNet, eliminationOrder);
    }
    
    
    // constructor for preprocessing
    FactorTree::FactorTree(BayesianNetwork *myNet, std::vector<int> &eliminationOrder, string elimOrderType){
        
        deleteThings = true;
        
        forDP = true;
        elimOrderOption = elimOrderType;
    
        for (int i = 0; i < eliminationOrder.size(); i++) {
            elimOrderMap.push_back(std::vector<Factor*>());
        }
        
        nrFactors = 0;
        nrSecondaryFactors = 0;
        
        totalObtainedBenefit = 0;
        totalComputationCost = 0;
        
        createInitialFactorsFromCPTs(myNet, eliminationOrder);
        buildFactorEliminationTree(myNet, eliminationOrder);
    }
    
    FactorTree::~FactorTree(){}
    
    // creates organic factors from CPTs
    void FactorTree::createInitialFactorsFromCPTs(BayesianNetwork *myNet, std::vector<int> &eliminationOrder) {
        
        for (int i = 0; i < myNet->nodesVector.size(); i++) {
            
            if(myNet->nodesVector[i]->isRoot()) { // if this BN node is a root in BN
                Factor *f = new Factor(myNet);
                nrFactors++;
//                f->factorLabel = "P(" + myNet->nodesVector[i]->getNodeLabel() + ")";
                f->isLeaf = true;
                f->scope.push_back(i);
                f->var_to_index[i] = 0;
                f->assignDomainInfo();
                if (forDP) {
                    f->usefulAloneProb = 1.0;
                    f->probSubtreeInQuery = 0;
                }
                f->factorTableLength = myNet->nodesVector[i]->nrStates;
                f->cost = f->factorTableLength;
                f->utility = f->factorTableLength;
//                f->utility = 0;
                int orderIndex = f->getElimParentOrderId(eliminationOrder, 0); //inference
                f->factorId = myNet->nodesVector.size() + i;
                f->factorsInSubtree.insert(f->factorId);
                elimOrderMap[orderIndex].push_back(f);
            }
            
            else {
                Factor *f = new Factor(myNet);
                nrFactors++;
                f->isLeaf = true;
                f->scope.push_back(i); // first node in scope corresponds to the child in this CPT
                f->var_to_index[i] = 0;
                f->factorTableLength = myNet->nodesVector[i]->nrStates;
//                std::string ss = "P(" + myNet->nodesVector[i]->getNodeLabel() + " | ";
                for (int j = 0; j < myNet->nodesVector[i]->parents.size(); j++) {
                    int parentId = myNet->nodesVector[i]->parents[j];
//                    ss += myNet->nodesVector[parentId]->getNodeLabel() + ", ";
                    f->scope.push_back(parentId); // scope is aligned with the order of parents read from input file
                    f->var_to_index[parentId] = j+1;
                    f->factorTableLength *= myNet->nodesVector[parentId]->nrStates;
                }
//                ss += ")";
//                f->factorLabel = ss;
                if (forDP) {
                    f->usefulAloneProb = 1.0;
                    f->probSubtreeInQuery = 0;
                }
                f->cost = f->factorTableLength;
                f->utility = f->factorTableLength;
                f->assignDomainInfo();
                int orderIndex = f->getElimParentOrderId(eliminationOrder, 0);
                f->factorId = myNet->nodesVector.size() + i;
                f->factorsInSubtree.insert(f->factorId);
                elimOrderMap[orderIndex].push_back(f);
            }
        }
    }
    
    // creates the secondary factors and the elimination tree
    void FactorTree::buildFactorEliminationTree(BayesianNetwork *myNet, std::vector<int> &eliminationOrder) {
        
        for (int i = 0; i < eliminationOrder.size(); i++) { // create secondary factors
            
            Factor *f = new Factor(myNet);
            nrFactors++;
            nrSecondaryFactors++;
            
            if(forDP) {
                f->var_to_sum_out = eliminationOrder[i];
            }
            else {
                if(std::find(queryVariables.begin(), queryVariables.end(), eliminationOrder[i]) == queryVariables.end()) { // if not a query var
                    f->var_to_sum_out = eliminationOrder[i];
                }
            }
            
            f->factorId = eliminationOrder[i];
            f->isLeaf = false;
            
            double forUtility = 0.0;
            f->factorTableLength = 1;
            f->subtreeSize = 1; // counts only sec. factors in size for DP
            int nrFactorsSubtree = 1; // counts all factors for inference
//            std::string ss = "";
            int scpCnter = 0;
            f->probSubtreeInQuery = 0;
            // this will be updated below with children's usefulAloneProbs
//            if (forDP) {
//                f->probSubtreeInQuery = myNet->nodesVector[eliminationOrder[i]]->probInQuery;
//            }
            
            for (int j = 0; j < elimOrderMap[i].size(); j++) { // iterating over child factors
                
                f->childFactors.push_back(elimOrderMap[i][j]); // assign the pointers to its child factors here
                
                if (forDP) {
                    if(!elimOrderMap[i][j]->isLeaf) { //assign the pointers to its secondary-factor type children
//                        f->probSubtreeInQuery += elimOrderMap[i][j]->probSubtreeInQuery;
                        f->secChildFactors.push_back(elimOrderMap[i][j]);
                        f->subtreeSize += elimOrderMap[i][j]->subtreeSize;
                    }
                }
                
                nrFactorsSubtree += elimOrderMap[i][j]->factorsInSubtree.size();
                
                forUtility += elimOrderMap[i][j]->utility; // sum utility of child factors'
                
                // set the scope of this factor from children's scopes
                for (int z = 0; z < elimOrderMap[i][j]->scope.size(); z++) {
                    // var_to_sum_out is not added to the scope
                    // i.e., if this node is a query variable, then it is added to the scope as its var_to_sum_out = -1
                    int varZ = elimOrderMap[i][j]->scope[z];
                    if(varZ == f->var_to_sum_out)
                        continue;
                    
                    if(!(f->inScope(varZ))) {
                        f->scope.push_back(varZ);
                        f->var_to_index[varZ] = scpCnter;
                        ++scpCnter;
                        f->factorTableLength *= myNet->nodesVector[varZ]->nrStates; // for size after summing out
                        
//                        ss +=  myNet->nodesVector[varZ]->getNodeLabel() + ", "; // for printing pretty
//                        ss +=  intToStr(varZ) + ", "; // for printing pretty
                    }
                }
                
            } // end of elimordermap[i] traversal
            
            f->assignDomainInfo();
      
            f->factorsInSubtree.reserve(nrFactorsSubtree);
            
            f->factorsInSubtree.insert(f->factorId);
            
            for (int j = 0; j < f->childFactors.size() ; j++) {
                f->factorsInSubtree.insert(f->childFactors[j]->factorsInSubtree.begin(), f->childFactors[j]->factorsInSubtree.end());
            }

            f->factorLabel = "G_"+ intToStr(f->factorId); // shorter version with var ids and no scope
//            f->factorLabel = "G_"+ intToStr(f->factorId) + "(" + ss + ")";; // shorter version with var ids
//            f->factorLabel = "G_"+ myNet->nodesVector[eliminationOrder[i]]->getNodeLabel(); // shorter version wo scope info
            //f->factorLabel = "G_"+ myNet->nodesVector[eliminationOrder[i]]->getNodeLabel() + "(" + ss + ")"; // for printing pretty
            
            
            if (forDP) {
                f->cost = 2 * (f->factorTableLength * myNet->nodesVector[eliminationOrder[i]]->nrStates);
                f->utility = forUtility + f->cost;
            }
            
            else if (!forDP && f->var_to_sum_out == myNet->nodesVector.size()) { // no var to sum out since query var
                f->cost = f->factorTableLength; // table already contains this query node
                f->utility = forUtility + f->cost;
            }
            
            else if (!forDP && f->var_to_sum_out != myNet->nodesVector.size()) {
                f->cost = 2 * (f->factorTableLength * myNet->nodesVector[eliminationOrder[i]]->nrStates);
                f->utility = forUtility + f->cost;
            }
            
            bool forSecondaryLeaf = true;
            for (int  j = 0; j < f->childFactors.size(); j++) {
                if (!(f->childFactors[j]->isLeaf)) {
                    forSecondaryLeaf = false;
                    break;
                }
            }
           if (forSecondaryLeaf)
                f->isSecondaryLeaf = true;
  
            int orderIndex = f->getElimParentOrderId(eliminationOrder, i + 1); // to check if this factor is root or not
            if(orderIndex != eliminationOrder.size()) { // if this factor is not the root
                elimOrderMap[orderIndex].push_back(f);
            }
            else {
                f->isRoot = true;
                f->usefulAloneProb = 0;
                f->utility = 0;
                f->factorTableLength = 0;
                rootFactor = f;
            }
            
        } // end of elimordermap traversal
        
        // assign the postOrder
        postOrder.reserve(nrFactors);
        getPostorder(rootFactor, postOrder);
        
        // load useful factor tables
        if (!forDP && indexSize > 0) {
            for (int i = 0; i < indexFactorIds.size(); i++) {
                postOrder[factor_to_poPosition[indexFactorIds[i]]]->isIndex = true;
            }

            findUsefulIndexFactors(rootFactor);
        }
        
        // compute and write tree stats -- executed only for preprocessing
        if (forDP) {
            // create epsilon ancestor to represent no ancestor
            Factor *eps = new Factor(myNet);
            eps->factorLabel = "epsilon";
            eps->usefulAloneProb = 0;
            eps->factorId = -1;
            
            // factor table length: nrLines
            int maxNrChildren = 0, minNrChildren =  myNet->nrNodes;
            double avgNrChildren = 0.0;
            uint64 minFactorTableLength = LLONG_MAX, maxFactorTableLength = 0, avgFactorTableLength = 0 ;
            uint64 minUnrolledLength = LLONG_MAX, maxUnrolledLength = 0, avgUnrolledLength = 0 ;
            
            int nrSecondaryLeaves = 0;
            
            // binary conversion operations for DP
            splitToBinary(myNet, rootFactor, rootFactor->secChildFactors);
            getBinaryPostorder(rootFactor, binaryPostOrder);
//            std::cout << "binary po size " << binaryPostOrder.size() << std::endl;
//            exit(1);
            
            for (int i = 0; i < binaryPostOrder.size(); i++) {
                Factor *f = binaryPostOrder[i];
                // set binary ancestors for new DP
                getBinaryAncestors(rootFactor, f, f->ancestors);
                // getAncestors(rootFactor, f, f->ancestors);
                // add epsilon ancestor representing no ancestor
                f->ancestors.push_back(eps); // needed for DP
            }
            
            // post order and some stats
            for (int i = 0; i < postOrder.size(); i++) {
                
                Factor *f = postOrder[i];
//                std::cout << i << " - " << f->factorLabel << std::endl;
                
                // for tree stats
                if(f->isLeaf)
                    continue;
                
                if(f->isSecondaryLeaf)
                    nrSecondaryLeaves++;
                
//                avgNrChildren += f->secChildFactors.size();
//                if(maxNrChildren < f->secChildFactors.size())
//                    maxNrChildren = f->secChildFactors.size();
//                if(minNrChildren > f->secChildFactors.size())
//                    minNrChildren = f->secChildFactors.size();
                
                avgNrChildren += f->childFactors.size();
                if(maxNrChildren < f->childFactors.size())
                    maxNrChildren = f->childFactors.size();
                if(minNrChildren > f->childFactors.size())
                    minNrChildren = f->childFactors.size();
                
                if (f->isRoot)
                    continue;
                
                // don't take root into account for below stats since it has zero lines after summing out when there is no query node
                avgFactorTableLength += f->factorTableLength;
                if (maxFactorTableLength < f->factorTableLength)
                    maxFactorTableLength = f->factorTableLength;
                if (minFactorTableLength > f->factorTableLength)
                    minFactorTableLength = f->factorTableLength;
                
                uint64 unrolledLength = f->cost / 2;
                avgUnrolledLength += unrolledLength;
                if (maxUnrolledLength < unrolledLength)
                    maxUnrolledLength = unrolledLength;
                if (minUnrolledLength > unrolledLength)
                    minUnrolledLength = unrolledLength;
                
//                std::cout << "f table length for factor " << f->factorLabel << " : " << f->factorTableLength << std::endl;
//                if (f->factorTableLength == -1) {
//                    std::cout << "here -1 length: " << f->factorId << std::endl;
//                }
               
//                forFTSize = f->factorTableLength * f->scope.size();
                
//                avgFactorTableSize += forFTSize;
//                if (maxFactorTableSize < forFTSize)
//                    maxFactorTableSize = forFTSize;
//                if (minFactorTableSize > forFTSize) {
//                    minFactorTableSize = forFTSize;
//                }
            }
            
//            avgNrChildren /= nrSecondaryFactors;
            avgNrChildren /= nrFactors;
            avgFactorTableLength /= nrSecondaryFactors;
            avgUnrolledLength /= nrSecondaryFactors;
            
            int treeHeight = height(rootFactor);
            int secondaryTreeHeight = secondaryHeight(rootFactor); //without organic leaves
            int binaryTreeHeight = binaryHeight(rootFactor) + 1; // plus one for organic leaves
            
            // write these stats to file
            string command1 = string("mkdir -p treeStats");
            system(command1.c_str());
            
            string filename = std::string("treeStats") + OS_SEP + std::string("treeStats_") + elimOrderOption + std::string(".txt");
            ofstream outMasterStream;
            
            if(outMasterStream.is_open())
                outMasterStream.close();
            
            outMasterStream.open(filename);
            
            if (outMasterStream.is_open() == false) {
                cerr << "ERROR: cannot open " << filename << " for writing " << endl;
                exit(1);
            }
            
            outMasterStream << "total number of nodes in elimination tree : " << nrFactors << "\n";
            outMasterStream << "there are " << nrSecondaryFactors << " secondary factors " << " and " << nrSecondaryLeaves << " secondary leaves " << "\n";
            outMasterStream << "min & max & avg nr. of children among secondary factors: " << minNrChildren << " & " << maxNrChildren << " & " << avgNrChildren << "\n";
            outMasterStream << "min & max & avg factorTableLength among secondary factors: " << minFactorTableLength << " & " << maxFactorTableLength << " & " << avgFactorTableLength << "\n";
            outMasterStream << "min & max & avg unrolled factorTableLength among secondary factors: " << minUnrolledLength << " & " << maxUnrolledLength << " & " << avgUnrolledLength << "\n";
            outMasterStream << "height of tree : " << treeHeight << "\n";
            outMasterStream << "secondary height of tree : " << secondaryTreeHeight << "\n";
            outMasterStream << "binary height of tree : " << binaryTreeHeight << "\n";
            
            outMasterStream.close();
        }
    }
    
    // just to print to screen for controls
    void FactorTree::printPostorder(Factor *f) {
        
        if (f == NULL)
            return;
        
        for (int i = 0; i < f->childFactors.size(); i++) {
            printPostorder(f->childFactors[i]);
        }
        
        std::cout << f->factorLabel << " - ";
    }
    
    void FactorTree::printPretty(Factor *f, std::string indent, bool last) {
        std::cout << indent;
        if (last) {
            std::cout << "\\-";
            indent += "  ";
        }
        else {
            std::cout << "|-";
            indent += "| ";
        }
        std::cout << f->factorLabel << std::endl;
        
        for (int i = 0; i < f->childFactors.size(); i++) {
            printPretty(f->childFactors[i], indent, i == f->childFactors.size() - 1);
        }
    }
    
    void FactorTree::printPrettySecondary(Factor *f, std::string indent, bool last) {
        std::cout << indent;
        if (last) {
            std::cout << "\\-";
            indent += "  ";
        }
        else {
            std::cout << "|-";
            indent += "| ";
        }
        std::cout << f->factorLabel << std::endl;
        
        for (int i = 0; i < f->secChildFactors.size(); i++) {
            printPrettySecondary(f->secChildFactors[i], indent, i == f->secChildFactors.size() - 1);
        }
    }
    
    
    void FactorTree::printPrettySecWithLenthRatios(Factor *f, std::string indent, bool last) {
        std::cout << indent;
        if (last) {
            std::cout << "\\-";
            indent += "  ";
        }
        else {
            std::cout << "|-";
            indent += "| ";
        }
        std::cout << f->factorLabel << std::endl;
        
        for (int i = 0; i < f->secChildFactors.size(); i++) {
            printPrettySecondary(f->secChildFactors[i], indent, i == f->secChildFactors.size() - 1);
        }
    }
    
    
    void FactorTree::printPrettyBinary(Factor *f, std::string indent, bool last) {
        std::cout << indent;
        if (last) {
            std::cout << "\\-";
            indent += "  ";
        }
        else {
            std::cout << "|-";
            indent += "| ";
        }
        std::cout << f->factorLabel << std::endl;
        
        for (int i = 0; i < f->binaryChildFactors.size(); i++) {
            printPrettyBinary(f->binaryChildFactors[i], indent, i == f->binaryChildFactors.size() - 1);
        }
    }
    
    
    int FactorTree::height(Factor *f) {
        
        if (f == NULL)
            return 0;
        
        else {
            std::vector<int> heights;
            for (int i = 0; i < f->childFactors.size(); i++) {
                heights.push_back(height(f->childFactors[i]));
            }
            
            int maxHeight = 0;
            for (int i = 0; i < heights.size(); i++) {
                if(maxHeight < heights[i])
                    maxHeight = heights[i];
            }
            return maxHeight + 1;
        }
    }
    
    int FactorTree::secondaryHeight(Factor *f) {
        
        if (f == NULL)
            return 0;
        
        else {
            std::vector<int> heights;
            for (int i = 0; i < f->secChildFactors.size(); i++) {
                heights.push_back(secondaryHeight(f->secChildFactors[i]));
            }
            
            int maxHeight = 0;
            for (int i = 0; i < heights.size(); i++) {
                if(maxHeight < heights[i])
                    maxHeight = heights[i];
            }
            return maxHeight + 1;
        }
    }
    
    int FactorTree::binaryHeight(Factor *f) {
        
        if (f == NULL)
            return 0;
        
        else {
            std::vector<int> heights;
            for (int i = 0; i < f->binaryChildFactors.size(); i++) {
                heights.push_back(binaryHeight(f->binaryChildFactors[i]));
            }
            
            int maxHeight = 0;
            for (int i = 0; i < heights.size(); i++) {
                if(maxHeight < heights[i])
                    maxHeight = heights[i];
            }
            return maxHeight + 1;
        }
    }
    
    // level order for secondary tree
    void FactorTree::getLevelOrder(Factor *f, std::vector< std::vector<Factor*>> &levels) {
        int h = secondaryHeight(f);
        for (int i = 1; i <= h; i++) {
            std::vector<Factor*> levelOrder;
            getGivenLevel(f, i, levelOrder);
            levels.push_back(levelOrder);
        }
    }
    // for level order of secondary tree
    void FactorTree::getGivenLevel(Factor *f, int level, std::vector<Factor*> &levelOrder) {
        
        if (f == NULL)
            return;
        
        if (level == 1) {
            levelOrder.push_back(f);
            
        }
        
        else if (level > 1) {
            for (int i = 0; i < f->secChildFactors.size(); i++) {
                getGivenLevel(f->secChildFactors[i], level - 1, levelOrder);
            }
            
        }
    }
    
    // level order for secondary tree
    void FactorTree::getBinaryLevelOrder(Factor *f, std::vector< std::vector<Factor*>> &levels) {
        int h = binaryHeight(f);
        for (int i = 1; i <= h; i++) {
            std::vector<Factor*> levelOrder;
            getBinaryGivenLevel(f, i, levelOrder);
            levels.push_back(levelOrder);
        }
    }
    // for level order of secondary tree
    void FactorTree::getBinaryGivenLevel(Factor *f, int level, std::vector<Factor*> &levelOrder) {
        
        if (f == NULL)
            return;
        
        if (level == 1) {
            levelOrder.push_back(f);
            
        }
        
        else if (level > 1) {
            for (int i = 0; i < f->binaryChildFactors.size(); i++) {
                getGivenLevel(f->binaryChildFactors[i], level - 1, levelOrder);
            }
            
        }
    }

    
//    void FactorTree::getLevelOrder(Factor *f, std::vector<Factor*> &levelOrder) {
//        int h = height(f);
//        for (int i = 1; i <= h; i++) {
//            getGivenLevel(f, i, levelOrder);
//        }
//    }
//
//    void FactorTree::getGivenLevel(Factor *f, int level, std::vector<Factor*> &levelOrder) {
//
//        if (f == NULL)
//            return;
//
//        if (level == 1) {
//            levelOrder.push_back(f);
//        }
//
//        else if (level > 1) {
//            for (int i = 0; i < f->childFactors.size(); i++) {
//                getGivenLevel(f->childFactors[i], level - 1, levelOrder);
//            }
//
//        }
//    }
    
    void FactorTree::getPostorder(Factor *f, std::vector<Factor*> &postOrder) {
        
        if (f == NULL)
            return;
        
        for (int i = 0; i < f->childFactors.size(); i++) {
            getPostorder(f->childFactors[i], postOrder);
        }
        
        postOrder.push_back(f);
        
        factor_to_poPosition[f->factorId] = postOrder.size()-1;
        
    }
    
    void FactorTree::getBinaryPostorder(Factor *f, std::vector<Factor*> &binaryPostOrder) {
        
        if (f == NULL)
            return;
        
        for (int i = 0; i < f->binaryChildFactors.size(); i++) {
            getBinaryPostorder(f->binaryChildFactors[i], binaryPostOrder);
        }
        
        binaryPostOrder.push_back(f);
    }
    
    
    void FactorTree::splitToBinary(BayesianNetwork *myNet, Factor *root, vector<Factor*> remainingSecChildren) {
        
        if (root->isSecondaryLeaf) {
            return;
        }
        
        if (remainingSecChildren.size() == 1) { // to account for factors that have only one secondary child
            root->binaryChildFactors = remainingSecChildren;
            splitToBinary(myNet, root->binaryChildFactors[0], root->binaryChildFactors[0]->secChildFactors);
        }
        
        else if (remainingSecChildren.size() == 2) {
            root->binaryChildFactors = remainingSecChildren;
            // start over from scratch for these 2
            splitToBinary(myNet, root->binaryChildFactors[0], root->binaryChildFactors[0]->secChildFactors);
            splitToBinary(myNet, root->binaryChildFactors[1], root->binaryChildFactors[1]->secChildFactors);
        }
        
        else if(remainingSecChildren.size() > 2) {
            root->binaryChildFactors.reserve(2);
            Factor *assignedRealChild = remainingSecChildren.back();
            root->binaryChildFactors.push_back(assignedRealChild); // assign a real child as binary child
            // // call from scratch for the assigned node
            splitToBinary(myNet, assignedRealChild, assignedRealChild->secChildFactors);
            remainingSecChildren.pop_back(); // remove it from remaining children
            
            Factor *f_dummy = new Factor(myNet);  // dummy child as binary-sibling of assignedRealChild
            f_dummy->factorLabel = "dummy";
            f_dummy->isDummy = true;
            f_dummy->utility = 0;
//            f_dummy->usefulAloneProb = 0; // setting 0 is wrong, instead returning probUseful = 0 wrt this ancestor
            f_dummy->factorId = assignedRealChild->factorId + 100000;
            if (root->isDummy)  // root = dummy factor case
                f_dummy->subtreeSize = root->subtreeSize - assignedRealChild->subtreeSize;
            
            else // root = secondary factor case
                f_dummy->subtreeSize = root->subtreeSize - 1 - assignedRealChild->subtreeSize ; // - 1 is not to count the root node as its subtree size contains itself
            
            root->binaryChildFactors.push_back(f_dummy);
            // recurse on f_dummy
            splitToBinary(myNet,f_dummy, remainingSecChildren);
        }
        
        // call it for assigned binary nodes children
            

//        // we don't need to care about organicChildren
//        int nrSecChildren = root->secChildFactors.size();
        
        
    }
    
    bool FactorTree::getAncestors(Factor *root, Factor *f, std::vector<Factor*> &ancestors) {
        
        if (root == NULL)
            return false;
        
        if (root == f)
            return true;
        
        // if target is present in either of the subtrees of the children of this node, then print this node
        for (int i = 0; i < root->childFactors.size(); i++) {
            if(getAncestors(root->childFactors[i], f, ancestors)) {
                ancestors.push_back(root);
//                std::cout << "printing anc: " << root->factorLabel << std::endl;
                return true;
            }
        }

        return false;
    }

    // only for DP
    bool FactorTree::getBinaryAncestors(Factor *root, Factor *f, std::vector<Factor*> &ancestors) {
        
        if (root == NULL)
            return false;
        
        if (root == f)
            return true;
        
        // if target is present in either of the subtrees of the children of this node, then print this node
        for (int i = 0; i < root->binaryChildFactors.size(); i++) {
            if(getBinaryAncestors(root->binaryChildFactors[i], f, ancestors)) {
                ancestors.push_back(root);
                //                std::cout << "printing anc: " << root->factorLabel << std::endl;
                return true;
            }
        }
        
        return false;
    }

    
    void FactorTree::printLevelOrder(Factor *f) {
        int h = height(f);
        for (int i = 1; i <= h; i++) {
            printGivenLevel(f, i);
            std::cout << std::endl;
        }
    }
    
    void FactorTree::printGivenLevel(Factor *f, int level) {
        if (f == NULL)
            return;
        
        if (level == 1) {
            std::cout << f->factorLabel << " ";
//            std::cout << "cost, utility, prob, table size, subtree size for " << f->factorLabel << ": " << f->cost << " & " << f->utility << " & " << f->usefulAloneProb << " & " << f->factorTableLength << " & " << f->subtreeSize << std::endl;
        }
        
        else if (level > 1) {
            for (int i = 0; i < f->childFactors.size(); i++) {
                printGivenLevel(f->childFactors[i], level - 1);
            }
        }
    }
    
    bool FactorTree::plainVariableElimination(BayesianNetwork *myNet, long double &tableCreationTimes) {
        
        bool withinTimeLimit = true;
        
        for (int i = 0; i < postOrder.size(); i++) {
            // cide: rt disina cikar sonra
            if (postOrder[i]->isLeaf) {
                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createOrganicFactorTable();
                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
                totalComputationCost += postOrder[i]->cost;
                ++nrCreated;
                tableCreationTimes += getElapsedWallTime(create_begin, create_end);
            }
            else{
                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createFactorTableFromSumProduct();
                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
                totalComputationCost += postOrder[i]->cost;
                ++nrCreated;
                tableCreationTimes += getElapsedWallTime(create_begin, create_end);
            }
            
            if (tableCreationTimes > 1200.0) {
                withinTimeLimit = false;
                break;
            }
            
            if (deleteThings) {
                for (int j = 0; j < postOrder[i]->childFactors.size(); j++) {
                    postOrder[i]->childFactors[j]->factorTableProbs.clear();
                    postOrder[i]->childFactors[j]->factorTableProbs.shrink_to_fit();
                }
            }
        }
        
        return withinTimeLimit;
        
    }
    
    bool FactorTree::indexedVariableElimination(BayesianNetwork *myNet, long double &tableCreationTimes) {
        
        bool withinTimeLimit = true;

        for (int i = 0; i < postOrder.size(); i++) {
            
            if (postOrder[i]->canSkip) {
                ++nrSkipped;
                totalObtainedBenefit += postOrder[i]->cost;
                totalComputationCost += postOrder[i]->cost;
                continue;
            }
            if (postOrder[i]->isUseful) {
                postOrder[i]->readFactorTablesFromIndex(folderName);
                ++nrUseful;
                totalObtainedBenefit += postOrder[i]->cost;
                totalComputationCost += postOrder[i]->cost;
                continue;
            }
            
            if (postOrder[i]->isLeaf) {
//                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createOrganicFactorTable();
//                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
//                tableCreationTimes += getElapsedWallTime(create_begin, create_end);
                totalComputationCost += postOrder[i]->cost;
                ++nrCreated;
            }
            else{
//                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createFactorTableFromSumProduct();
//                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
//                tableCreationTimes += getElapsedWallTime(create_begin, create_end);
                totalComputationCost += postOrder[i]->cost;
                ++nrCreated;
            }
            
//            if (tableCreationTimes > 1200.0) {
//                withinTimeLimit = false;
//                break;
//            }
//
//            if (deleteThings) {
//                for (int j = 0; j < postOrder[i]->childFactors.size(); j++) {
//                    postOrder[i]->childFactors[j]->factorTableProbs.clear();
//                    postOrder[i]->childFactors[j]->factorTableProbs.shrink_to_fit();
//                }
//            }
        }
        
        return withinTimeLimit;
        
    }
    
    void FactorTree::findUsefulIndexFactors(Factor *f) {
        
        if (f->isLeaf)
            return;
        
        if (!f->isIndex) { // if it is not an index factor
            for (int i = 0; i < f->childFactors.size(); i++) {
                findUsefulIndexFactors(f->childFactors[i]);
            }
        }
        
        else {
            bool isUseful = true;
            // check whether its subtree contains any query variable
            for (int i = 0; i < queryVariables.size(); i++) {
                if (f->factorsInSubtree.find(queryVariables[i]) != f->factorsInSubtree.end()) {
                    isUseful = false;
                    break;
                }
            }
            if (isUseful) {
                f->isUseful = true; // this is assigned true only to index factors which are useful
                //                std::cout << "subtree of useful " << f->factorLabel << std::endl;
                std::unordered_set<int>::iterator it;
                for (it = f->factorsInSubtree.begin(); it != f->factorsInSubtree.end(); it++) {
                    if (*it == f->factorId)
                        continue;
                    postOrder[factor_to_poPosition[*it]]->canSkip = true;
                }
                return;
            }
            else {
                for (int i = 0; i < f->childFactors.size(); i++) {
                    findUsefulIndexFactors(f->childFactors[i]);
                }
            }
        }
        
    }
    
    void FactorTree::readIndexFactorIds(string folderName, std::vector<int> &indexFactorIds) {
        
        string filename = folderName + OS_SEP + "indexIDs.txt";
        ifstream input_file(filename);
        
        if (input_file.is_open()) {
            std::string line;
            getline(input_file, line);
            splitint(line, ",", indexFactorIds);
            input_file.close();
        }
        
        else {
            cerr << "ERROR: cannot open " << filename << " for reading " << endl;
            exit(1);
        }
    }
    
    void FactorTree::outputCostsForCorrelationAnalysis(BayesianNetwork *myNet) {
        
        string costsFilename = std::string("treeStats") + OS_SEP + std::string("costs.txt");
        
        // output also costs for the correlation analysis
        ofstream outMasterStream2;
        
        if(outMasterStream2.is_open())
            outMasterStream2.close();
        
        outMasterStream2.open(costsFilename);
        
        if (outMasterStream2.is_open() == false) {
            cerr << "ERROR: cannot open " << costsFilename << " for writing " << endl;
            exit(1);
        }
        
        for (int i = 0; i < postOrder.size(); i++) {
            
            if (postOrder[i]->isLeaf) {
                postOrder[i]->createOrganicFactorTable();
            }
            
            else {
                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createFactorTableFromSumProduct();
                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
                postOrder[i]->tableCreationTime = getElapsedWallTime(create_begin, create_end);
                
                double cost1 = 0, cost4 = 0;
                double cost2 = postOrder[i]->factorTableLength * myNet->nodesVector[postOrder[i]->var_to_sum_out]->nrStates; // size b4 summing out
                double cost3 = 2 * postOrder[i]->factorTableLength * myNet->nodesVector[postOrder[i]->var_to_sum_out]->nrStates;
                
                // compute different type of costs here
                for (int j = 0; j < postOrder[i]->childFactors.size(); j++) {
                    
                    cost1 += postOrder[i]->childFactors[j]->factorTableLength;
                    cost2 += postOrder[i]->childFactors[j]->factorTableLength;
                    cost4 += postOrder[i]->childFactors[j]->factorTableLength * log(1 + postOrder[i]->childFactors[j]->factorTableLength);
                }
                
                outMasterStream2 << std::fixed << std::setprecision(5) << postOrder[i]->tableCreationTime << " " << cost1 << " " << std::fixed << std::setprecision(0) << cost2 << " " << cost3 << " " << std::fixed << std::setprecision(0) << cost4 << "\n";
                
            }
            
            if (deleteThings) {
                for (int j = 0; j < postOrder[i]->childFactors.size(); j++) {
                    postOrder[i]->childFactors[j]->factorTableProbs.clear();
                    postOrder[i]->childFactors[j]->factorTableProbs.shrink_to_fit();
                }
            }
        }
        
        outMasterStream2.close();
        
    }

    
    void FactorTree::fillAndWriteIndexTables(string folderName, int indexSize) {
        
        // if index_files_folder doesn't exist, create it
        string command1 = string("mkdir -p ") + folderName;
        system(command1.c_str());
        string command2 = string("rm -r ") + folderName; // remove older index files from dir if exists any
        system(command2.c_str());
        system(command1.c_str());
        
        std::vector<int> indexFactors;
        
        for (int i = 0; i < postOrder.size(); i++) {
            
            // for now comment this for measuring the costs correlation with table creation time
            if (indexFactors.size() == indexSize) // no more index tables left to write so can stop creation
                break;
            
            if (postOrder[i]->isLeaf) {
                postOrder[i]->createOrganicFactorTable();
            }

            else {
                std::chrono::steady_clock::time_point create_begin = std::chrono::steady_clock::now();
                postOrder[i]->createFactorTableFromSumProduct();
                std::chrono::steady_clock::time_point create_end = std::chrono::steady_clock::now();
                postOrder[i]->tableCreationTime = getElapsedWallTime(create_begin, create_end);
            }
            
            // write the table if it is an index factor
            if (postOrder[i]->isIndex) {
                postOrder[i]->writeIndexFactorTables(folderName);
                indexFactors.push_back(postOrder[i]->factorId);
            }

//            if (deleteThings) {
//                for (int j = 0; j < postOrder[i]->childFactors.size(); j++) {
//                    postOrder[i]->childFactors[j]->factorTableProbs.clear();
//                    postOrder[i]->childFactors[j]->factorTableProbs.shrink_to_fit();
//                }
//            }
            
        }
        
        // write a file that contains the id of index factors only -- so that we only import tables of useful factors during inference
        string filename = folderName + OS_SEP + std::string("indexIDs.txt");
        ofstream outMasterStream;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(filename);
        
        if (outMasterStream.is_open() == false) {
            cerr << "ERROR: cannot open " << filename << " for writing " << endl;
            exit(1);
        }
        
        for (int i = 0; i < indexFactors.size(); i++) {
            if (i < indexFactors.size() - 1)
                outMasterStream << indexFactors[i] << ",";
            else
                outMasterStream << indexFactors[i] << "\n";
        }
        
        outMasterStream.close();
    }
    
    void FactorTree::assignUsefulnessProbs(BayesianNetwork *myNet, string workloadType) {
        
        if (workloadType.compare("unif") == 0) {
            assignUniformProbNodeInQuery(myNet);
        }
        else if (workloadType.compare("biasr") == 0) { // bias towards root
            assignRootBiasedProbNodeInQuery(myNet);
        }
//        else if (workloadType.compare("biasl") == 0) { // bias towards leaves
//
//        }
        
        else {
            std::cout << "unrecognized workload type option, exiting." << std::endl;
            exit(1);
        }
        
        //
        for (int i = 0; i < postOrder.size(); i++) {
            
            Factor *f = postOrder[i];
            
            if (f->isLeaf)
                continue;
            
            std::unordered_set<int>::iterator it;
            for (it = f->factorsInSubtree.begin(); it != f->factorsInSubtree.end(); it++) {
                if (*it < myNet->nodesVector.size()) {
                    f->probSubtreeInQuery += myNet->nodesVector[*it]->probInQuery;
                }
            }
           
            f->usefulAloneProb = 1.0 - f->probSubtreeInQuery;
            if (f->usefulAloneProb < 0) {
                f->usefulAloneProb = 0;
            }
            if (f->isRoot) {
                f->usefulAloneProb = 0;
            }
//            std::cout << f->factorId << " " << f->usefulAloneProb << std::endl;
        }
        
    }
    
    void FactorTree::assignUniformProbNodeInQuery(BayesianNetwork *myNet) {

//        std::cout << "here prob in query : " << prb << std::endl;
        double prb = (double) (1.0 / ((double) (myNet->nrNodes)));
        
        for (int i = 0; i < myNet->nrNodes; i++) {
            myNet->nodesVector[i]->probInQuery = prb;
        }
        
    }
    
    void FactorTree::assignRootBiasedProbNodeInQuery(BayesianNetwork *myNet) {
        
        // nodes that are closer to the root are less likely to be in query
        std::vector<double> weights(myNet->nrNodes);
        
        std::vector< std::vector<Factor*>> levels;
        getLevelOrder(rootFactor, levels);
        int nrLevels = levels.size();
        double totWeight = 0.0;
        
        for (int i = 0; i < levels.size(); i++) {
            for (int j = 0; j < levels[i].size(); j++) {
                Factor *f = levels[i][j];
                weights[f->factorId] = (double) (nrLevels - i);
                totWeight += weights[f->factorId];
            }
        }

//        std::cout << "controlling weight assignments : " << std::endl;
        for (int i = 0; i < myNet->nrNodes; i++) {
//            std::cout << i << " -- " << weights[i] << std::endl;
            double w = weights[i] / totWeight;
            myNet->nodesVector[i]->probInQuery = w;

            // cide -- here for outputting biasr query generation probs
//            std::cout << myNet->nodesVector[i]->probInQuery << std::endl;
        }
    }
    
    void FactorTree::writeProbNodeInQueryToFile(BayesianNetwork *myNet, string workloadType) {
        
        string folderName = "treeStats";
        
        string filename = folderName + OS_SEP + std::string("query_probs_") + workloadType + std::string(".txt");
        ofstream outMasterStream;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(filename);
        
        if (outMasterStream.is_open() == false) {
            cerr << "ERROR: cannot open " << filename << " for writing " << endl;
            exit(1);
        }
        else {
            for (int i = 0; i < myNet->nrNodes; i++) {
                outMasterStream << myNet->nodesVector[i]->probInQuery << std::endl;
            }
            outMasterStream.close();
        }
    
    }
    
    
} //
