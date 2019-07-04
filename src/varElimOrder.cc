#include "varElimOrder.h"


namespace qtm {
    
    void computeEliminationOrder(BayesianNetwork *myNet, std::vector<int> &elimOrderVec, string elimOrderOpt) {
        
        if (elimOrderOpt.compare("ts") == 0) {
            myNet->getTopologicalOrder(elimOrderVec);
        }
        
        else if (elimOrderOpt.compare("rnd") == 0) {
            for (int i = 0; i < myNet->nrNodes; i++) {
                elimOrderVec.push_back(i);
            }
            std::random_shuffle(elimOrderVec.begin(), elimOrderVec.end());
        }
        
        else if (elimOrderOpt.compare("my") == 0)
            my_min_fill(myNet,elimOrderVec);
        
        else {
            std::unordered_set<int> vars;
            vars.reserve(myNet->nrNodes);
            std::vector<std::vector<int> > neighbors(myNet->nrNodes);
            int maxNrStates = 0;
            int best_var;
            
            for (int i = 0; i < myNet->nrNodes; i++) {
                vars.insert(i);
                neighbors[i].reserve(myNet->nodesVector[i]->children.size() + myNet->nodesVector[i]->parents.size());
                for (int j = 0; j < myNet->nodesVector[i]->children.size(); j++) {
                    neighbors[i].push_back(myNet->nodesVector[i]->children[j]);
                }
                for (int j = 0; j < myNet->nodesVector[i]->parents.size(); j++) {
                    neighbors[i].push_back(myNet->nodesVector[i]->parents[j]);
                }
                // merry the parents
                for (int j = 0; j < myNet->nodesVector[i]->parents.size(); j++) {
                    for (int z = 0; z < myNet->nodesVector[i]->parents.size(); z++) {
                        int parent1 = myNet->nodesVector[i]->parents[j];
                        int parent2 = myNet->nodesVector[i]->parents[z];
                        if (j < z && !(myNet->hasEdge(parent1, parent2))) {
                            neighbors[parent1].push_back(parent2);
                            neighbors[parent2].push_back(parent1);
                        }
                    }
                }
                
                if (myNet->nodesVector[i]->nrStates > maxNrStates) {
                    maxNrStates = myNet->nodesVector[i]->nrStates;
                }
            }
            
            while (!vars.empty()) {
                //# choose an elimination ordering: ts, rnd, wmf, mf, mn, mw
                if (elimOrderOpt.compare("mn") == 0)
                    best_var = minimum_neighbors(neighbors,vars);
                
                else if (elimOrderOpt.compare("mw") == 0)
                    best_var = minimum_weight(myNet,neighbors,vars,maxNrStates);
                
                else if (elimOrderOpt.compare("mf") == 0)
                    best_var = min_fill(neighbors,vars);
                
                else if (elimOrderOpt.compare("wmf") == 0)
                    best_var = weighted_min_fill(myNet,neighbors,vars,maxNrStates);
                
                
                else {
                    std::cerr << "ERROR: unrecognized elimination ordering option, exiting." << std::endl;
                    exit(1);
                }
                
                elimOrderVec.push_back(best_var);
                vars.erase(best_var);
                
                // erase best_var from graph
                for (int i = 0; i < neighbors[best_var].size(); i++) {
                    int id = neighbors[best_var][i];
                    auto it = std::find(neighbors[id].begin(), neighbors[id].end(), best_var);
                    if(it != neighbors[id].end())
                        neighbors[id].erase(it);
                }
                
                // add fill-in edges
                for (int i = 0; i < neighbors[best_var].size(); i++) {
                    for (int j = 0; j < neighbors[best_var].size(); j++) {
                        if (i < j && !isConnected(neighbors,neighbors[best_var][i],neighbors[best_var][j])) {
                            neighbors[neighbors[best_var][i]].push_back(neighbors[best_var][j]);
                            neighbors[neighbors[best_var][j]].push_back(neighbors[best_var][i]);
                        }
                    }
                }
            }
        }
    }
    
    int minimum_neighbors(std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars) {
        
        int best_var = *(vars.begin());
        int minCost = neighbors.size() + 1;
        std::unordered_set<int>::iterator it;
        
        for (it = vars.begin(); it != vars.end(); it++) {
            int cost = neighbors[*it].size();
            if (cost < minCost) {
                best_var = *it;
                minCost = cost;
            }
        }
        return best_var;
    }
    
    int minimum_weight(BayesianNetwork *myNet, std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars, int maxNrStates) {
        
        int best_var = *(vars.begin());
        int minCost = neighbors.size() * maxNrStates + 1;
        std::unordered_set<int>::iterator it;
        
        for (it = vars.begin(); it != vars.end(); it++) {
            int cost = 1;
            for (int i = 0; i < neighbors[*it].size(); i++) {
                cost *= myNet->nodesVector[neighbors[*it][i]]->nrStates;
            }
            if (cost < minCost) {
                best_var = *it;
                minCost = cost;
            }
        }
        return best_var;
    }
    
    // bunu bozmadim di mi la
    int min_fill(std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars) {
        
        int best_var = *(vars.begin());
        int minCost = neighbors.size() + 1;
        std::unordered_set<int>::iterator it;

        for (it = vars.begin(); it != vars.end(); it++) {
            int cost = 0;
            
            for (int i = 0; i < neighbors[*it].size(); i++) {
                for (int j = 0; j < neighbors[*it].size(); j++) {
                    if (i < j && !isConnected(neighbors,neighbors[*it][i],neighbors[*it][j])) {
                        ++cost;
                    }
                }
            }
            
            if (cost < minCost) {
                best_var = *it;
                minCost = cost;
            }
            
            else if (cost == minCost) { // use minimum neighbor as tie breaker
                
                if (neighbors[*it].size() < neighbors[best_var].size()) {
                    best_var = *it;
                }
            }
        }
        return best_var;
    }
    
    
    int weighted_min_fill(BayesianNetwork *myNet, std::vector<std::vector<int> > &neighbors, std::unordered_set<int> &vars, int maxNrStates) {
        
        int best_var = *(vars.begin());
        int minCost = neighbors.size() * maxNrStates + 1;
        std::unordered_set<int>::iterator it;
        
        for (it = vars.begin(); it != vars.end(); it++) {
            int cost = 0;
            for (int i = 0; i < neighbors[*it].size(); i++) {
                for (int j = 0; j < neighbors[*it].size(); j++) {
                    if (i < j && !isConnected(neighbors,neighbors[*it][i],neighbors[*it][j])) {
                        cost += (myNet->nodesVector[neighbors[*it][i]]->nrStates * myNet->nodesVector[neighbors[*it][j]]->nrStates);
                    }
                }
            }
            
            if (cost < minCost) {
                best_var = *it;
                minCost = cost;
            }
            
            else if (cost == minCost) { // use minimum neighbor as tie breaker
                
                if (neighbors[*it].size() < neighbors[best_var].size()) {
                    best_var = *it;
                }
            }
        }
        return best_var;
    }
    
    void my_min_fill(BayesianNetwork *myNet, std::vector<int> &elimOrderVec) {
        
        std::unordered_set<int> vars;
        vars.reserve(myNet->nrNodes);
        std::vector<std::vector<int> > neighbors(myNet->nrNodes);
        int maxNrStates = 0;
        int best_var;
        
        for (int i = 0; i < myNet->nrNodes; i++) {
            vars.insert(i);
            neighbors[i].reserve(myNet->nodesVector[i]->children.size() + myNet->nodesVector[i]->parents.size());
            for (int j = 0; j < myNet->nodesVector[i]->children.size(); j++) {
                neighbors[i].push_back(myNet->nodesVector[i]->children[j]);
            }
            for (int j = 0; j < myNet->nodesVector[i]->parents.size(); j++) {
                neighbors[i].push_back(myNet->nodesVector[i]->parents[j]);
            }
            
            // merry the parents
            for (int j = 0; j < myNet->nodesVector[i]->parents.size(); j++) {
                for (int z = 0; z < myNet->nodesVector[i]->parents.size(); z++) {
                    int parent1 = myNet->nodesVector[i]->parents[j];
                    int parent2 = myNet->nodesVector[i]->parents[z];
                    if (j != z && !(myNet->hasEdge(parent1, parent2))) {
                        neighbors[parent1].push_back(parent2);
                        neighbors[parent2].push_back(parent1);
                    }
                }
            }
            
            if (myNet->nodesVector[i]->nrStates > maxNrStates) {
                maxNrStates = myNet->nodesVector[i]->nrStates;
            }
        }
        
        std::vector<std::vector<int> > scopes(myNet->nrNodes,std::vector<int>());
        std::vector<int> elimParents(myNet->nrNodes);
        std::vector<std::vector<int> > tree_rr(myNet->nrNodes,std::vector<int>());
        
        while (!vars.empty()) {
            best_var = *(vars.begin());
            int minCost = neighbors.size() * maxNrStates + 1;
            std::unordered_set<int>::iterator it;
            std::vector<uint64> costs;
            std::unordered_set<int> unrolled_scope;
            
            for (it = vars.begin(); it != vars.end(); it++) {
                int cost = 1;
                for (int i = 0; i < neighbors[*it].size(); i++) {
                    unrolled_scope.insert(neighbors[*it][i]);
                }
                unrolled_scope.insert(*it);
//                cost *= myNet->nodesVector[neighbors[*it][i]]->nrStates;
//                cost *= myNet->nodesVector[*it]->nrStates;
                
//                std::unordered_set<int> unrolled_scope;
//                for (int i = 0; i < neighbors[*it].size(); i++) {
//                    for (int j = 0; j < neighbors[*it].size(); j++) {
//                        if (i < j && !isConnected(neighbors,neighbors[*it][i],neighbors[*it][j])) {
//                            unrolled_scope.insert(neighbors[*it][i]);
//                            unrolled_scope.insert(neighbors[*it][j]);
//                        }
//                    }
//                }
//                unrolled_scope.insert(*it);
//                std::unordered_set<int>::iterator it2;
//                for (it2 = unrolled_scope.begin(); it2 != unrolled_scope.end(); it2++) {
//                    cost *= myNet->nodesVector[*it2]->nrStates;
//                }
                ///
                
                for (int i = 0; i < tree_rr[*it].size(); i++) {
                    bool elimChildrenCand = true;
                    for (int j = i; j < tree_rr[*it].size(); j++) {
                        if(elimParents[tree_rr[*it][i]] == tree_rr[*it][j]) {
                            elimChildrenCand = false;
                            break;
                        }
                    }
                    
                    if (elimChildrenCand) {
                        for (int j = 0; j < scopes[elimChildrenCand].size(); j++) {
                            unrolled_scope.insert(scopes[elimChildrenCand][j]);
                        }
                    }
                }
      
                std::unordered_set<int>::iterator it2;
                for (it2 = unrolled_scope.begin(); it2 != unrolled_scope.end(); it2++) {
                    cost *= myNet->nodesVector[*it2]->nrStates;
                }
                
                ///
                if (cost < minCost) {
                    best_var = *it;
                    minCost = cost;
                }
                
                else if (cost == minCost) { // use minimum neighbor as tie breaker
                    if (neighbors[*it].size() < neighbors[best_var].size()) {
                        best_var = *it;
                    }
                }
            } // select best_var
            
            elimOrderVec.push_back(best_var);
            for (int i = 0; i < tree_rr[best_var].size(); i++) { // assign this as elim parent to children
                bool elimChildrenCand = true;
                for (int j = i; j < tree_rr[best_var].size(); j++) {
                    if(elimParents[tree_rr[best_var][i]] == tree_rr[best_var][j]) {
                        elimChildrenCand = false;
                        break;
                    }
                }
                
                if (elimChildrenCand) {
                    elimParents[tree_rr[best_var][i]] = best_var;
                }
            }
            
            for (int i = 0; i < neighbors[best_var].size(); i++) {
                tree_rr[neighbors[best_var][i]].push_back(best_var);
                scopes[best_var].push_back(neighbors[best_var][i]);
            }
                
            
//            std::unordered_set<int> scope;
//            for (int i = 0; i < neighbors[best_var].size(); i++) {
//                for (int j = 0; j < neighbors[best_var].size(); j++) {
//                    if (i < j && !isConnected(neighbors,neighbors[best_var][i],neighbors[best_var][j])) {
//                        scope.insert(neighbors[best_var][i]);
//                    }
//                }
//            }
            vars.erase(best_var);
            
            for (int i = 0; i < neighbors[best_var].size(); i++) {
                int id = neighbors[best_var][i];
                auto it = std::find(neighbors[id].begin(), neighbors[id].end(), best_var);
                if(it != neighbors[id].end())
                    neighbors[id].erase(it);
            }
            
            // add fill-in edges
            for (int i = 0; i < neighbors[best_var].size(); i++) {
                for (int j = 0; j < neighbors[best_var].size(); j++) {
                    if (i != j && !isConnected(neighbors,neighbors[best_var][i],neighbors[best_var][j])) {
                        neighbors[neighbors[best_var][i]].push_back(neighbors[best_var][j]);
                        neighbors[neighbors[best_var][j]].push_back(neighbors[best_var][i]);
                    }
                }
            }
            
        }
        

        
    }
    
    bool isConnected(std::vector<std::vector<int> > &neighbors, int node1, int node2) {
        
        if (neighbors[node1].size() < neighbors[node2].size()) {
            if (std::find(neighbors[node1].begin(), neighbors[node1].end(), node2) != neighbors[node1].end())
                return true;
            else
                return false;
        }
        
        else {
            if (std::find(neighbors[node2].begin(), neighbors[node2].end(), node1) != neighbors[node2].end())
                return true;
            else
                return false;
        }
    }
    
    
    
    
    
} //




