#include "BayesianNode.h"

namespace qtm {
    
    BayesianNode::BayesianNode(int varId, std::string stringLabel, int nrOfStates, std::vector<std::string> &stateLabels) {
        
        nrStates = nrOfStates,
        nodeLabel = lowercase(stringLabel);
        _id = varId;
        
        if(stateLabels.size() != nrStates) {
            std::cerr << "ERROR: there is a problem with nr of states." << std::endl;
            exit(1);
        }
        
        nodeStateLabels.reserve(nrStates);
        for (int i = 0; i < nrStates; i++) {
            nodeStateLabels.push_back(lowercase(stateLabels[i]));
        }
        
        currentColor = WHITE;
    }
    
    BayesianNode::~BayesianNode() {}
    
    std::string BayesianNode::getNodeLabel(){
        return nodeLabel;
    }
    
    int BayesianNode::getStateIdFromLabel(std::string label) {
        
        int stateId = 10000;
        trim(label);
        for (int i = 0; i < nrStates; i++) {
            if(lowercase(label).compare(nodeStateLabels[i]) == 0) {
                stateId = i;
                break;
            }
        }
        return stateId;
    }
    
    void BayesianNode::setProbabilities(std::vector<int> parent_vector, std::vector<double> prob_vector) {
        
        for (int i = 0; i < prob_vector.size(); i++) { // is aligned with node's state nr
            std::vector<int> temp_vec;
            temp_vec.push_back(i);
            
            if (!isRoot()) {
                for (int j = 0; j < parent_vector.size(); j++)
                    temp_vec.push_back(parent_vector[j]);
            }
            
            condTable.insert(std::make_pair(temp_vec, prob_vector[i]));
//            condTableStates.push_back(temp_vec);
//            condTableProbs.push_back(prob_vector[i]);
        }
        
    }
    
    bool BayesianNode::isRoot() {
        if(parents.size() == 0)
            return true;
        return false;
    }
    
    void BayesianNode::setColor(color colorToSet){
        currentColor = colorToSet;
    }
    
    BayesianNode::color BayesianNode::getColor(){
        return currentColor;
    }
    
    
    
}




