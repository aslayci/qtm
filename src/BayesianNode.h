#ifndef _BAYESIAN_NODE_H
#define _BAYESIAN_NODE_H

#include "utils.h"

namespace qtm {
    
    class BayesianNode {
        
    public:
        
        // color used for obtaining topological orders
          enum color {
           WHITE, // undiscovered
           GREY, // adjacent node that is not already discovered
           BLACK, // discovered and are adjacent to only other black or gray vertices
          };
        
        BayesianNode(int varId, std::string stringLabel, int nrOfStates, std::vector<std::string> &stateLabels);
        ~BayesianNode();
        
        int nrStates;
        int _id; 
        std::string nodeLabel;
        std::vector<std::string> nodeStateLabels;
        double probInQuery; 
        
        std::vector<int> children;
        std::vector<int> parents;

//        std::vector<std::vector<int> > condTableStates;
//        std::vector<double> condTableProbs;
        CPT condTable;
        void setProbabilities(std::vector<int>, std::vector<double>);
        
        bool isRoot();
        std::string getNodeLabel();
        void getNodeStateLabels(std::vector<std::string> &stateLabels);
        int getStateIdFromLabel(std::string label);
        
        // for dfs and topological ordering
        color currentColor;
        void setColor(color);
        color getColor();
        
    };
    
}

#endif // _BAYESIAN_NODE_H
