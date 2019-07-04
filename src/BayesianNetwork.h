#ifndef _BAYESIAN_NET_H
#define _BAYESIAN_NET_H

#include "BayesianNode.h"


namespace qtm {
        
    class BayesianNetwork{
        
        friend std::ostream& operator<<(std::ostream& os, const Vi& vi);
        friend std::ostream& operator<<(std::ostream& os, const Vvi& vvi);
        
    public:
        
        BayesianNetwork();
        ~BayesianNetwork();
       
        int nrNodes;
        int nrEdges;
        int long long nrParameters;
        
        void setNrParameters();
        int long long getNrParameters();
        
        std::vector<BayesianNode*> nodesVector;
        void addNode(BayesianNode *myNode);
        int getNodeIdFromLabel(std::string label);
        BayesianNode* operator[](int index);
        
        bool addEdge(int fromNode, int toNode);
        bool hasEdge(int fromNode, int toNode);
        
        int returnNrEdges();
        int returnNrOutEdges(int index);
        int returnNrInEdges(int index);
        int returnNrRoots();
        
        void resetAllColors();
        void getTopologicalOrder(std::vector<int> &elimOrderVec);
        std::list<int> depthFirstSearch(int startingNode);
        std::pair<std::list<int>, int> rawDepthFirstSearch(int index);
        
        
    };
    
}

#endif // _BAYESIAN_NET_H
