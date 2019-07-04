#include "BayesianNetwork.h"

namespace qtm {
    
    BayesianNetwork::BayesianNetwork() {nrNodes = 0; nrEdges = 0; nrParameters = 0;}
    
    BayesianNetwork::~BayesianNetwork(){}
    
    int BayesianNetwork::getNodeIdFromLabel(std::string label) {
        int nodeId = 10000;
        trim(label);
        for (int i = 0; i < nrNodes; i++) {
            if(lowercase(label).compare(nodesVector[i]->nodeLabel) == 0) {
                nodeId = i;
                break;
            }
        }
        return nodeId;
    }
    
    BayesianNode* BayesianNetwork::operator[](int index){
        if (index >= nodesVector.size()) throw std::domain_error("ERROR: index is out of range.");
        return nodesVector[index];
    }
    
    bool BayesianNetwork::addEdge(int fromNode, int toNode){
        if(fromNode > nodesVector.size() || toNode > nodesVector.size()){
            std::cerr << "ERROR: index is out of range." << std::endl;
            return false;
        }
        if(fromNode == toNode) return false;
        nodesVector[fromNode]->children.push_back(toNode);
        nodesVector[toNode]->parents.push_back(fromNode);
        nrEdges++;
        return true;
    }
    
    void BayesianNetwork::addNode(BayesianNode *myNode) {
        nodesVector.push_back(myNode);
        nrNodes++;
    }
    
    void BayesianNetwork::setNrParameters() {
        nrParameters = 0;
        for (int i = 0; i < nodesVector.size(); i++) {
            nrParameters += nodesVector[i]->condTable.size();
        }
    }
    
    int long long BayesianNetwork::getNrParameters() {
        return nrParameters;
    }
    
    void BayesianNetwork::getTopologicalOrder(std::vector<int> &elimOrderVec) {
        
        std::multimap<int, int> index_time_map;
        
        //iterate tru each node and apply DFS algorithm to compute finish time
        for(int nodes_counter = 0; nodes_counter < nodesVector.size(); nodes_counter++){
            auto deep_vec = depthFirstSearch(nodes_counter);
            index_time_map.insert(std::pair<int, int>(deep_vec.size(), nodes_counter));
        }
        
        for(auto it_map=index_time_map.rbegin(); it_map!=index_time_map.rend(); ++it_map){
            elimOrderVec.push_back(it_map->second);
        }
    }
    
    bool BayesianNetwork::hasEdge(int node1, int node2) {
        if (std::find(nodesVector[node1]->children.begin(), nodesVector[node1]->children.end(), node2) != nodesVector[node1]->children.end())
            return true;
        else if (std::find(nodesVector[node1]->parents.begin(), nodesVector[node1]->parents.end(), node2) != nodesVector[node1]->parents.end())
            return true;
        else
            return false;
    }
    
    void BayesianNetwork::resetAllColors(){
        for(int i = 0; i < nodesVector.size();i++){
            nodesVector[i]->setColor(BayesianNode::color::WHITE);
        }
    }
    
    std::list<int> BayesianNetwork::depthFirstSearch(int startingNode) {
        resetAllColors();
        auto deep_pair = rawDepthFirstSearch(startingNode);
        auto deep_list = deep_pair.first;
        deep_list.reverse();
        return deep_list;
    }
    
    std::pair<std::list<int>, int> BayesianNetwork::rawDepthFirstSearch(int startingNode) {
        nodesVector[startingNode]->setColor(BayesianNode::color::GREY);
        std::list<int> list_to_return;
        int multiconnections = 0;
        
        // chain iteration to all the children
        for(int i = 0; i < nodesVector[startingNode]->children.size(); i++)
        {
            int childId = nodesVector[startingNode]->children[i];
            if(nodesVector[childId]->getColor() == BayesianNode::color::BLACK){
                multiconnections += 1;
            }
            if(nodesVector[childId]->getColor() == BayesianNode::color::WHITE){
                nodesVector[childId]->setColor(BayesianNode::color::GREY);
                auto children_pair = rawDepthFirstSearch(childId);
                multiconnections += children_pair.second;
                list_to_return.merge(children_pair.first);
            }
        }
        
        // startingNode visit complete
        nodesVector[startingNode]->setColor(BayesianNode::color::BLACK);
        list_to_return.push_back(startingNode);
        return std::make_pair(list_to_return, multiconnections);
    }
    
    
    
    
} //




