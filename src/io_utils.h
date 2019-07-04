#ifndef IO_UTILS_H
#define IO_UTILS_H

#include "utils.h"
#include "BayesianNetwork.h"

namespace qtm {
    
    // bn file operations
    bool get_variable_name(ifstream &input_file, string &token);
    void get_variable_states(ifstream &input_file, int &nrStates);
    void read_variables(BayesianNetwork *myNet, ifstream &input_file);
    void get_bn_edges(BayesianNetwork *myNet, ifstream &input_file);
    void read_bif_format(BayesianNetwork *myNet, string filename);
//    void get_nr_states_for_inits(string filename, vector<int> &varStatesVec);
    
    // query argument operations
    void process_queryById_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates);
    void process_queryByLabel_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates);
    int process_queryFromFile_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates);
    
    
}



#endif
