#include "io_utils.h"

namespace qtm {
    
    void process_queryById_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates) {
        string delims = ",";
        vector<string> queryAtoms;
        split(queryString, delims, queryAtoms);
        int noEvidenceState = -1;
        
        for (int i = 0; i < queryAtoms.size(); i++) {
            std::size_t eqPos = queryAtoms[i].find("=");
            if (eqPos != std::string::npos) { // there is evidence assigned to this var
                int varId = strToInt(queryAtoms[i].substr(0, eqPos-1));
                int labelId = strToInt(queryAtoms[i].substr(eqPos));
                if (varId > myNet->nrNodes || labelId > myNet->nodesVector[varId]->nrStates) {
                    cerr << "ERROR: invalid query arguments, exiting " << endl;
                    exit(1);
                }
                queryVars.push_back(varId);
                evidenceStates.push_back(labelId);
                //std::cout << "found!" << '\n';
            }
            else {
                queryVars.push_back(strToInt(queryAtoms[i]));
                evidenceStates.push_back(noEvidenceState);
            }
        }
    }
    
    void process_queryByLabel_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates) {
        string delims = ",";
        vector<string> queryAtoms;
        split(queryString, delims, queryAtoms);
        int noEvidenceState = -1;
        
        for (int i = 0; i < queryAtoms.size(); i++) {
            std::size_t eqPos = queryAtoms[i].find("=");
            if (eqPos != std::string::npos) { // if there is evidence assigned to this var
                string varLabel = queryAtoms[i].substr(0, eqPos);
                string evLabel = queryAtoms[i].substr(eqPos+1);
                int varId = myNet->getNodeIdFromLabel(varLabel);
                queryVars.push_back(varId);
                int stateId = myNet->nodesVector[varId]->getStateIdFromLabel(evLabel);
                evidenceStates.push_back(stateId);
            }
            else {
                int varId = myNet->getNodeIdFromLabel(queryAtoms[i]);
                queryVars.push_back(varId);
                evidenceStates.push_back(noEvidenceState);
            }
        }
    }
    
    int process_queryFromFile_arguments(BayesianNetwork *myNet, string queryString, std::vector<int> &queryVars, std::vector<int> &evidenceStates) {
        std::size_t pos1 = queryString.find("-");
        string qFileName = queryString.substr(0, pos1);
        int qId = strToInt(queryString.substr(pos1+1));
        string queryAtomsString;
        
        ifstream input_file(qFileName);
        
        if (input_file.is_open()) {
            int counter = 0;
            while(!input_file.eof()) {
                std::string line;
                getline(input_file, line);
                if(line.empty())
                    continue;
                if (counter == qId) {
                    queryAtomsString = line;
                    break;
                }
                ++counter;
            }
            input_file.close();
        }

        else {
            cerr << "ERROR: couldn't read the input query file " << qFileName << endl;
            exit(1);
        }
        
        process_queryById_arguments(myNet, queryAtomsString, queryVars, evidenceStates);
        return qId;
        
    }
    
    bool get_variable_name(ifstream &input_file, string &token) {
        while (input_file) {
            input_file >> token;
            if (token.compare("variable") == 0) { //return true
                input_file >> token;
                token = trim(token);
                //                cout << "variable name : " << token << " ";
                return true;
            }
            if (token.compare("probability") == 0) {
                return false;
            }
        }
        return false;
    }
    
    void get_variable_states(ifstream &input_file, int &nrStates, vector<string> &stateLabels) {
        string delims = ",";
        string token, line;
        while (input_file) {
            input_file >> token;
            if (token.compare("type") == 0) {
                int cnter = 0;
                while (cnter++ < 2) { // to get rid off "discrete ["
                    input_file >> token;
                }
                input_file >> nrStates;
                //                cout << "nr states : " << nrStates << endl;
                getline(input_file, line);
                std::string::size_type pos1 = line.find_first_of("{");
                std::string::size_type pos2 = line.find_first_of("}");
                token = line.substr(pos1+1, pos2-pos1-1);
                split(token, delims, stateLabels);
                break;
            }
        }
    }
    
    void read_variables(BayesianNetwork *myNet, ifstream &input_file) {
        string varName;
        int nrStates;
        int cnter = 0;
        while(input_file) {
            if (get_variable_name(input_file, varName)) {
                //                cout << "variable read " << varName << endl;
                vector<string> stateLabels;
                get_variable_states(input_file, nrStates, stateLabels);
                //                cout << "nr states read " << nrStates << endl;
                BayesianNode *myNode = new BayesianNode(cnter, varName, nrStates, stateLabels);
                ++cnter;
                myNet->addNode(myNode);
            }
            else
                break;
        }
    }
    
    void get_bn_edges(BayesianNetwork *myNet, ifstream &input_file) {
        string line, token;
        bool initial = true;
        while(!input_file.eof()) {
            
            if (!initial) {
                input_file >> token;
            }
            
            if (token.compare("probability") == 0 || initial) {
                initial = false;
                getline(input_file, line);
                std::string::size_type pos1 = line.find_first_of("(");
                std::string::size_type pos2 = line.find_first_of(")");
                token = line.substr(pos1+1, pos2-pos1-1);
                token = trim(token);
                
                if (token.find("|") != string::npos) { // node has at least one parent
                    // get the id of the variable before |
                    token = trim(token);
                    pos1 = 0;
                    pos2 = token.find_first_of("|");
                    string childLabel = token.substr(pos1, pos2-1);
                    int childId = myNet->getNodeIdFromLabel(childLabel);
                    int parentId;
                    token = token.substr(pos2+1);
                    vector<int> parentNodeIds;
                    vector<string> parentNodeLabels;
                    int nrNextLinesToRead = 1;
                    
                    if (token.find(",") != string::npos) { // if exists more than one parent for this node
                        split(token, ",", parentNodeLabels);
                        for (int i = 0 ; i < parentNodeLabels.size(); i++) {
                            // add the edges from each parent
                            parentId = myNet->getNodeIdFromLabel(parentNodeLabels[i]);
                            parentNodeIds.push_back(parentId);
                            myNet->addEdge(parentId, childId);
                            nrNextLinesToRead *= myNet->nodesVector[parentId]->nrStates;
                        }
                    }
                    else {
                        parentNodeLabels.push_back(trim(token));
                        parentId = myNet->getNodeIdFromLabel(parentNodeLabels[0]);
                        parentNodeIds.push_back(parentId);
                        myNet->addEdge(parentId, childId);
                        nrNextLinesToRead = myNet->nodesVector[parentId]->nrStates;
                    }
                    
                    int sayac = 0;
                    while (sayac++ < nrNextLinesToRead) { // read cond. prob. table line by line from the next nrNextLinesToRead lines
                        vector<int> parentStateIds;
                        vector<string> parentStateLabels;
                        
                        getline(input_file, line);
                        pos1 = line.find_first_of("(");
                        pos2 = line.find_first_of(")");
                        token = line.substr(pos1+1, pos2-pos1-1);
                        
                        if (token.find(",") != string::npos) {
                            split(token, ",", parentStateLabels);
                        }
                        else {
                            parentStateLabels.push_back(trim(token));
                        }
                        for (int i = 0 ; i < parentStateLabels.size(); i++) { // get state id of parents from the label of their states
                            parentId = parentNodeIds[i];
                            parentStateIds.push_back(myNet->nodesVector[parentId]->getStateIdFromLabel(parentStateLabels[i]));
                        }
                        token = line.substr(pos2-pos1+3);
                        vector<double> condTableLine;
                        splitDouble(token,",", condTableLine);
                        myNet->nodesVector[childId]->setProbabilities(parentStateIds,condTableLine);
                    }
                } // end if node has at least one parent
                
                else { // node is a root of BN
                    int nid = myNet->getNodeIdFromLabel(token);
                    input_file >> token; // getting rid of the word "table"
                    getline(input_file, token);
                    vector<double> condTableLine;
                    splitDouble(token,",", condTableLine);
                    myNet->nodesVector[nid]->setProbabilities({},condTableLine);
                }
                continue;
            }
        }
    }
    
    void read_bif_format(BayesianNetwork *myNet, string filename) {
        
        ifstream input_file(filename);
        
        if (input_file.is_open()) {
            read_variables(myNet, input_file);
            get_bn_edges(myNet, input_file);
            input_file.close();
        }
        
        else {
            cerr << "ERROR: couldn't read the file " << filename << endl;
            exit(1);
        }
        
    }
    
    // this is to read the nr of states before creating the BN object
    //    void get_nr_states_for_inits(string filename, vector<int int> &varStatesVec) {
    //
    //        ifstream input_file(filename);
    //
    //        if (input_file.is_open()) {
    //            string delims = ",";
    //            string token, line;
    //            int counter = 0, hede = -1;
    //
    //            while (input_file) {
    //                input_file >> token;
    //                if (token.compare("type") == 0) {
    //                    int cnter = 0;
    //                    while (cnter++ < 2) { // to get rid off "discrete ["
    //                        input_file >> token;
    //                    }
    //                    input_file >> hede;
    //                    counter++;
    //                    varStatesVec.push_back(hede);
    ////                    if(counter == nrVariables)
    ////                        break;
    //                }
    //
    //                if (token.compare("probability") == 0)
    //                    break;
    //            }
    //            input_file.close();
    //        }
    //
    //        else {
    //            cerr << "ERROR: couldn't read the file" << filename << endl;
    //            exit(1);
    //        }
    //    }
    
    
    
}

