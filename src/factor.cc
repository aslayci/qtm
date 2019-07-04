#include "factor.h"

namespace qtm {
    
    Factor::Factor(BayesianNetwork *myNet) {
        isRoot = false;
        isLeaf = false;
        isSecondaryLeaf = false;
        factorLabel = "";
        var_to_sum_out = myNet->nodesVector.size();
        inQuery = false;
        canSkip = false;
        isIndex = false;
        isDummy = false;
        isUseful = false;
        factorTableLength = 0;
        subtreeSize = 0; // this doesn't include the "organic" nodes (i.e., factors from BN) sincde it is used for DP
        factorId = myNet->nodesVector.size();
        this->myNet = myNet;
    }
    
    Factor::Factor(BayesianNetwork *myNet, Factor *f1, Factor *f2) {
        
        this->myNet = myNet;
        
        scope = f1->scope;
        int width2 = f2->width;
        for (int i = 0; i < width2; ++i) {
            int v = f2->scope[i];
            if (!f1->inScope(v)) {
                scope.push_back(v);
            }
        }
        width = scope.size();
    
        size = 1;
        if (width > 0) {
            offset.reserve(width);
            for (int i = width-1; i >= 0; --i) {
                offset[i] = size;
                size *= myNet->nodesVector[scope[i]]->nrStates;
                var_to_index[scope[i]] = i;
            }
        }
       
    }
    
    Factor::Factor(BayesianNetwork *myNet, double value) {
        width = 0;
        size = 1;
        factorTableProbs.push_back(value);
        this->myNet = myNet;
    }
    
    void Factor::assignDomainInfo() {
        width = scope.size();
        size = 1;
        if (width > 0) {
            offset.reserve(width);
            for (int i = width-1; i >= 0; --i) {
                offset[i] = size;
                size *= myNet->nodesVector[scope[i]]->nrStates;
            }
        }
    }

    
    Factor::~Factor(){
      
    }
    
    int Factor::getElimParentOrderId(std::vector<int> &eliminationOrder, int elimIndex) { //inference
        int elimParentOrderId = eliminationOrder.size();
        for (int i = elimIndex; i < eliminationOrder.size(); i++) { //inference, starts from elimIndex
            for (int j = 0; j < scope.size(); j++) {
                if(eliminationOrder[i] == scope[j])
                       return i;
            }
            
        }
        return elimParentOrderId;
    }
    
    bool Factor::inScope(int nodeId) {
        unordered_map<int,int>::const_iterator it = var_to_index.find(nodeId);
        return (it != var_to_index.end());
    }
    
    int Factor::position_consistent_valuation(vector<int> valuation, Factor *f) {
        if (width == 0) { return 0; }
        else {
            int pos = 0;
            int index = 0;
            for (auto v : scope) {
                unordered_map<int,int>::const_iterator it_index = var_to_index.find(v);
                unordered_map<int,int>::const_iterator it_index2 = f->var_to_index.find(v);
                if (it_index != var_to_index.end() && it_index2 != var_to_index.end()) {
                    pos +=  offset[it_index->second] * valuation[it_index2->second];
                }
                index++;
            }
            return pos;
        }
    }
    
    
    Factor* Factor::product(Factor *f1, Factor *f2) {
        Factor *f3 = new Factor(myNet, f1, f2);
        vector<int> valuation(f3->width, 0);
        f3->factorTableProbs.reserve(f3->size);
//        std::cout << "f3 scope: " << f3->scope << std::endl;
        
        for (int i = 0; i < f3->size; ++i) {
            int pos1 = f1->position_consistent_valuation(valuation, f3);
            int pos2 = f2->position_consistent_valuation(valuation, f3);
            
            double value = f1->factorTableProbs[pos1] * f2->factorTableProbs[pos2];
            f3->factorTableProbs.push_back(value);
            f3->next_valuation(valuation);
        }
        
        return f3;
    }


    int Factor::position_consistent_valuation(vector<int> valuation, Factor *f, int var, int value) {
        int pos = position_consistent_valuation(valuation, f);
        unordered_map<int,int>::const_iterator it_index = var_to_index.find(var);
        if (it_index != var_to_index.end()) {
            pos += offset[it_index->second] * value;
        }
        return pos;
    }
    
    void Factor::sum_out(Factor *f) {
        
        int variable_size = myNet->nodesVector[var_to_sum_out]->nrStates;
        vector<double> values(size, 0.0);
        factorTableProbs = values;
        values.clear();
        values.shrink_to_fit();
        
        vector<int> valuation(width, 0);
        for (int i = 0; i < size; ++i) {
            for (int val = 0; val < variable_size; ++val) {
                int pos = f->position_consistent_valuation(valuation, this, var_to_sum_out, val);
                double value = f->factorTableProbs[pos];
                factorTableProbs[i] += value;
            }
            next_valuation(valuation);
        }
    }
    
    void Factor::next_valuation(vector<int> &valuation) {
        int j;
        for (j = valuation.size()-1; j >= 0 && valuation[j] == myNet->nodesVector[scope[j]]->nrStates-1; --j) {
            valuation[j] = 0;
        }
        if (j >= 0) {
            valuation[j]++;
        }
    }
    
    // executed only for secondary factors as organic factors' tables are ready from CPTs
    void Factor::createFactorTableFromSumProduct() {
        
        Factor *prod = new Factor(myNet, 1.0);
        std::vector<Factor*> prods(childFactors.size()+1);
        prods[0]=prod;
        for (int i = 1; i < childFactors.size()+1; i++) {
            prods[i] = product(prods[i-1], childFactors[i-1]);
            prods[i-1]->factorTableProbs.clear();
            prods[i-1]->factorTableProbs.shrink_to_fit();
        }
        
        if (var_to_sum_out != myNet->nodesVector.size()) {
            sum_out(prods[prods.size()-1]);
        
        }
        else {
            factorTableProbs = prods[prods.size()-1]->factorTableProbs;// cigdem: burasi this = prod olabilir mi?
            scope = prods[prods.size()-1]->scope;
            offset = prods[prods.size()-1]->offset;
            width = prods[prods.size()-1]->width;
            var_to_index = prods[prods.size()-1]->var_to_index;
            size = prods[prods.size()-1]->size;
        }
    }
    
    
    void Factor::writeIndexFactorTables(string folderName) {
        
        string filename = folderName + OS_SEP + "index_factor_" + intToStr(factorId) + ".txt";
        string delim = ",";
        ofstream outMasterStream;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(filename);
//        outMasterStream.open(filename.c_str()); c_str() is required for versions older than c++11
        
        if (outMasterStream.is_open() == false) {
            cerr << "ERROR: cannot open index file " << filename << " for writing " << endl;
            exit(1);
        }
        
        for (uint64 i = 0; i < factorTableProbs.size(); i++) {
            outMasterStream << factorTableProbs[i] << delim << " ";
        }
        
        outMasterStream.close();
    }
    
    void Factor::readFactorTablesFromIndex(string folderName) {
        
        string filename = folderName + OS_SEP + "index_factor_" + intToStr(factorId) + ".txt";
        //        std::cout << "reading " << filename << std::endl; // bura
        ifstream input_file(filename);
        
        factorTableProbs.reserve(size);
        if (input_file.is_open()) {
            for (int j = 0; j < size; ++j) {
                double value;
                read_next_double(input_file, value);
                factorTableProbs.push_back(value);
            }
            input_file.close();
        }
       
        else {
            cerr << "ERROR: couldn't read index file " << filename << endl;
            exit(1);
        }
    }
    
    // this is only for organic leaves -- don't forget to account for evidence later for Y_q case!!
    void Factor::createOrganicFactorTable() {
        // to-do: add evidence queries
        factorTableProbs.reserve(factorTableLength);
        CPT::iterator it = myNet->nodesVector[scope[0]]->condTable.begin();
        while(it != myNet->nodesVector[scope[0]]->condTable.end()) {
            //            factorTableStates.push_back(it->first);
            factorTableProbs.push_back(it->second);
            it++;
        }
        
    }
    
    
    
    
} //




