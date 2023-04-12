/*
 


 */

/* 
 * File:   TreeInitializer.cpp
 * Author: virgolin
 * 
 * Created on June 27, 2018, 4:30 PM
 */

#include "GPGOMEA/Variation/TreeInitializer.h"

using namespace std;
using namespace arma;

Node * TreeInitializer::InitializeRandomTree(TreeInitShape init_shape, size_t max_height, const std::vector<Operator *> & functions, const std::vector<Operator *> & terminals, size_t nr_trees) const {

    size_t height = max_height;
    if (tree_init_type == TreeInitType::TreeInitRHH)
        height = randu() * (max_height + 1 - 2) + 2;

    Node * multitree = new Multitree();
    for (size_t i = 0; i<nr_trees;i++){
        SingleNode * n;

        if (init_shape == TreeInitFULL) {
            n = GenerateTreeFull(height, functions, terminals);
        } else if (init_shape == TreeInitGROW) {
            n = GenerateTreeGrow(height, 0, functions, terminals);
        } else {
            throw std::runtime_error("GOMEATreeInitializer::InitializeRandomTree invalid tree initialization type");
        }
        if (nr_trees == 1)
            return n;
        else multitree->AppendSingleNode(n);
    }
    return multitree;

}

SingleNode* TreeInitializer::GenerateTreeFull(size_t height_left, const std::vector<Operator*>& functions, const std::vector<Operator*>& terminals) const {
    SingleNode * curr;
    double_t rnd = arma::randu();
    size_t index;

    if (height_left > 0) {
        index = rnd * functions.size();
        curr = new SingleNode(functions[index]->Clone());
    } else {
        index = rnd * terminals.size();
        curr = new SingleNode(terminals[index]->Clone());
    }

    for (size_t i = 0; i < curr->GetArity(); i++) {
        curr->AppendChild(GenerateTreeFull(height_left - 1, functions, terminals));
    }

    return curr;
}


SingleNode* TreeInitializer::GenerateTreeGrow(size_t height_left, size_t cur_height, const std::vector<Operator*>& functions, const std::vector<Operator*>& terminals) const {

    SingleNode * curr = NULL;
    double_t rnd = arma::randu();

    size_t index;
    if (height_left > 0) {
        if (rnd >= 0.5) { // || height_left > 2 ) {
            rnd = arma::randu();
            index = rnd * functions.size();
            curr = new SingleNode(functions[index]->Clone());
        } else {
            rnd = arma::randu();
            index = rnd * terminals.size();
            curr = new SingleNode((terminals[index])->Clone());
        }

    } else {
        rnd = arma::randu();
        index = rnd * terminals.size();
        curr = new SingleNode((terminals[index])->Clone());
    }

    for (size_t i = 0; i < curr->GetArity(); i++) {
        curr->AppendChild(GenerateTreeGrow(height_left - 1, cur_height + 1, functions, terminals));
    }

    return curr;
}
