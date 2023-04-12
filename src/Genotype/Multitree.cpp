/*



 */

/*
 * File:   Multitree.cpp
 * Author: Sijben
 *
 * Created on April, 6
 */

#include "GPGOMEA/Genotype/Multitree.h"

using namespace std;
using namespace arma;

Multitree::Multitree() : Node(NodeType::Multi)  {
}

Multitree::~Multitree() {
}

arma::mat Multitree::GetOutput(const arma::mat & x, bool use_caching) {
    if (use_caching && !cached_output.empty())
        return cached_output;
    size_t output_size = nodes[0]->GetOutput(x,false).n_elem;
    size_t nr_trees = nodes.size();
    arma::mat out =  arma::mat(output_size, nr_trees, arma::fill::none);
    for (size_t i = 0; i<nodes.size(); i++){
        arma::vec output_per_tree = nodes[i]->GetOutput(x , false);
        out.col(i) = output_per_tree;
    }
    if (use_caching)
        cached_output = out;

    return out;
}


void Multitree::ClearCachedOutput(bool also_ancestors) {
    cached_output.clear();
}


Node* Multitree::CloneSubtree() const {
    Node * new_node = new Multitree();
    for (SingleNode * n: nodes)
        new_node->AppendSingleNode(n->CloneSubtree());
    new_node->cached_objectives = this->cached_objectives;
    new_node->cached_objectives_test = this->cached_objectives_test;
    return new_node;
}

std::string Multitree::GetSubtreeExpression(bool only_active_nodes) {
    string expr;
    size_t index = 1;
    for (SingleNode * n : nodes){
        expr += "Expression" + std::to_string(index);
        expr += n->GetSubtreeExpression(only_active_nodes);
        expr += "\n ";
        index++;
    }
    return expr;

}

std::string Multitree::GetSubtreeHumanExpression() {
    string expr;
    size_t index = 1;
    for (SingleNode * n : nodes){
        expr += "Expression " + std::to_string(index) + ": ";
        expr += n->GetSubtreeHumanExpression();
//        expr += "\n ";
        index++;

    }
    return expr;
}

std::string Multitree::GetPythonExpression() {
    string expr;
    size_t index = 1;
    for (SingleNode * n : nodes){
        expr += "Expression " + std::to_string(index) + ": ";
        expr += n->GetPythonExpression();
//        expr += "\n ";
        index++;

    }
    return expr;
}

std::string Multitree::GetValue() const{
    throw NotImplementedException("Multitree::GetValue not implemented.");
}

void Multitree::AppendChild(Node * c){
    throw NotImplementedException("Multitree::AppendChild not implemented.");
}

std::vector<Node *>::iterator Multitree::DetachChild(Node * c){
    throw NotImplementedException("Multitree::DetachChild not implemented.");
}

Node *  Multitree::DetachChild(size_t index){
    throw NotImplementedException("Multitree::DetachChild not implemented.");
}

void  Multitree::InsertChild(Node * c, std::vector<Node *>::iterator){
    throw NotImplementedException("Multitree::InsertChild not implemented.");
}

std::vector<Node*> Multitree::GetSubtreeNodes(bool only_active) {
    std::vector<Node *> subtree_nodes;
    for (SingleNode* n : nodes){
        std::vector<Node *> subtree_nodes_per_tree = n->GetSubtreeNodes(only_active);
        subtree_nodes.insert(subtree_nodes.end(), subtree_nodes_per_tree.begin(), subtree_nodes_per_tree.end());
    }
    return subtree_nodes;
}

Operator *  Multitree::ChangeOperator(const Operator & other){
    throw NotImplementedException("Multitree::ChangeOperator not implemented.");
}


size_t Multitree::GetHeight(bool only_active) const{
    throw NotImplementedException("Multitree::GetHeight not implemented.");
}

size_t Multitree::GetDepth() const {
    throw NotImplementedException("Multitree::GetDepth not implemented.");
}

bool Multitree::IsActive(){
    throw NotImplementedException("Multitree::IsActive not implemented.");
}

std::vector<arma::vec> Multitree::GetDesiredOutput(const std::vector<arma::vec> & Y, const arma::mat & X, bool use_caching){
    throw NotImplementedException("Multitree::GetDesiredOutput not implemented.");
}
void Multitree::ClearSubtree() {
    for (SingleNode * n : nodes){
        n->ClearSubtree();
    }
    delete this;
}

Node * Multitree::GetParent(){
    throw NotImplementedException("Multitree::GetParent not implemented.");
}

std::vector<Node * > Multitree::GetChildren(){
    throw NotImplementedException("Multitree::GetChildren not implemented.");
}


Operator * Multitree::GetOperator(){
    throw NotImplementedException("Multitree::GetOperator not implemented.");
}
void Multitree::AppendSingleNode(Node *  n){
    this->nodes.push_back((SingleNode * ) n);
}

int Multitree::Count_N_NaComp(int count) {
    int c;
    for (SingleNode * n: nodes){
        c += n->Count_N_NaComp(count);
    }
    return c;
}