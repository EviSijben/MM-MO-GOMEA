/*



 */

/*
 * File:   Node.cpp
 * Author: Sijben
 *
 * Created on April 7, 2021
 */

#include "GPGOMEA/Genotype/Node.h"

using namespace std;
using namespace arma;

Node::Node(NodeType type) : type(type) {
}

Node::~Node() {
}

arma::mat Node::GetOutput(const arma::mat & x, bool use_caching){
    throw NotImplementedException("Node::GetOutput not implemented.");
}

void Node::ClearCachedOutput(bool also_ancestors){
    throw NotImplementedException("Node::ClearCachedOutput not implemented.");
}

Node* Node::CloneSubtree() const{
    throw NotImplementedException("Node::CloneSubtree not implemented.");
}

std::string Node::GetSubtreeExpression(bool only_active){
    throw NotImplementedException("Node::GetSubtreeExpression not implemented.");
}

std::string Node::GetSubtreeHumanExpression(){
    throw NotImplementedException("Node::GetSubtreeHumanExpression not implemented.");
}

std::string Node::GetPythonExpression(){
    throw NotImplementedException("Node::GetDuckDBReadableExpression not implemented.");
}

std::vector<Node *> Node::GetSubtreeNodes(bool only_active){
    throw NotImplementedException("Node::GetSubtreeNodes not implemented.");
}

void Node::ClearSubtree(){
    throw NotImplementedException("Node::ClearSubtree not implemented.");
}
size_t Node::GetArity() const {
    throw NotImplementedException("Node::GetArity not implemented.");
}

std::string Node::GetValue() const{
    throw NotImplementedException("Node::GetValue not implemented.");
}

void Node::AppendChild(Node * c){
    throw NotImplementedException("Node::AppendChild not implemented.");
}

std::vector<Node *>::iterator Node::DetachChild(Node * c){
    throw NotImplementedException("Node::DetachChild not implemented.");
}

Node *  Node::DetachChild(size_t index){
    throw NotImplementedException("Node::DetachChild not implemented.");
}

void  Node::InsertChild(Node * c, std::vector<Node *>::iterator){
    throw NotImplementedException("Node::InsertChild not implemented.");
}

Operator *  Node::ChangeOperator(const Operator & other){
    throw NotImplementedException("Node::ChangeOperator not implemented.");
}


size_t Node::GetHeight(bool only_active) const{
    throw NotImplementedException("Node::GetHeight not implemented.");
}

size_t Node::GetDepth() const {
    throw NotImplementedException("Node::GetDepth not implemented.");
}
size_t Node::GetRelativeIndex() {
    throw NotImplementedException("Node::GetRelativeIndex not implemented.");

}
bool Node::IsActive(){
    throw NotImplementedException("Node::IsActive not implemented.");
}

std::vector<arma::vec> Node::GetDesiredOutput(const std::vector<arma::vec> & Y, const arma::mat & X, bool use_caching){
    throw NotImplementedException("Node::GetDesiredOutput not implemented.");
}

Node * Node::GetParent(){
    throw NotImplementedException("Node::GetParent not implemented.");
}

std::vector<Node * > Node::GetChildren(){
    throw NotImplementedException("Node::GetChildren not implemented.");
}

Operator * Node::GetOperator(){
    throw NotImplementedException("Node::GetOperator not implemented.");
}

void Node::AppendSingleNode(Node *  n){
    throw NotImplementedException("Node::AppendSingleNode not implemented.");
}

bool Node::Dominates(Node * o) {
    bool strictly_better_somewhere = false;
    for(size_t i = 0; i < this->cached_objectives.n_elem; i++) {
        if (this->cached_objectives[i] < o->cached_objectives[i])
            strictly_better_somewhere = true;
        else if (this->cached_objectives[i] > o->cached_objectives[i])
            return false;
    }
    return strictly_better_somewhere;
}
int Node::Count_N_NaComp(int count) {
    throw NotImplementedException("Node::Count_N_NaComp not implemented.");
}
