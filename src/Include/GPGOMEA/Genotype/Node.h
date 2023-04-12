/*
 


 */

/* 
 * File:   Node.h
 * Author: virgolin
 *
 * Created on June 27, 2018, 11:57 AM
 */

#ifndef NODE_H
#define NODE_H

#include "GPGOMEA/Utils/Exceptions.h"
#include "GPGOMEA/Operators/Operator.h"

#include <vector>
#include <armadillo>
#include <iostream>
#include <assert.h>
#include <queue>

enum SubtreeParseType {
    SubtreeParsePREORDER, SubtreeParsePOSTORDER, SubtreeParseLEVELORDER
};

enum NodeType {
    Single, Multi
};


class Node {
public:

    Node(NodeType type);

    virtual ~Node();

    virtual arma::mat GetOutput(const arma::mat & x, bool use_caching);

    virtual void ClearCachedOutput(bool also_ancestors);

    virtual Node * CloneSubtree() const;

    virtual std::string GetSubtreeExpression(bool only_active);

    virtual std::string GetSubtreeHumanExpression();

    virtual std::string GetPythonExpression();

    virtual size_t GetArity() const;

    virtual std::string GetValue() const;

    virtual void AppendChild(Node * c);

    virtual std::vector<Node *>::iterator DetachChild(Node * c);

    virtual Node * DetachChild(size_t index);

    virtual void InsertChild(Node * c, std::vector<Node *>::iterator);

    virtual std::vector<Node *> GetSubtreeNodes(bool only_active);

    virtual Operator * ChangeOperator(const Operator & other);

    virtual size_t GetHeight(bool only_active = true) const;
    virtual size_t GetDepth() const;

    virtual size_t GetRelativeIndex();

    virtual bool IsActive();

    virtual std::vector<arma::vec> GetDesiredOutput(const std::vector<arma::vec> & Y, const arma::mat & X, bool use_caching);

    virtual void ClearSubtree();

    virtual Node* GetParent();

    virtual std::vector<Node * > GetChildren();

    virtual Operator * GetOperator();

    virtual void AppendSingleNode(Node * n);

    bool Dominates(Node * o);

    virtual int Count_N_NaComp(int count=-1);

    double_t cached_fitness = arma::datum::inf;


    arma::mat cached_output;
    NodeType type;


    arma::vec cached_objectives;
    arma::vec cached_objectives_test;


    double_t crowding_distance;
    size_t rank;
    size_t NIS = 0;


private:
    void ComputeHeightRecursive(size_t & h, bool only_active) const;

    void GetSubtreeNodesPreorderRecursive(std::vector<Node*> & v, bool only_active);
    void GetSubtreeNodesPostorderRecursive(std::vector<Node*> & v, bool only_active);
    void GetSubtreeNodesLevelOrder(std::vector<Node*> & v,bool only_active);
    void GetSubtreeExpressionRecursive(std::string & expr, bool only_active);
    void GetSubtreeHumanExpressionRecursive(std::string & expr);
};



#endif /* NODE_H */

