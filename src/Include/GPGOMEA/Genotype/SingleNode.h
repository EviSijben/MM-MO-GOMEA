/*
 


 */

/* 
 * File:   Node.h
 * Author: virgolin
 *
 * Created on June 27, 2018, 11:57 AM
 * 
 * ##########################################
 * 
 * New filename: SingleNode.h
 * Adjusted by: Sijben 
 * 
 * Adjusted on: April 7
 */

#ifndef SINGLENODE_H
#define SINGLENODE_H

#include "GPGOMEA/Genotype/Node.h"

class SingleNode : public Node {
public:

    SingleNode(Operator * op);
    SingleNode(const SingleNode& orig);
    ~SingleNode() override;

    arma::mat GetOutput(const arma::mat & x, bool use_caching) override;

    void ClearCachedOutput(bool also_ancestors) override;

    Node * CloneSubtree() const override;

    std::string GetSubtreeExpression(bool only_active) override;

    std::string GetSubtreeHumanExpression() override;

    std::string GetPythonExpression() override;

    void  GetPythonExpR(std::string& expr);

    size_t GetArity() const override;

    std::string GetValue() const override;

    void AppendChild(Node * c) override;

    std::vector<Node *>::iterator DetachChild(Node * c) override;

    Node * DetachChild(size_t index) override;

    void InsertChild(Node * c, std::vector<Node *>::iterator) override;

    std::vector<Node *> GetSubtreeNodes(bool only_active) override;

    Operator * ChangeOperator(const Operator & other) override;

    size_t GetHeight(bool only_active) const override; // default true
    size_t GetDepth() const override;

    size_t GetRelativeIndex() override;

    bool IsActive()override;

    std::vector<arma::vec> GetDesiredOutput(const std::vector<arma::vec> & Y, const arma::mat & X, bool use_caching) override;

    void ClearSubtree() override;

    Node * GetParent() override;

    std::vector<Node *> GetChildren() override;

    Operator * GetOperator() override;

    void AppendSingleNode(Node * n) override;

    int Count_N_NaComp(int count)override; // default -1

private:
    Node * parent;
    std::vector<Node *> children;

    Operator * op;

    void ComputeHeightRecursive(size_t & h, bool only_active) const;

    void GetSubtreeNodesPreorderRecursive(std::vector<Node*> & v, bool only_active);
    void GetSubtreeNodesPostorderRecursive(std::vector<Node*> & v, bool only_active);
    void GetSubtreeNodesLevelOrder(std::vector<Node*> & v,bool only_active);
    void GetSubtreeExpressionRecursive(std::string & expr, bool only_active);
    void GetSubtreeHumanExpressionRecursive(std::string & expr);
};



#endif /* SINGLENODE_H */

