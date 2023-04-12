/*
 


 */

/* 
 * File:   Multitree.h
 * Author: Sijben
 *
 * Created on April 6
 */

#ifndef MULTITREE_H
#define MULTITREE_H

#include "GPGOMEA/Genotype/Node.h"
#include "GPGOMEA/Genotype/SingleNode.h"


class Multitree : public Node {
public:

    Multitree();

    ~Multitree() override;

    arma::mat GetOutput(const arma::mat & x, bool use_caching) override;

    void ClearCachedOutput(bool also_ancestors) override;

    Node * CloneSubtree() const override;

    std::string GetSubtreeExpression(bool only_active) override;

    std::string GetSubtreeHumanExpression() override;

    std::string GetPythonExpression() override;


    std::string GetValue() const override;

    void AppendChild(Node * c) override;

    std::vector<Node *>::iterator DetachChild(Node * c) override;

    Node * DetachChild(size_t index) override;

    void InsertChild(Node * c, std::vector<Node *>::iterator) override;

    std::vector<Node *> GetSubtreeNodes(bool only_active) override;

    Operator * ChangeOperator(const Operator & other) override;

    size_t GetHeight(bool only_active = true) const override;
    size_t GetDepth() const override;

    bool IsActive()override;

    std::vector<arma::vec> GetDesiredOutput(const std::vector<arma::vec> & Y, const arma::mat & X, bool use_caching) override;

    void ClearSubtree() override;

    Node * GetParent() override;

    std::vector<Node * > GetChildren() override;

    Operator * GetOperator() override;

    void AppendSingleNode(Node * n) override;

    int Count_N_NaComp(int count=-1)override;

    std::vector<SingleNode *> nodes;

private:

};



#endif /* MULTITREE_H */

