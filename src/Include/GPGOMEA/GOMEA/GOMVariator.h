/*
 


 */

/* 
 * File:   GOMVariator.h
 * Author: virgolin
 *
 * Created on June 28, 2018, 6:07 PM
 */

#ifndef GOMVARIATOR_H
#define GOMVARIATOR_H

#include "GPGOMEA/Genotype/Node.h"
#include "GPGOMEA/Operators/Operator.h"
#include "GPGOMEA/Fitness/Fitness.h"
#include "GPGOMEA/GOMEA/GOMEATreeInitializer.h"
#include "GPGOMEA/Utils/Exceptions.h"
#include "GPGOMEA/Semantics/SemanticLibrary.h"
#include "GPGOMEA/Semantics/SemanticBackpropagator.h"
#include "GPGOMEA/Variation/SubtreeVariator.h"

#include <vector>
#include <armadillo>

class MOArchive;

class GOMVariator {
public:

    GOMVariator() {
    };
    static Node * GOMMO(const Node &sol,const std::vector<Operator *> &functions, const std::vector<Operator *> &terminals, const std::vector<std::vector<size_t> > &FOS,
                        Fitness &fit, bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem,
                        size_t NIS_const);

    static Node *
    GOMMO(Node *&sol, const std::vector<Node *> &donors, const std::vector<std::vector<size_t>> &FOS, Fitness &fit,
          bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem, size_t NIS_const);

    static Node *GOMMOFI(Node *&sol, const std::vector<std::vector<size_t> > &FOS,
                         Fitness &fit, bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem,
                         size_t NIS_const);

    static Node *
    GOM(const Node &sol, const std::vector<Node *> &donors, const std::vector<std::vector<size_t>> &FOS, Fitness &fit,
        bool use_caching);

    static Node *
    GOM(const Node &sol, const std::vector<Operator *> &functions, const std::vector<Operator *> &terminals,
        const std::vector<std::vector<size_t>> &FOS, Fitness &fit, bool use_caching);

    static std::pair<bool, bool>
    check_changes_SO(Node *offspring, arma::vec back_obj, bool FI, size_t obj);

    static std::pair<bool, bool>
    check_changes_MO(Node *offspring, arma::vec back_obj, bool FI, MOArchive *mo_archive);

    static Node *MakeBiggerGOMEATree(const Node &original, size_t orig_height, size_t desired_height, size_t max_arity,
                                     const GOMEATreeInitializer &tree_init, const std::vector<Operator *> &functions,
                                     const std::vector<Operator *> &terminals);

private:


};

#endif /* GOMVARIATOR_H */

