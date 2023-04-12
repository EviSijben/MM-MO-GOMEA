//
// Created by evi on 4/26/21.
//
#ifndef GOMEAMOGENERATIONHANDLER_H
#define GOMEAMOGENERATIONHANDLER_H

#include "GPGOMEA/Evolution/GenerationHandler.h"
#include "GPGOMEA/GOMEA/GOMEAFOS.h"
#include "GPGOMEA/GOMEA/GOMVariator.h"
#include "GPGOMEA/Evolution/PopulationInitializer.h"
#include "GPGOMEA/Evolution/MOArchive.h"

#include <armadillo>

class GOMEAMOGenerationHandler : public GenerationHandler {
public:

    GOMEAMOGenerationHandler(ConfigurationOptions * conf, TreeInitializer * tree_initializer, Fitness * fitness, SemanticLibrary * semlib = NULL, SemanticBackpropagator * semback = NULL)
            : GenerationHandler(conf, tree_initializer, fitness, semlib, semback) {
       linkage_normalization_matrix = new arma::mat();
    };

    virtual ~GOMEAMOGenerationHandler() {
//        for (arma::mat * lnm : linkage_normalization_matrix)
//            if (lnm)
//                delete lnm;
    }

    void PerformGeneration(std::vector<Node *> & population) override;


    void InitLinkageMatrix(std::vector<Node*> & population);

    bool CheckPopulationConverged(const std::vector<Node*>& population) override;

    std::pair<std::pair<std::vector<std::vector<Node*>>,std::vector<std::vector<Node*>>>,std::vector<size_t>>
    K_leader_means(std::vector<Node *> &population);

    arma::mat *  linkage_normalization_matrix;
    bool gomea_converged = false;

    MOArchive * mo_archive;



private:

};

#endif /* GOMEAMOGENERATIONHANDLER_H */
