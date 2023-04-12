/*
 


 */

/* 
 * File:   EvolutionRun.h
 * Author: virgolin
 *
 * Created on June 28, 2018, 11:28 AM
 */

#ifndef EVOLUTIONRUN_H
#define EVOLUTIONRUN_H


#include "GPGOMEA/Evolution/EvolutionState.h"
#include "GPGOMEA/Evolution/PopulationInitializer.h"
#include "GPGOMEA/Utils/Logger.h"
#include "GPGOMEA/GOMEA/GOMEAGenerationHandler.h"
#include "GPGOMEA/Evolution/GenerationHandler.h"
#include "GPGOMEA/Evolution/NSGA2GenerationHandler.h"
#include "GPGOMEA/Fitness/MOFitness.h"
#include "GPGOMEA/Evolution/MOArchive.h"

#include <iostream>
#include <vector>

class EvolutionRun {
public:

    EvolutionRun(EvolutionState & st) {
        config = new ConfigurationOptions(*st.config); // clone of the configuration settings
        if (st.semantic_library)
            semantic_library = new SemanticLibrary(*st.semantic_library); // clone semantic library

        if (st.config->gomea) {
            if (dynamic_cast<MOFitness*>(st.fitness)){
                generation_handler = new GOMEAMOGenerationHandler(*((GOMEAMOGenerationHandler*) st.generation_handler)); // clone generation handler
                ((GOMEAMOGenerationHandler*) generation_handler)->linkage_normalization_matrix= new arma::mat();// detach pointer to previous linkage normalization matrix
                ((GOMEAMOGenerationHandler*) generation_handler)->gomea_converged = false;
                ((GOMEAMOGenerationHandler*) generation_handler)->mo_archive = &mo_archive;
                is_gomea_and_multiobj = true;
            }else{
                generation_handler = new GOMEAGenerationHandler(*((GOMEAGenerationHandler*) st.generation_handler)); // clone generation handler
                if (((GOMEAGenerationHandler*) st.generation_handler)->linkage_normalization_matrix)
                    ((GOMEAGenerationHandler*) generation_handler)->linkage_normalization_matrix = new arma::mat(); // detach pointer to previous linkage normalization matrix
                ((GOMEAGenerationHandler*) generation_handler)->gomea_converged = false;
            }
        } else {
            if (dynamic_cast<MOFitness*>(st.fitness)){
                generation_handler = new NSGA2GenerationHandler(*dynamic_cast<NSGA2GenerationHandler*>(st.generation_handler));
                is_nsga_and_multiobj = true;
            }
            else
                generation_handler = new GenerationHandler(*st.generation_handler);
            if (st.semantic_library)
                generation_handler->semlib = semantic_library;
        }

        // set correct handler to config
        generation_handler->conf = config;

        tree_initializer = st.tree_initializer; // share same tree initializer
        fitness = st.fitness; // share same fitness
        mo_archive.fitness = fitness;
        mo_archive.config = config;
        if (dynamic_cast<MOFitness*>(st.fitness)){
            mo_archive.so_archive = std::vector<Node*>(dynamic_cast<MOFitness*>(st.fitness)->sub_fitness_functions.size(),
                                                       nullptr);
        }


        nr_trees = st.nr_trees;
    };

    virtual ~EvolutionRun() {
        for (Node * n : population) {
            n->ClearSubtree();
        }

        if (elitist)
            elitist->ClearSubtree();

        delete config;
        if (semantic_library)
            delete semantic_library;
        delete generation_handler;
    }

    void Initialize();
    void DoGeneration();

    std::vector<Node*> population;
    ConfigurationOptions * config = NULL;
    TreeInitializer * tree_initializer = NULL;
    GenerationHandler * generation_handler = NULL;
    Fitness * fitness = NULL;
    SemanticLibrary * semantic_library = NULL;
    size_t nr_trees;

    arma::vec pop_fitnesses;

    Node * elitist = NULL;
    double_t elitist_fit = arma::datum::inf;
    size_t elitist_size;

    size_t gen = 0;

    bool is_nsga_and_multiobj = false;
    bool is_gomea_and_multiobj = false;
    MOArchive mo_archive;
    std::vector<Node*> mo_archive_n;


private:

};

#endif /* EVOLUTIONRUN_H */
