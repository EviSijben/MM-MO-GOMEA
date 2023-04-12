/*
 


 */

/*
 * File:   SymbolicRegressionFitness.cpp
 * Author: virgolin
 *
 * Created on June 27, 2018, 6:23 PM
 */


#include "GPGOMEA/Fitness/SymbolicRegressionFitness.h"

double_t SymbolicRegressionFitness::ComputeFitness(Node* n, bool use_caching) {

    evaluations++;

    arma::mat P = n->GetOutput(TrainX, use_caching);
    double fit;
    if (n->type == NodeType::Multi){
        size_t nr_trees = size(P)[1];
        fit = 0;
        for (size_t i = 0; i < nr_trees; i++) {
            fit += ComputeMSE(P.col(i), TrainY);
        }
    }
    else{
        fit = ComputeMSE(P, TrainY);
    }
    if (std::isnan(fit)) {
        fit = arma::datum::inf;
    }
    n->cached_fitness = fit;
    return fit;

}

double_t SymbolicRegressionFitness::GetValidationFit(Node* n) {
    arma::mat P = n->GetOutput(ValidationX, false);
    double fit;
    if (n->type == NodeType::Multi){
        size_t nr_trees = size(P)[1];
        fit = 0;
        for (size_t i = 0; i < nr_trees; i++){
            fit += ComputeMSE(P.col(i), ValidationY);
        }
    }
    else{
        fit = ComputeMSE(P, ValidationY);
    }
    return fit;
}


double_t SymbolicRegressionFitness::GetTestFit(Node* n) {
    arma::mat P = n->GetOutput(TestX, false);
    double fit;
    if (n->type == NodeType::Multi){
        size_t nr_trees = size(P)[1];
        fit = 0;
        for (size_t i = 0; i < nr_trees; i++){
            fit += ComputeMSE(P.col(i), TestY);
        }
    }
    else{
        fit = ComputeMSE(P, TestY);
    }
    return fit;

}

double_t SymbolicRegressionFitness::ComputeMSE(const arma::vec& P, const arma::vec& Y) {
    arma::vec res = Y - P;
    double_t mse = arma::mean( arma::square(res) );
    return mse;
}

