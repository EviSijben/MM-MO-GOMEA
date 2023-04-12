/*



 */




#include "GPGOMEA/Fitness/DiversifiedSymbolicRegressionFitness.h"

double_t DiversifiedSymbolicRegressionFitness::ComputeFitness(Node* n, bool use_caching) {

    evaluations++;

    arma::mat P = n->GetOutput(TrainX, use_caching);
    double fit;
    if (n->type == NodeType::Multi) {
        for (size_t i = 0; i < size(P)[1]; i++) {
            P.col(i) = arma::abs(P.col(i) - TrainY);
        }
        fit = arma::mean(arma::square((arma::min(P, 1))));
    }else{
        throw std::runtime_error("DiversifiedSymbolicRegressionFitness::Should not use Diversified fitness without multitree");
    }
    if (std::isnan(fit)) {
        fit = arma::datum::inf;
    }
    n->cached_fitness = fit;
    return fit;

}

double_t DiversifiedSymbolicRegressionFitness::GetValidationFit(Node* n) {
    arma::mat P = n->GetOutput(ValidationX, false);
    double fit;
    for (size_t i = 0; i < size(P)[1]; i++) {
        P.col(i) = abs(P.col(i) - ValidationY);
    }
    fit = ComputeMSE(arma::min(P, 1));

    return fit;
}


double_t DiversifiedSymbolicRegressionFitness::GetTestFit(Node* n) {
    arma::mat P = n->GetOutput(TestX, false);
    double fit;

    for (size_t i = 0; i < size(P)[1]; i++) {
        P.col(i) = abs(P.col(i) - TestY);
    }
    fit = ComputeMSE(arma::min(P, 1));


    return fit;

}


double_t DiversifiedSymbolicRegressionFitness::ComputeMSE(const arma::vec& res) {
    double_t mse = arma::mean( arma::square(res) );
    return mse;
}
