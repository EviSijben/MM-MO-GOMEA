/*
 


 */

/* 

 */

#ifndef ENSAMBLEFITNESS_H
#define ENSAMBLEFITNESS_H

#include "GPGOMEA/Fitness/Fitness.h"

#include <armadillo>

class DiversifiedSymbolicRegressionFitness : public Fitness {
public:

    double_t ComputeFitness(Node * n, bool use_caching) override;
    
    double_t GetTestFit(Node * n) override;
    
    double_t GetValidationFit(Node * n) override;



private:

    double_t ComputeMSE(const arma::vec & res);





};

#endif /* ENSAMBLEFITNESS_H */

