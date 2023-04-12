//
// Created by evi on 11/3/21.
//

#ifndef GPGOMEA_MOARCHIVE_H
#define GPGOMEA_MOARCHIVE_H

#include <boost/thread.hpp>
#include "GPGOMEA/Genotype/Node.h"
#include "GPGOMEA/Genotype/Multitree.h"

class MOArchive {
public:
    Fitness *fitness;
    ConfigurationOptions *config;

    std::vector<Node *> mo_archive;
    std::vector<Node *> so_archive;

private:
    boost::shared_mutex mo_lock;
    boost::shared_mutex so_lock;

public:
    void UpdateSOArchive(Node *offspring) ;
    void UpdateMOArchive(Node *offspring) ;

    Node *ReturnCopyRandomMOMember() ;

    Node *ReturnCopySOMember(size_t idx) ;

    bool NonDominated(Node *offspring) ;

    bool DiversityAdded(Node *offspring, size_t idx) ;

    void InitSOArchive(std::vector<Node *> population) ;
    void InitMOArchive(std::vector<Node *> population) ;
    void SaveResults(size_t gen) ;


private:
};

#endif //GPGOMEA_MOARCHIVE_H



