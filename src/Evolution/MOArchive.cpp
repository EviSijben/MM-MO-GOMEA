//
// Created by evi on 11/3/21.
//

#include "GPGOMEA/Evolution/MOArchive.h"

void MOArchive::UpdateSOArchive(Node *offspring) {
    boost::unique_lock<boost::shared_mutex> write_guard(so_lock);
    for (size_t i = 0; i < offspring->cached_objectives.size(); i++) {
        if (so_archive[i] == nullptr) {
            Node *new_node = offspring->CloneSubtree();
            new_node->cached_fitness = offspring->cached_fitness;
            new_node->cached_objectives = offspring->cached_objectives;
            so_archive[i] = new_node;    // clone it

        } else if (so_archive[i]->cached_objectives[i] > offspring->cached_objectives[i]) {
            so_archive[i]->ClearSubtree();
            so_archive[i] = nullptr;
            Node *new_node = offspring->CloneSubtree();
            new_node->cached_fitness = offspring->cached_fitness;
            new_node->cached_objectives = offspring->cached_objectives;
            so_archive[i] = new_node;    // clone it
        }
    }

}

void MOArchive::UpdateMOArchive(Node *offspring) {
    boost::unique_lock<boost::shared_mutex> write_guard(mo_lock);
    bool solution_is_dominated = false;
    bool diversity_added = false;
    bool identical_objectives_already_exist;
    for (size_t i = 0; i < mo_archive.size(); i++) {
        // check domination
        solution_is_dominated = mo_archive[i]->Dominates(offspring);
        if (solution_is_dominated)
            break;

        identical_objectives_already_exist = true;
        for (size_t j = 0; j < offspring->cached_objectives.n_elem; j++) {
            if (offspring->cached_objectives[j] != mo_archive[i]->cached_objectives[j]) {
                identical_objectives_already_exist = false;
                break;
            }
        }
        if (identical_objectives_already_exist) {
            if (DiversityAdded(offspring, i)) {
                diversity_added = true;
                mo_archive[i]->ClearSubtree();
                mo_archive[i] = nullptr;
            }
            break;
        }

        if (offspring->Dominates(mo_archive[i])) {
            mo_archive[i]->ClearSubtree();
            mo_archive[i] = nullptr;
        }
    }
    mo_archive.erase(
            std::remove_if(mo_archive.begin(), mo_archive.end(), [](Node *node) { return node == nullptr; }),
            mo_archive.end());

    if ((!solution_is_dominated && !identical_objectives_already_exist) || (diversity_added)) {
        Node *new_node = offspring->CloneSubtree();
        fitness->ComputeFitness(new_node, false);
        mo_archive.push_back(new_node);    // clone it
    }

}


Node * MOArchive::ReturnCopyRandomMOMember() {
    boost::shared_lock<boost::shared_mutex> read_guard(mo_lock);
    size_t index = arma::randu() * mo_archive.size();
    Node *copy = mo_archive[index]->CloneSubtree();
    copy->cached_fitness = mo_archive[index]->cached_fitness;
    copy->cached_objectives = mo_archive[index]->cached_objectives;
    return copy;
}

Node * MOArchive::ReturnCopySOMember(size_t idx) {
    boost::shared_lock<boost::shared_mutex> read_guard(so_lock);
    Node *copy = so_archive[idx]->CloneSubtree();
    copy->cached_fitness = so_archive[idx]->cached_fitness;
    copy->cached_objectives = so_archive[idx]->cached_objectives;
    return copy;
}


bool MOArchive::NonDominated(Node *offspring) {
    bool solution_is_dominated = false;
    bool identical_objectives_already_exist;
    bool diversity_added = false;
    {
        boost::shared_lock<boost::shared_mutex> read_guard(mo_lock);
        if (mo_archive.empty())
            return true;
        for (size_t i = 0; i < mo_archive.size(); i++) {
            // check domination
            solution_is_dominated = mo_archive[i]->Dominates(offspring);
            if (solution_is_dominated)
                break;

            identical_objectives_already_exist = true;
            for (size_t j = 0; j < offspring->cached_objectives.n_elem; j++) {
                if (offspring->cached_objectives[j] != mo_archive[i]->cached_objectives[j]) {
                    identical_objectives_already_exist = false;
                    break;
                }
            }
            if (identical_objectives_already_exist) {
                if (DiversityAdded(offspring, i)) {
                    diversity_added = true;
                }
                break;
            }
        }
    }
    if (diversity_added) {
        UpdateMOArchive(offspring);
    }
    return (!solution_is_dominated && !identical_objectives_already_exist) || (diversity_added);
}

bool MOArchive::DiversityAdded(Node *offspring, size_t idx) {
    if (arma::mean(arma::mean(arma::square(
            offspring->GetOutput(fitness->TrainX, false) - mo_archive[idx]->GetOutput(fitness->TrainX, false)))) ==
        0) {
        return false;
    }else return true;
}


void MOArchive::InitSOArchive(std::vector<Node *> population) {
    for (Node *n: population) {
        UpdateSOArchive(n);
    }

}

void MOArchive::InitMOArchive(std::vector<Node *> population) {
    for (Node *n: population) {
        UpdateMOArchive(n);
    }

}
void MOArchive::SaveResults(size_t gen) {
    boost::unique_lock<boost::shared_mutex> read_guard(mo_lock);
    std::string results_path= config->results_path;
    if (std::to_string(gen).size() == 1)
        results_path += "/mo_archive_gen00" + std::to_string(gen) + ".csv";
    if (std::to_string(gen).size() == 2)
        results_path += "/mo_archive_gen0" + std::to_string(gen) + ".csv";
    if (std::to_string(gen).size() == 3)
        results_path += "/mo_archive_gen" + std::to_string(gen) + ".csv";
    std::ofstream myfile(results_path, std::ios::trunc);
    myfile << "nr" << "|" << "obj1_train" << "|" << "obj2_train" << "|" << "obj1_test" << "|" << "obj2_test" << "|" << "exp1" << "|" << "exp2" << std::endl;
    for (size_t p = 0; p < mo_archive.size(); p++) {
        fitness->GetTestFit(mo_archive[p]);
        if ((mo_archive[0])->type == NodeType::Multi) {
            myfile << p << "| " << mo_archive[p]->cached_objectives[0] << "| " << mo_archive[p]->cached_objectives[1]<< "| " << mo_archive[p]->cached_objectives_test[0]<< "| " << mo_archive[p]->cached_objectives_test[1];
            for (Node *k: (((Multitree *) mo_archive[p])->nodes)) {
                myfile << "|" << k->GetPythonExpression();
            }
            myfile << std::endl;

        } else {
            myfile << p << "| " << mo_archive[p]->cached_objectives[0] << "| " << mo_archive[p]->cached_objectives[1]<< "| " << mo_archive[p]->cached_objectives_test[0]<< "| " << mo_archive[p]->cached_objectives_test[1]
                   << mo_archive[p]->GetPythonExpression() << std::endl;

        }


    }

}



