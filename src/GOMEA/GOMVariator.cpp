/*



 */

/*
 * File:   GOMVariator.cpp
 * Author: virgolin
 *
 * Created on June 28, 2018, 6:07 PM
 */



#include "GPGOMEA/GOMEA/GOMVariator.h"

#include <random>
#include "GPGOMEA/GOMEA/GOMEAFOS.h"
#include "GPGOMEA/Evolution/MOArchive.h"

using namespace std;
using namespace arma;

static bool MeaningfulChange(const vector<size_t> &changed_indices, const vector<Node *> &offspring_nodes) {
    for (size_t j: changed_indices) {
        Node *off_node = offspring_nodes[j];
        if (off_node->IsActive()) {
            return true;
        }
    }
    return false;
}

static void ClearCache(const vector<size_t> &changed_indices, vector<Node *> &offspring_nodes) {
    for (size_t j: changed_indices) {
        offspring_nodes[j]->ClearCachedOutput(true);
    }
}

static inline void RevertChanges(const vector<size_t> &changed_indices, vector<Node *> &offspring_nodes,
                                 vector<Operator *> &backup_operators, bool use_caching) {
    for (size_t j = 0; j < changed_indices.size(); j++) {
        auto off_node = offspring_nodes[changed_indices[j]];
        Operator *op = off_node->ChangeOperator(*backup_operators[j]);
        delete op;
        if (use_caching)
            off_node->ClearCachedOutput(true);
    }
}

static void ChangeTree(const vector<size_t> &F, vector<Node *> &offspring_nodes, const vector<Node *> &donor_nodes,
                       vector<Operator *> &backup_operators, vector<size_t> &changed_indices) {
    for (size_t j: F) {
        auto off_node = offspring_nodes[j];
        auto don_node = donor_nodes[j];
        if (off_node->GetValue().compare(don_node->GetValue()) == 0) {
            continue;
        } else {
            Operator *replaced_op = off_node->ChangeOperator(*don_node->GetOperator());
            backup_operators.push_back(replaced_op);
            changed_indices.push_back(j);
        }
    }
}

static void ChangeTreeWithoutDonor(const vector<size_t> &F, vector<Node *> &offspring_nodes,
                                   vector<Operator *> &backup_operators, vector<size_t> &changed_indices,const std::vector<Operator *> &functions, const std::vector<Operator *> &terminals) {
    for (size_t j: F) {
        auto off_node = offspring_nodes[j];
        size_t h = off_node->GetHeight(false);

        Operator *don_op;
        if (h > 0 && randu() < 0.5) {
            don_op = functions[randu() * functions.size()];
        } else {
            don_op = terminals[randu() * terminals.size()];
        }

        if (don_op->name.compare(off_node->GetValue()) == 0) {
            continue;
        }

        Operator *replaced_op = off_node->ChangeOperator(*don_op);
        backup_operators.push_back(replaced_op);
        changed_indices.push_back(j);
    }
}

static vector<size_t> GetPermutedFOSIndices(size_t fos_size) {
    vector<size_t> permut_fos_indices;
    for (size_t i = 0; i < fos_size; i++) {
        permut_fos_indices.push_back(i);
    }
    std::random_shuffle(permut_fos_indices.begin(), permut_fos_indices.end());
    return permut_fos_indices;
}

Node *GOMVariator::GOMMO(const Node &sol,const std::vector<Operator *> &functions, const std::vector<Operator *> &terminals, const std::vector<std::vector<size_t> > &FOS,
                         Fitness &fit, bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem,
                         size_t NIS_const) {
    Node *offspring = sol.CloneSubtree();
    arma::vec back_obj = sol.cached_objectives;

    bool improved = false;
    bool changed = false;

    auto permut_fos_indices = GetPermutedFOSIndices(FOS.size());

    vector<Node *> offspring_nodes = offspring->GetSubtreeNodes(false);
    for (size_t i = 0; i < FOS.size(); i++) {
        vector<size_t> F = FOS[permut_fos_indices[i]];

        vector<Operator *> backup_operators;
        vector<size_t> changed_indices;

        // change the nodes of the offspring with random operators and terminals
        ChangeTreeWithoutDonor(F,offspring_nodes,backup_operators,changed_indices,functions,terminals);

        // check for meaningful changes (done later because situation may change during iterations)
        bool meaningful_change = MeaningfulChange(changed_indices, offspring_nodes);

        // clean the cache of node outputs (even in case of not meaningful changes, for the future)
        if (use_caching) {
            ClearCache(changed_indices, offspring_nodes);
        }

        // determine whether to keep or not the changes
        if (meaningful_change) {
            fit.ComputeFitness(offspring, use_caching);
            std::pair<bool, bool> check_changes = extreem ? check_changes_SO(offspring, back_obj, false, obj)
                                                          : check_changes_MO(offspring, back_obj, false, mo_archive);
            if(check_changes.second){
                improved = true;
            }
            if (check_changes.first) {
                // keep changes
                back_obj = offspring->cached_objectives;
                changed = true;
                mo_archive->UpdateMOArchive(offspring);
                mo_archive->UpdateSOArchive(offspring);
            } else {
                // revert
                RevertChanges(changed_indices, offspring_nodes, backup_operators, use_caching);
                offspring->cached_fitness = back_obj[0];
                offspring->cached_objectives = back_obj;
            }
        }else{
            offspring->cached_fitness = back_obj[0];
            offspring->cached_objectives = back_obj;
        }

        // discard backup
        for (Operator *op: backup_operators) {
            delete op;
        }
    }
    offspring->NIS = improved ? 0 : offspring->NIS + 1;

    if ((!changed || offspring->NIS > NIS_const) && !mo_archive->mo_archive.empty()) {
        offspring->NIS = 0;
        return GOMVariator::GOMMOFI(offspring,  FOS, fit, use_caching, mo_archive, obj, extreem, NIS_const);
    }

    return offspring;
}

Node *GOMVariator::GOMMO(Node *&sol, const std::vector<Node *> &donors, const std::vector<std::vector<size_t> > &FOS,
                         Fitness &fit, bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem,
                         size_t NIS_const) {
    Node *offspring = sol->CloneSubtree();

    arma::vec back_obj = sol->cached_objectives;

    auto permut_fos_indices = GetPermutedFOSIndices(FOS.size());

    vector<Node *> offspring_nodes = offspring->GetSubtreeNodes(false);

    bool changed = false;
    bool improved = false;

    for (size_t i = 0; i < FOS.size(); i++) {
        const vector<size_t> &F = FOS[permut_fos_indices[i]];
        auto donor_nodes = donors[randu() * donors.size()]->GetSubtreeNodes(false);

        // change the nodes of the offspring with the ones of the donor
        vector<Operator *> backup_operators;
        vector<size_t> changed_indices;
        ChangeTree(F, offspring_nodes, donor_nodes, backup_operators, changed_indices);

        // check for meaningful changes (done later because situation may change during iterations)
        bool meaningful_change = MeaningfulChange(changed_indices, offspring_nodes);

        // clean the cache of node outputs (even in case of not meaningful changes, for the future)
        if (use_caching) {
            ClearCache(changed_indices, offspring_nodes);
        }

        if (meaningful_change) {
            fit.ComputeFitness(offspring, use_caching);
            std::pair<bool, bool> check_changes = extreem ? check_changes_SO(offspring, back_obj,  false, obj)
                                                          : check_changes_MO(offspring, back_obj, false, mo_archive);
            if(check_changes.second){
                improved = true;
            }
            if (check_changes.first) {
                // keep changes
                back_obj = offspring->cached_objectives;
                changed = true;
                mo_archive->UpdateMOArchive(offspring);
                mo_archive->UpdateSOArchive(offspring);
            } else {
                RevertChanges(changed_indices,offspring_nodes,backup_operators,use_caching);
                offspring->cached_fitness = back_obj[0];
                offspring->cached_objectives = back_obj;
            }
        }else{
            offspring->cached_fitness = back_obj[0];
            offspring->cached_objectives = back_obj;

        }



        // discard backup
        for (Operator *op: backup_operators) {
            delete op;
        }

    }

    offspring->NIS = improved ? 0 : offspring->NIS + 1;

    if ((!changed || offspring->NIS > NIS_const) && !mo_archive->mo_archive.empty()) {
        offspring->NIS = 0;
        return GOMVariator::GOMMOFI(offspring, FOS, fit, use_caching, mo_archive, obj, extreem, NIS_const);
    }
    return offspring;
}

Node *GOMVariator::GOMMOFI(Node *& offspring, const std::vector<std::vector<size_t> > &FOS,
                           Fitness &fit, bool use_caching, MOArchive *mo_archive, size_t obj, bool extreem,
                           size_t NIS_const) {
    arma::vec back_obj = offspring->cached_objectives;

    auto permut_fos_indices = GetPermutedFOSIndices(FOS.size());

    vector<Node *> offspring_nodes = offspring->GetSubtreeNodes(false);

    bool changed = false;

    for (size_t i = 0; i < FOS.size(); i++) {
        vector<size_t> F = FOS[permut_fos_indices[i]];
        auto donor = extreem ? mo_archive->ReturnCopySOMember(obj) : mo_archive->ReturnCopyRandomMOMember();
        auto donor_nodes = donor->GetSubtreeNodes(false);

        // change the nodes of the offspring with the ones of the donor
        vector<Operator *> backup_operators;
        vector<size_t> changed_indices;
        ChangeTree(F, offspring_nodes, donor_nodes, backup_operators, changed_indices);

        // check for meaningful changes (done later because situation may change during iterations)
        bool meaningful_change = MeaningfulChange(changed_indices, offspring_nodes);

        // clean the cache of node outputs (even in case of not meaningful changes, for the future)
        if (use_caching) {
            ClearCache(changed_indices, offspring_nodes);
        }

        if (meaningful_change) {
            fit.ComputeFitness(offspring, use_caching);
            std::pair<bool, bool> check_changes = extreem ? check_changes_SO(offspring, back_obj, true, obj)
                                                          : check_changes_MO(offspring, back_obj, true, mo_archive);
            if (check_changes.first) {
                // keep changes
                back_obj = offspring->cached_objectives;
                changed = true;
                mo_archive->UpdateMOArchive(offspring);
                mo_archive->UpdateSOArchive(offspring);
            } else {
                // revert
                RevertChanges(changed_indices, offspring_nodes, backup_operators, use_caching);
                offspring->cached_fitness = back_obj[0];
                offspring->cached_objectives = back_obj;
            }
        }else{
            offspring->cached_fitness = back_obj[0];
            offspring->cached_objectives = back_obj;
        }


        // discard backup
        for (Operator *op: backup_operators) {
            delete op;
        }
        donor->ClearSubtree();

        if (changed) {
            break;
        }
    }

    if ((!changed ) && !mo_archive->mo_archive.empty()) {
        offspring->ClearSubtree();
        if (extreem)
            return mo_archive->ReturnCopySOMember(obj);
        else return mo_archive->ReturnCopyRandomMOMember();
    }

    return offspring;

}

std::pair<bool, bool>
GOMVariator::check_changes_SO(Node *offspring, arma::vec back_obj, bool FI, size_t obj) {
    if (offspring->cached_objectives[obj] < back_obj[obj] ) {
        return std::make_pair(true,true);
    }
    if (offspring->cached_objectives[obj] == back_obj[obj] ) {
        if (FI) {
            return std::make_pair(false,false);
        } else {
            return std::make_pair(true,false);
        }
    }

    return std::make_pair(false, false);
}

std::pair<bool, bool>
GOMVariator::check_changes_MO(Node *offspring, arma::vec back_obj,  bool FI,MOArchive *mo_archive) {
    bool dominates = false;

    for (size_t j = 0; j < offspring->cached_objectives.n_elem; j++) {
        if (offspring->cached_objectives[j] < back_obj[j])
            dominates = true;
        else if (offspring->cached_objectives[j] > back_obj[j]) {
            dominates = false;
            break;
        }
    }
    if (dominates) {
        return std::make_pair(true, true);
    }
    // or they have the same objective values (in case of not FI)
    if (!FI) {
        bool same = true;
        for (size_t j = 0; j < offspring->cached_objectives.n_elem; j++) {
            if (offspring->cached_objectives[j] != back_obj[j]) {
                same = false;
                break;
            }
        }
        if (same)
            return std::make_pair(true,false);
    }

    // or is not dominated by mo_archive
    if (mo_archive->NonDominated(offspring)) {
        mo_archive->UpdateMOArchive(offspring);
        return std::make_pair(true,true);

    }
    return std::make_pair(false,false);
}

Node *GOMVariator::GOM(const Node &sol, const std::vector<Node *> &donors, const std::vector<std::vector<size_t>> &FOS,
                       Fitness &fit, bool use_caching) {
    Node *offspring = sol.CloneSubtree();
    double_t back_fit = sol.cached_fitness;

    auto permut_fos_indices = GetPermutedFOSIndices(FOS.size());

    vector<Node *> offspring_nodes = offspring->GetSubtreeNodes(false);
    for (size_t i = 0; i < FOS.size(); i++) {
        vector<size_t> F = FOS[permut_fos_indices[i]];
        auto donor_nodes = donors[randu() * donors.size()]->GetSubtreeNodes(false);
        // change the nodes of the offspring with the ones of the donor
        vector<Operator *> backup_operators;
        vector<size_t> changed_indices;
        ChangeTree(F, offspring_nodes, donor_nodes, backup_operators, changed_indices);

        // check for meaningful changes (done later because situation may change during iterations)
        bool meaningful_change = MeaningfulChange(changed_indices, offspring_nodes);

        // clean the cache of node outputs (even in case of not meaningful changes, for the future)
        if (use_caching) {
            ClearCache(changed_indices, offspring_nodes);
        }

        // determine whether to keep or not the changes
        bool keep_changes = true;
        if (meaningful_change) {
            // check fitness it not worse
            double_t new_fit = fit.ComputeFitness(offspring, use_caching);
            keep_changes =  new_fit <= back_fit ;
        }

        if (keep_changes) {
            back_fit = offspring->cached_fitness;
        } else {
            RevertChanges(changed_indices, offspring_nodes, backup_operators, use_caching);
            offspring->cached_fitness = back_fit;
        }

        // discard backup
        for (Operator *op: backup_operators) {
            delete op;
        }
    }

    return offspring;
}

Node *
GOMVariator::GOM(const Node &sol, const std::vector<Operator *> &functions, const std::vector<Operator *> &terminals,
                 const std::vector<std::vector<size_t> > &FOS, Fitness &fit, bool use_caching) {

    Node *offspring = sol.CloneSubtree();
    double_t back_fit = sol.cached_fitness;

    auto permut_fos_indices = GetPermutedFOSIndices(FOS.size());

    vector<Node *> offspring_nodes = offspring->GetSubtreeNodes(false);
    for (size_t i = 0; i < FOS.size(); i++) {
        vector<size_t> F = FOS[permut_fos_indices[i]];

        // change the nodes of the offspring with the ones of the donor
        vector<Operator *> backup_operators;
        vector<size_t> changed_indices;

        ChangeTreeWithoutDonor(F,offspring_nodes,backup_operators,changed_indices,functions,terminals);

        // check for meaningful changes (done later because situation may change during iterations)
        bool meaningful_change = MeaningfulChange(changed_indices, offspring_nodes);

        // clean the cache of node outputs (even in case of not meaningful changes, for the future)
        if (use_caching) {
            ClearCache(changed_indices, offspring_nodes);
        }

        // determine whether to keep or not the changes
        bool keep_changes = true;
        if (meaningful_change) {
            // check fitness it not worse
            double_t new_fit = fit.ComputeFitness(offspring, use_caching);
            keep_changes = new_fit <= back_fit;
        }


        if (keep_changes) {
            // keep changes
            back_fit = offspring->cached_fitness;
        } else {
            // revert
            RevertChanges(changed_indices, offspring_nodes, backup_operators, use_caching);
            offspring->cached_fitness = back_fit;
        }

        // discard backup
        for (Operator *op: backup_operators) {
            delete op;
        }
    }
    return offspring;
}


Node *
GOMVariator::MakeBiggerGOMEATree(const Node &original, size_t orig_height, size_t desired_height, size_t max_arity,
                                 const GOMEATreeInitializer &tree_init, const std::vector<Operator *> &functions,
                                 const std::vector<Operator *> &terminals) {

    assert(desired_height > orig_height);

    Node *bigger = original.CloneSubtree();

    vector<Node *> nodes = bigger->GetSubtreeNodes(false);

    int height_left = desired_height - orig_height - 1;

    for (Node *n: nodes) {
        // take leaves
        if (n->GetDepth() == orig_height) {
            // append stuff
            Node *to_append;
            for (size_t i = 0; i < max_arity; i++) {
                if (height_left == 0) {
                    to_append = new SingleNode(terminals[randu() * terminals.size()]->Clone());
                } else {
                    to_append = tree_init.GenerateTreeGrow(height_left, height_left, 0, max_arity, functions,
                                                           terminals);
                }
                n->AppendChild(to_append);
            }
        }
    }

    size_t n = bigger->GetSubtreeNodes(false).size();
    size_t m = pow(2, desired_height + 1) - 1;
    assert(n == m);

    return bigger;

}
