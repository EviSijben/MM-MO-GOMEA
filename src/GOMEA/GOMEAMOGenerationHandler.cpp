//
// Created by evi on 4/26/21.
//

#include "GPGOMEA/GOMEA/GOMEAMOGenerationHandler.h"
#include <math.h>
//#include "matplotlib.h"
//
//namespace plt = matplotlibcpp;

void GOMEAMOGenerationHandler::PerformGeneration(std::vector<Node *> &population) {
    // if using GOMEA and it converged, re-initialize it
    if (!conf->use_IMS && gomea_converged) {
        // this is going to be handled by the IMS if it is active
        for (Node *n: population) {
            n->ClearSubtree();
        }
        population.clear();
        population = PopulationInitializer::InitializeTreePopulation(*conf, *tree_initializer, *fitness);
        gomea_converged = false;
    }

    std::vector<Node *> offspring;
    size_t nr_obj = population[0]->cached_objectives.size();
    size_t NIS = 1 + log10(conf->population_size);

    //Make clusters
    std::pair<std::pair<std::vector<std::vector<Node *>>, std::vector<std::vector<Node *>>>, std::vector<size_t>> output = K_leader_means(population);
    std::vector<std::vector<Node *>> clustered_population = output.first.first;
    std::vector<std::vector<Node *>> clustered_population_equal = output.first.second;
    std::vector<size_t> clusternr = output.second;


    // Build FOS for each cluster in population
    std::vector<std::vector<std::vector<size_t >>> FOSs;
    for (size_t i = 0; i < clustered_population.size(); i++) {
        FOSs.push_back(GOMEAFOS::GenerateFOS(clustered_population_equal[i], conf->fos_type, conf->gomfos_noroot,
                                             linkage_normalization_matrix));
    }

    // Variate solutions

    std::vector<std::pair<size_t, size_t>> idx;
    for (size_t i = 0; i < clustered_population.size(); i++) {
        for (size_t j = 0; j < clustered_population[i].size(); j++) {
            idx.emplace_back(i, j);
        }
    }

    std::mutex lock;
#pragma omp parallel for schedule(static)
    for (size_t x = 0; x <idx.size(); x++) {
        size_t &i = idx[x].first;
        size_t &j = idx[x].second;
        Node *off;
        if (clustered_population[i].size() > 1) {
            off = GOMVariator::GOMMO(clustered_population[i][j], clustered_population_equal[i], FOSs[i], *fitness,
                                     conf->caching, mo_archive, clusternr[i], clusternr[i] < nr_obj, NIS);
        } else{
            off = GOMVariator::GOMMO(*clustered_population[i][j], conf->functions, conf->terminals, FOSs[i], *fitness,
                                     conf->caching, mo_archive, clusternr[i], clusternr[i] < nr_obj, NIS);
        }

        std::lock_guard<std::mutex> plock(lock);
        offspring.push_back(off);
    }

    assert(offspring.size() == population.size());

    // Replace population
    for (size_t i = 0; i < population.size(); i++) {
        population[i]->ClearSubtree();
    }
    population = offspring;


}

void GOMEAMOGenerationHandler::InitLinkageMatrix(std::vector<Node*> & population){
    GOMEAFOS::GenerateFOS(population, conf->fos_type, conf->gomfos_noroot,linkage_normalization_matrix);
}


bool GOMEAMOGenerationHandler::CheckPopulationConverged(const std::vector<Node *> &population) {
    size_t idx = arma::randu() * population.size();
    std::vector<Node *> rand_sol_nodes = population[idx]->GetSubtreeNodes(false);
    std::vector<bool> active_sol_nodes;
    active_sol_nodes.reserve(rand_sol_nodes.size());
    for (size_t i = 0; i < rand_sol_nodes.size(); i++)
        active_sol_nodes.push_back(rand_sol_nodes[i]->IsActive());

    for (size_t i = 0; i < population.size(); i++) {
        if (i == idx)
            continue;

        std::vector<Node *> nodes = population[i]->GetSubtreeNodes(false);

        for (size_t j = 0; j < nodes.size(); j++) {
            if ((nodes[j]->IsActive() != active_sol_nodes[j]) ||
                nodes[j]->GetValue().compare(rand_sol_nodes[j]->GetValue()) != 0)
                return false;
        }
    }

    return true;

}

std::pair<std::pair<std::vector<std::vector<Node *>>, std::vector<std::vector<Node *>>>, std::vector<size_t>>
GOMEAMOGenerationHandler::K_leader_means(std::vector<Node *> &population) {
    size_t k = 7;

    size_t pop_size = population.size();
    // keep track of possible leaders
    std::vector<size_t> range(pop_size);
    std::iota(range.begin(), range.end(), 0);
    // all solutions are still possible as leaders
    std::unordered_set<size_t> remaining_solutions(range.begin(), range.end());

    // get objective value data and normalize
    size_t nr_objs = population[0]->cached_objectives.size();
    arma::Mat<double> norm_data = arma::Mat<double>(nr_objs, pop_size);
    for (size_t i = 0; i < pop_size; i++) {
        for (size_t j = 0; j < nr_objs; j++) {
            norm_data.at(j, i) = population[i]->cached_objectives[j];
        }
    }

    for (size_t i = 0; i < nr_objs; i++) {
        double min;
        bool min_initialized = false;
        double max;
        bool max_initialized = false;
        // get min and max for objective
        for (size_t x: remaining_solutions) {
            double_t value = norm_data.at(i, x);
            if (!min_initialized || value < min) {
                min = value;
                min_initialized = true;
            }
            if (!max_initialized || value > max) {
                max = value;
                max_initialized = true;
            }
        }
        //normalize each value in population for that objective
        for (size_t j = 0; j < pop_size; j++) {
            norm_data.at(i, j) = (norm_data.at(i, j) - min) / ((max - min)+0.0000000001);
        }
    }

    // get first leaders based on smallest objective value
    size_t initalised_k = 0;
    std::vector<size_t> idx_leaders;
    // take a random objective
    size_t random_obj = std::round(arma::randu()*(nr_objs-1));

    size_t idx_min = *remaining_solutions.begin();
    for (size_t x: remaining_solutions) {
        if(norm_data.at(random_obj,idx_min)>norm_data.at(random_obj,x)){
            idx_min = x;
        }
    }
    size_t first_leader = idx_min;
    idx_leaders.push_back(first_leader);
    remaining_solutions.erase(first_leader);
    initalised_k++;


    //find other leaders
    //first get distance of all remaining possible leaders to first leader
    arma::vec dists = arma::vec(pop_size, arma::fill::zeros);
    for (size_t x: remaining_solutions) {
        dists[x] = arma::mean(arma::square(norm_data.col(x) - norm_data.col(idx_leaders[0])));
    }


    while (initalised_k < k && !remaining_solutions.empty()) {
        // select leader with longest distance to other leaders
        size_t new_leader = dists.index_max();
        idx_leaders.push_back(new_leader);
        initalised_k++;
        // minimum distance of new leader to other leaders now becomes 0 and the new leader is not a remaining possible leader
        dists[new_leader] = 0;
        remaining_solutions.erase(new_leader);
        // update distances of possible new leaders: if the new solution is closer to the individual then the one currently closest we update the distance
        for (size_t x: remaining_solutions) {
            double minimum = std::min(arma::mean(arma::square(norm_data.col(x) - norm_data.col(new_leader))), dists[x]);
            dists[x] = minimum;
        }
    }

    // intialize cluster tags for k-means
    std::vector<size_t> clustertags_kmeans(pop_size);
    // intialize bool for whether there was a change in cluster assignments
    bool cluster_change = true;
    // initialize matrix for centers of the clusters (start with location of the leaders)
    arma::mat cluster_centers = arma::mat(nr_objs, idx_leaders.size(), arma::fill::zeros);
    for (size_t i = 0; i < idx_leaders.size(); i++) {
        cluster_centers.col(i) = norm_data.col(idx_leaders[i]);
    }

// k-means
    size_t iter = 0;
    while (cluster_change && iter < 50) {
        cluster_change = false;
        arma::mat cluster_centers_temp = arma::mat(nr_objs, idx_leaders.size(), arma::fill::zeros);
        std::vector<size_t> counts = std::vector<size_t>(idx_leaders.size(), 0);

        for (size_t i = 0; i < population.size(); i++) {
            bool lowest_dist_init = false;
            double lowest_dist;
            size_t index_min;
            // check the center that is closest for this indivdual
            for (size_t j = 0; j < idx_leaders.size(); j++) {
                if (!lowest_dist_init ||
                    lowest_dist > arma::mean(arma::square(norm_data.col(i) - cluster_centers.col(j)))) {
                    lowest_dist = arma::mean(arma::square(norm_data.col(i) - cluster_centers.col(j)));
                    index_min = j;
                    lowest_dist_init = true;
                }
            }
            // if it is different than before we update it
            if (index_min != clustertags_kmeans[i]) {
                clustertags_kmeans[i] = index_min;
                cluster_change = true;
            }
            // cluster_centers_temp are updated as well as the number of individuals belonging to a center, for calculating the new clustercenters later on
            cluster_centers_temp.col(clustertags_kmeans[i]) += norm_data.col(i);
            counts[clustertags_kmeans[i]] += 1;
        }
        for (size_t nr_lead = 0; nr_lead < idx_leaders.size(); nr_lead++) {
            // if the center doesn't have any individuals, reset to position leader
            if (counts[nr_lead] == 0) {
                cluster_centers.col(nr_lead) = norm_data.col(idx_leaders[nr_lead]);
                clustertags_kmeans[idx_leaders[nr_lead]] = nr_lead;
            }
                // else finalize new cluster center by calculating sum values of all individuals divided by the number of individuals of that cluster
            else {
                for (size_t obj = 0; obj < nr_objs; obj++) {
                    cluster_centers.at(obj, nr_lead) = cluster_centers_temp.at(obj, nr_lead) / counts[nr_lead];
                }
            }
        }
        iter++;

    }

    std::vector<std::vector<size_t>> clustertags_equal = std::vector<std::vector<size_t >>(pop_size,std::vector<size_t>());

    // store the solutions of each cluster in a separate vector for FOS in clustered_population
    std::vector<std::vector<Node *>> clustered_population = std::vector<std::vector<Node * >>(initalised_k,
                                                                                              std::vector<Node *>());
    std::vector<std::vector<Node *>> clustered_population_equal = std::vector<std::vector<Node * >>(initalised_k,
                                                                                                    std::vector<Node *>());
    // assign to each cluster the closests 2*pop_size/cluster solutions
    for (size_t x = 0; x < initalised_k; x++) {
        arma::vec distances = arma::vec(pop_size, arma::fill::none);
        for (size_t i = 0; i < pop_size; i++) {
            distances[i] = arma::mean(arma::square(norm_data.col(i) - cluster_centers.col(x) ));
        }
        for (size_t times = 0; times < ((2 * pop_size) / initalised_k); times++) {
            clustered_population_equal[x].push_back(population[distances.index_min()]);
            clustertags_equal[distances.index_min()].push_back(x);
            distances[distances.index_min()] = arma::datum::inf;
        }

    }

    // define extreme clusters
    std::vector<size_t> clusternr = std::vector<size_t>(initalised_k, nr_objs + 1);

    for (size_t i = 0; i < nr_objs; i++)
        clusternr[cluster_centers.row(i).index_min()] = i;


    // if not jet assigned solution in equal size clustering, assign to closest. if multiple are assigned, assign to random of multiple center
    for (size_t i = 0; i < pop_size; i++) {
        if (!clustertags_equal[i].empty()) {
            size_t assign;
            if (clustertags_equal[i].size() == 1) {
                assign = clustertags_equal[i][0];
            } else {
                assign = clustertags_equal[i][arma::randu() * clustertags_equal[i].size()];
            }
            clustered_population[assign].push_back(population[i]);
        }else{
            clustered_population[clustertags_kmeans[i]].push_back(population[i]);
        }

    }



    return std::make_pair(std::make_pair(clustered_population, clustered_population_equal), clusternr);


}


