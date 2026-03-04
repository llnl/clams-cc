#pragma once

#include "clustering_metrics_utils.hpp"

/**
 * Process point_to_clusters_map to fill cluster_overlap_map for calculating
 * non-information theoretic clustering metrics
 *
 * @param point_to_clusters_map A YGM map for point_id -> (cluster_id in
 * clustering 1, cluster_id in clustering 2).
 * @param cluster_overlap_map An empty YGM map of cluster pair (i,j) ->
 * tuple(overlap size, size of cluster i, size of cluster j), where i is a
 * cluster_id in the first clustering, j is a cluster_id in the second
 * clustering, overlap size is the size of the intersection between clusters i
 * and j. This is the cluster overlap map we use when calculating clustering
 * metrics with information theoretic metrics.
 * @param cluster_size_map1 An empty YGM map for cluster_id -> cluster size in
 * clustering 1.
 * @param cluster_size_map2 An empty YGM map for cluster_id -> cluster size in
 * clustering 2.
 *
 */
void fill_cluster_overlap_and_size_maps_mi(
    ygm::container::map<point_id_type,
                        std::pair<cluster_id_type, cluster_id_type>>
        &point_to_clusters_map,
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                        std::tuple<uint64_t, uint64_t, uint64_t>>
                                                   &cluster_overlap_map,
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map1,
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map2) {
  ygm::comm &mpi_comm                = point_to_clusters_map.comm();
  auto       cluster_overlap_map_ptr = cluster_overlap_map.get_ygm_ptr();

  // Make sure maps are empty to start
  cluster_overlap_map.clear();
  cluster_size_map1.clear();
  cluster_size_map2.clear();

  /* Get cluster sizes for each clustering and the overlap sizes
     Only count points that are not labeled as noise in clustering 2 */
  auto process_point_lambda =
      [&cluster_overlap_map, &cluster_size_map1, &cluster_size_map2](
          [[maybe_unused]] const point_id_type              &point_id,
          const std::pair<cluster_id_type, cluster_id_type> &cluster_pair) {
        // Lambda to update the size of the overlap
        // Takes cluster1 and cluster2 and updates cluster1 -> (cluster,
        // overlap += 1)
        auto update_overlap_lambda =
            []([[maybe_unused]] const std::pair<cluster_id_type,
                                                cluster_id_type> &cluster_pair,
               std::tuple<uint64_t, uint64_t, uint64_t> &overlap_tuple) {
              std::get<0>(overlap_tuple) = std::get<0>(overlap_tuple) + 1;
            };

        // Lambda to increase the cluster size in cluster_size_map 1&2
        auto increment_cluster_size_lambda =
            []([[maybe_unused]] const cluster_id_type &cluster_id,
               uint64_t &cluster_size) { cluster_size += 1; };

        // Increment the overlap count for cluster pair (i,j)
        // Only consider clusters that are not noise (i.e., cluster_id >=
        // 0)
        if (!is_cluster_noise(cluster_pair.second)) {
          cluster_overlap_map.async_visit(cluster_pair, update_overlap_lambda);
          cluster_size_map1.async_visit(cluster_pair.first,
                                        increment_cluster_size_lambda);
          cluster_size_map2.async_visit(cluster_pair.second,
                                        increment_cluster_size_lambda);
        }
        // cluster_overlap_map.async_visit(cluster_pair,
        // update_overlap_lambda);
      };

  point_to_clusters_map.for_all(process_point_lambda);
  mpi_comm.barrier();

  auto get_cluster_size_lambda =
      [&cluster_size_map1, &cluster_size_map2, &cluster_overlap_map_ptr](
          const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
          std::tuple<uint64_t, uint64_t, uint64_t>          &overlap_tuple) {
        auto get_size1_lambda =
            []([[maybe_unused]] const cluster_id_type            &cluster1,
               const uint64_t                                    &cluster_size,
               const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
               auto cluster_overlap_map_ptr) {
              auto assign_size1_lambda =
                  []([[maybe_unused]] const std::pair<
                         cluster_id_type, cluster_id_type>    &cluster_pair,
                     std::tuple<uint64_t, uint64_t, uint64_t> &overlap_tuple,
                     const uint64_t                            cluster_size) {
                    std::get<1>(overlap_tuple) = cluster_size;
                  };
              cluster_overlap_map_ptr->async_visit(
                  cluster_pair, assign_size1_lambda, cluster_size);
            };
        cluster_size_map1.async_visit(cluster_pair.first, get_size1_lambda,
                                      cluster_pair, cluster_overlap_map_ptr);

        auto get_size2_lambda =
            []([[maybe_unused]] const cluster_id_type            &cluster2,
               const uint64_t                                    &cluster_size,
               const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
               auto cluster_overlap_map_ptr) {
              auto assign_size2_lambda =
                  []([[maybe_unused]] const std::pair<
                         cluster_id_type, cluster_id_type>    &cluster_pair,
                     std::tuple<uint64_t, uint64_t, uint64_t> &overlap_tuple,
                     const uint64_t                            cluster_size) {
                    std::get<2>(overlap_tuple) = cluster_size;
                  };
              cluster_overlap_map_ptr->async_visit(
                  cluster_pair, assign_size2_lambda, cluster_size);
            };
        cluster_size_map2.async_visit(cluster_pair.second, get_size2_lambda,
                                      cluster_pair, cluster_overlap_map_ptr);
      };
  cluster_overlap_map.for_all(get_cluster_size_lambda);
  mpi_comm.barrier();

  return;
}

/**
 * Calculate entropy for a clustering from its cluster sizes
 *
 * @param cluster_size_map A YGM map of cluster id -> cluster size for a
 * clustering.
 * @param num_points The number of points
 * @return The entropy of the clustering
 *
 */
double calculate_clustering_entropy(
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map,
    uint64_t                                        num_points) {
  ygm::comm &mpi_comm = cluster_size_map.comm();

  double local_entropy  = 0;
  auto   entropy_lambda = [&local_entropy, &num_points](
                            [[maybe_unused]] const cluster_id_type &cluster_id,
                            const uint64_t &cluster_size) {
    double fraction =
        static_cast<double>(cluster_size) / static_cast<double>(num_points);
    local_entropy = local_entropy - fraction * std::log(fraction);
  };

  cluster_size_map.for_all(entropy_lambda);
  mpi_comm.barrier();

  double entropy = ygm::sum(local_entropy, mpi_comm);
  return entropy;
}

/**
 * Calculate joint entropy for two clusterings from their overlap map
 *
 * @param cluster_overlap_map YGM map of cluster pair (i,j) ->
 * tuple(overlap size, size of cluster i, size of cluster j), where i is a
 * cluster_id in the first clustering, j is a cluster_id in the second
 * clustering, overlap size is the size of the intersection between clusters i
 * and j. This is the cluster overlap map we use when calculating clustering
 * metrics with information theoretic metrics.
 * @param num_points The number of points.
 * @return The joint entropy for the pair of clusterings.
 *
 */
double calculate_joint_entropy(
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                        std::tuple<uint64_t, uint64_t, uint64_t>>
            &cluster_overlap_map,
    uint64_t num_points) {
  ygm::comm &mpi_comm = cluster_overlap_map.comm();

  double local_joint_entropy = 0;
  auto   joint_entropy_lambda =
      [&local_joint_entropy, &num_points](
          [[maybe_unused]] const std::pair<cluster_id_type, cluster_id_type>
                                                         &cluster_pair,
          const std::tuple<uint64_t, uint64_t, uint64_t> &overlap_tuple) {
        uint64_t overlap_size = std::get<0>(overlap_tuple);
        double   fraction =
            static_cast<double>(overlap_size) / static_cast<double>(num_points);
        local_joint_entropy =
            local_joint_entropy - fraction * std::log(fraction);
      };
  cluster_overlap_map.for_all(joint_entropy_lambda);
  mpi_comm.barrier();

  double joint_entropy = ygm::sum(local_joint_entropy, mpi_comm);
  return joint_entropy;
}

/**
 * Calculate joint entropy for two clusterings from their overlap map
 *
 * @param cluster_overlap_map YGM map of cluster pair (i,j) ->
 * tuple(overlap size, size of cluster i, size of cluster j), where i is a
 * cluster_id in the first clustering, j is a cluster_id in the second
 * clustering, overlap size is the size of the intersection between clusters i
 * and j. This is the cluster overlap map we use when calculating clustering
 * metrics with information theoretic metrics.
 * @return The sum of squares of the the overlap size for each pair of clusters
 * i, j, where i is in clustering 1 and j is in clustering 2
 */
uint64_t calculate_sum_squares_overlap_mi(
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                        std::tuple<uint64_t, uint64_t, uint64_t>>
        &cluster_overlap_map) {
  ygm::comm &mpi_comm = cluster_overlap_map.comm();

  // Calculate the sums of overlap sizes squared
  uint64_t local_sum_squares_overlap = 0;
  auto     sum_squares_overlap_lambda =
      [&local_sum_squares_overlap](
          [[maybe_unused]] const std::pair<cluster_id_type, cluster_id_type>
                                                         &cluster_pair,
          const std::tuple<uint64_t, uint64_t, uint64_t> &overlap_tuple) {
        uint64_t overlap_size = std::get<0>(overlap_tuple);
        local_sum_squares_overlap =
            local_sum_squares_overlap + overlap_size * overlap_size;
      };
  cluster_overlap_map.for_all(sum_squares_overlap_lambda);
  mpi_comm.barrier();

  uint64_t sum_squares_overlap = ygm::sum(local_sum_squares_overlap, mpi_comm);
  return sum_squares_overlap;
}

/**
 * Fills a cluster_size_count_map with the numbers of clusters of each size by
 * reading a cluster_size_map of cluster id to size.
 *
 * @param cluster_size_map A YGM map of cluster id -> cluster size for a
 * clustering.
 * @param cluster_size_count_map An empty YGM map of cluster size -> number of
 * clusters with that size
 */
void fill_cluster_size_count_map(
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map,
    ygm::container::map<uint64_t, uint64_t>        &cluster_size_count_map) {
  ygm::comm &mpi_comm = cluster_size_map.comm();

  cluster_size_count_map.clear();

  auto get_size_counts_lambda =
      [&cluster_size_count_map](
          [[maybe_unused]] const cluster_id_type &cluster_id,
          const uint64_t                         &cluster_size) {
        cluster_size_count_map.async_visit(
            cluster_size,
            [](const uint64_t &cluster_size, uint64_t &count) { count += 1; });
      };
  cluster_size_map.for_all(get_size_counts_lambda);
  mpi_comm.barrier();
}

/**
 * Fills an empty size_pair_map with the pairs of cluster size options in the
 * two clusterings and the number of times those size pairs occur.
 *
 * @param size_pair_map An empty YGM map of pair(size1,size2) -> pair(number of
 * instances, adjustment_contribution), where (size1, size2) are the sizes for
 * pairs of clusters in different clusterings and we always pick size1 <= size2,
 * number of instances is how many distinct pairs of clusters in different
 * clusterings have sizes (size1, size2), and adjustment_contribution is the
 * contribution to the innermost sum in the expected value of mutual
 * information. This function will fill the keys (size1, size2) of the map and
 * the instances part of the values.
 * @param cluster_size_count_map1 A YGM map for clustering 1 of cluster size ->
 * number of clusters with that size
 * @param cluster_size_count_map1 A YGM map for clustering 2 of cluster size ->
 * number of clusters with that size
 */
void get_size_pairs_and_counts_for_size_pair_map(
    ygm::container::map<std::pair<uint64_t, uint64_t>,
                        std::pair<uint64_t, double>> &size_pair_map,
    ygm::container::map<uint64_t, uint64_t>          &cluster_size_count_map1,
    ygm::container::map<uint64_t, uint64_t>          &cluster_size_count_map2) {
  ygm::comm &mpi_comm          = size_pair_map.comm();
  auto       size_pair_map_ptr = size_pair_map.get_ygm_ptr();

  size_pair_map.clear();

  // We make a local copy of which ever clustering has fewer distinct sizes
  bool clustering1_has_fewer_size_counts = true;
  if (cluster_size_count_map2.size() < cluster_size_count_map1.size()) {
    bool clustering1_has_fewer_size_counts = false;
  }

  // Vector contains pairs of cluster size, number of clusters with that
  // size
  std::vector<std::pair<uint64_t, uint64_t>> local_vector_for_this_rank;
  local_vector_for_this_rank.clear();
  auto populate_vector_for_this_rank_lambda =
      [&local_vector_for_this_rank](
          [[maybe_unused]] const uint64_t &cluster_size,
          const uint64_t                  &count) {
        local_vector_for_this_rank.push_back(
            std::make_pair(cluster_size, count));
      };
  if (clustering1_has_fewer_size_counts) {
    cluster_size_count_map1.for_all(populate_vector_for_this_rank_lambda);
  } else {
    cluster_size_count_map2.for_all(populate_vector_for_this_rank_lambda);
  }
  mpi_comm.barrier();

  // Pool all the entries on this rank together so each rank gets a
  // full copy of vector
  static std::vector<std::pair<uint64_t, uint64_t>> local_size_count_vector;

  local_size_count_vector.clear();
  mpi_comm.async_bcast(
      [](std::vector<std::pair<uint64_t, uint64_t>> local_vector) {
        local_size_count_vector.insert(local_size_count_vector.end(),
                                       local_vector.begin(),
                                       local_vector.end());
      },
      local_vector_for_this_rank);
  mpi_comm.barrier();

  /* Get the count of each cluster-size pair */
  auto get_distinct_size_pairs_lambda = [&size_pair_map_ptr](
                                            [[maybe_unused]] const uint64_t
                                                           &cluster_size,
                                            const uint64_t &count) {
    for (auto local_map_pair : local_size_count_vector) {
      std::pair<uint64_t, uint64_t> size_pair =
          std::make_pair(std::min(cluster_size, local_map_pair.first),
                         std::max(cluster_size, local_map_pair.first));

      uint64_t increase = count * local_map_pair.second;

      auto increment_size_pair_count_lambda =
          [&count, &local_map_pair](
              [[maybe_unused]] const std::pair<uint64_t, uint64_t> &size_pair,
              std::pair<uint64_t, double> &pair_info, uint64_t increase) {
            pair_info.first = pair_info.first + increase;
          };

      size_pair_map_ptr->async_visit(
          size_pair, increment_size_pair_count_lambda, increase);
    };
  };

  if (clustering1_has_fewer_size_counts) {
    cluster_size_count_map2.for_all(get_distinct_size_pairs_lambda);
  } else {
    cluster_size_count_map1.for_all(get_distinct_size_pairs_lambda);
  }
  mpi_comm.barrier();
}

/**
 * Fills the expected mutual information contribution part of a
 * size_pair_map. This must be run after
 * get_size_pairs_and_counts_for_size_pair_map.
 *
 * @param size_pair_map A YGM map of pair(size1,size2) -> pair(number of
 * instances, adjustment_contribution), where (size1, size2) are the sizes for
 * pairs of clusters in different clusterings and we always pick size1 <= size2,
 * number of instances is how many distinct pairs of clusters in different
 * clusterings have sizes (size1, size2), and adjustment_contribution is the
 * contribution to the innermost sum in the expected value of mutual
 * information. This function will fill the adjustment_contribution part of the
 * map.
 * @param num_points The number of points.
 */
void fill_size_pair_map_expected_mi_contributions(
    ygm::container::map<std::pair<uint64_t, uint64_t>,
                        std::pair<uint64_t, double>> &size_pair_map,
    uint64_t                                          num_points) {
  ygm::comm &mpi_comm          = size_pair_map.comm();
  auto       size_pair_map_ptr = size_pair_map.get_ygm_ptr();

  /*
    Make lookup vectors to speed up expected mutual information
    calculation For this rank, look at all the cluster size pairs (ai,
    bj) with ai <= bj.
    - The small lgamma vector has length max(bj)+1 and gives us
    n_ij!, (ai - n_ij)!, and (bj - n_ij)! in the AMI calculation.
    local_lookup_vector_small_lgamma[n] = std::lgamma(n + 1) = log(n!)
    - The large lgamma vector has length max(ai + bj) - min(bj) + 1 and
    gives us (N - ai - bj + n_ij)! (here, N = number of points, n_ij =
    size of overlap) local_lookup_vector_large_lgamma[n] = std::lgamma(N
    - (n + min(bj)) + 1)
    - The log vector has length max(ai)+1
    local_lookup_vector_log[n] = n * log(n) / N
  */
  uint64_t local_max_ai            = 0;
  uint64_t local_max_bj            = 0;
  uint64_t local_min_bj            = num_points;
  uint64_t local_max_sum_size_pair = 0;

  auto get_lookup_vector_lengths_lambda =
      [&local_max_ai, &local_max_bj, &local_min_bj, &local_max_sum_size_pair](
          const std::pair<uint64_t, uint64_t>          &size_pair,
          [[maybe_unused]] std::pair<uint64_t, double> &pair_info) {
        if (size_pair.first > local_max_ai) {
          local_max_ai = size_pair.first;
        }
        if (size_pair.second > local_max_bj) {
          local_max_bj = size_pair.second;
        }
        if (size_pair.second < local_min_bj) {
          local_min_bj = size_pair.second;
        }
        uint64_t size_pair_sum = size_pair.first + size_pair.second;
        if (size_pair_sum > local_max_sum_size_pair) {
          local_max_sum_size_pair = size_pair_sum;
        }
      };
  size_pair_map.for_all(get_lookup_vector_lengths_lambda);
  mpi_comm.barrier();

  // If there's no size pairs store on this rank (e.g., in small tests),
  // then we need to make local_min_bj = 0 to avoid errors
  if (local_min_bj == num_points) {
    local_min_bj = 0;
  }

  // Populate our look-up vectors
  std::vector<double> local_lookup_vector_small_lgamma(local_max_bj + 1, 1.0);
  uint64_t            size_local_lookup_vector_large_lgamma =
      std::min(local_max_sum_size_pair - local_min_bj + 1,
               num_points - local_min_bj + 1);
  std::vector<double> local_lookup_vector_large_lgamma(
      size_local_lookup_vector_large_lgamma, 1.0);
  std::vector<double> local_lookup_vector_log(local_max_ai + 1, 0.0);

  local_lookup_vector_small_lgamma[0] = 0;
  local_lookup_vector_large_lgamma[0] =
      std::lgamma(num_points - local_min_bj + 1);
  local_lookup_vector_log[0] = 0;
  if (local_max_ai >= 1) {
    for (uint64_t n = 1; n <= local_max_ai; ++n) {
      local_lookup_vector_small_lgamma[n] = std::lgamma(n + 1);
      local_lookup_vector_log[n] = static_cast<double>(n) * std::log(n) /
                                   static_cast<double>(num_points);
    }
  }
  if (local_max_ai + 1 <= local_max_bj) {
    for (uint64_t n = local_max_ai + 1; n <= local_max_bj; ++n) {
      local_lookup_vector_small_lgamma[n] = std::lgamma(n + 1);
    }
  }
  for (uint64_t n = 0; n < size_local_lookup_vector_large_lgamma; ++n) {
    local_lookup_vector_large_lgamma[n] =
        std::lgamma(num_points - (n + local_min_bj) + 1);
  }
  mpi_comm.barrier();

  // For each cluster-size pair, calculate the contribution to the
  // innermost sum of the expected value for adjusted mutual information
  double log_num_points    = std::log(num_points);
  double lgamma_num_points = std::lgamma(num_points + 1);

  auto calculate_expectation_contribution_lambda =
      [&log_num_points, &lgamma_num_points, &num_points,
       &local_lookup_vector_small_lgamma, &local_lookup_vector_large_lgamma,
       &local_min_bj, &local_lookup_vector_log,
       &size_local_lookup_vector_large_lgamma](
          [[maybe_unused]] const std::pair<uint64_t, uint64_t> &size_pair,
          std::pair<uint64_t, double>                          &pair_info) {
        double expectation_contribution = 0.0;
        double shared_term1 = 1.0 / static_cast<double>(num_points) *
                              (log_num_points - std::log(size_pair.first) -
                               std::log(size_pair.second));
        double shared_term2 = std::lgamma(size_pair.first + 1) +
                              std::lgamma(size_pair.second + 1) +
                              std::lgamma(num_points - size_pair.first + 1) +
                              std::lgamma(num_points - size_pair.second + 1) -
                              lgamma_num_points;

        uint64_t start_index =
            std::max(static_cast<int64_t>(1),
                     static_cast<int64_t>(size_pair.first + size_pair.second -
                                          num_points));
        uint64_t end_index =
            size_pair.first;  // since size_pair.first <= size_pair.second

        double term1, term2;

        for (uint64_t n_ij = start_index; n_ij <= end_index; ++n_ij) {
          term1 = local_lookup_vector_log[n_ij] +
                  static_cast<double>(n_ij) * shared_term1;
          term2 = shared_term2 - local_lookup_vector_small_lgamma[n_ij] -
                  local_lookup_vector_small_lgamma[size_pair.first - n_ij] -
                  local_lookup_vector_small_lgamma[size_pair.second - n_ij] -
                  local_lookup_vector_large_lgamma[size_pair.first +
                                                   size_pair.second - n_ij -
                                                   local_min_bj];
          term2 = std::exp(term2);
          expectation_contribution += term1 * term2;
        }

        pair_info.second = expectation_contribution;
      };
  size_pair_map.for_all(calculate_expectation_contribution_lambda);
  mpi_comm.barrier();

  local_lookup_vector_small_lgamma.clear();
  local_lookup_vector_large_lgamma.clear();
  local_lookup_vector_log.clear();

  return;
}

/**
 * Calculates the expected mutual information from size_pair_map.
 *
 * @param size_pair_map A YGM map of pair(size1,size2) -> pair(number of
 * instances, adjustment_contribution), where (size1, size2) are the sizes for
 * pairs of clusters in different clusterings and we always pick size1 <= size2,
 * number of instances is how many distinct pairs of clusters in different
 * clusterings have sizes (size1, size2), and adjustment_contribution is the
 * contribution to the innermost sum in the expected value of mutual
 * information.
 * @param num_points The number of points.
 */
double calculate_expected_mutual_information(
    ygm::container::map<std::pair<uint64_t, uint64_t>,
                        std::pair<uint64_t, double>> &size_pair_map) {
  ygm::comm &mpi_comm = size_pair_map.comm();

  // Calculate the expected mutual information
  double local_expected_mutual_info = 0.0;
  auto   expected_mutual_info_lambda =
      [&local_expected_mutual_info](
          [[maybe_unused]] const std::pair<uint64_t, uint64_t> &size_pair,
          const std::pair<uint64_t, double>                    &pair_info) {
        local_expected_mutual_info +=
            static_cast<double>(pair_info.first) * pair_info.second;
      };
  size_pair_map.for_all(expected_mutual_info_lambda);
  mpi_comm.barrier();

  double expected_mutual_info = ygm::sum(local_expected_mutual_info, mpi_comm);
  return expected_mutual_info;
}