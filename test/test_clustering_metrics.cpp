// Test clustering metrics calculation including information-theoretic
// quantities

#undef NDEBUG

#include <clustering_metrics_mi.hpp>
#include <clustering_metrics_utils.hpp>

int main(int argc, char **argv) {
  ygm::comm world(&argc, &argv);

  // map of point id -> pair of cluster id in clustering 1, cluster id in
  // clustering 2
  ygm::container::map<point_id_type,
                      std::pair<cluster_id_type, cluster_id_type>>
      point_to_clusters_map(world);
  point_to_clusters_map.async_insert(0, std::make_pair(0, 1));
  point_to_clusters_map.async_insert(1, std::make_pair(0, 0));
  point_to_clusters_map.async_insert(2, std::make_pair(0, 0));
  point_to_clusters_map.async_insert(3, std::make_pair(1, 0));
  point_to_clusters_map.async_insert(4, std::make_pair(1, 1));
  point_to_clusters_map.async_insert(5, std::make_pair(2, 2));
  world.barrier();

  uint64_t num_points = point_to_clusters_map.size();

  /*
    Required YGM maps for calculation
  */

  // maps of cluster id in clustering -> cluster size
  ygm::container::map<cluster_id_type, uint64_t> cluster_size_map1(world);
  ygm::container::map<cluster_id_type, uint64_t> cluster_size_map2(world);

  // maps of cluster size -> number of clusters in clustering with that
  // cluster size
  ygm::container::map<uint64_t, uint64_t> cluster_size_count_map1(world);
  ygm::container::map<uint64_t, uint64_t> cluster_size_count_map2(world);

  // map of (i,j) -> tuple(overlap size, size of cluster i, size of cluster
  // j),where i is a cluster_id in the first clustering, j is a cluster_id in
  // the second clustering, and overlap size is the size of the intersection
  // between clusters i and j
  ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                      std::tuple<uint64_t, uint64_t, uint64_t>>
      cluster_overlap_map(world);

  // Process point_to_clusters_map to fill cluster_overlap_map,
  // cluster_size_map1, and cluster_size_map2
  fill_cluster_overlap_and_size_maps_mi(point_to_clusters_map,
                                        cluster_overlap_map, cluster_size_map1,
                                        cluster_size_map2);

  // Overlap options are (0,0), (0,1), (1,0), (1,1), (2,2)
  assert(cluster_overlap_map.size() == 5);

  /*
  Calculate entropies
  */
  double entropy1 = calculate_clustering_entropy(cluster_size_map1, num_points);
  double entropy2 = calculate_clustering_entropy(cluster_size_map2, num_points);
  double joint_entropy =
      calculate_joint_entropy(cluster_overlap_map, num_points);

  // std::cout << std::setprecision(10) << entropy1 << std::endl;
  // std::cout << std::setprecision(10) << entropy2 << std::endl;
  // std::cout << std::setprecision(10) << joint_entropy << std::endl;

  assert(std::abs(entropy1 - 1.0114042) < 1e-6);
  assert(std::abs(entropy2 - 1.0114042) < 1e-6);
  assert(std::abs(joint_entropy - 1.5607104) < 1e-6);

  /*
  Calculate the sums of squares required
  */

  uint64_t sum_squares_cluster1, sum_squares_cluster2, sum_squares_overlap;
  sum_squares_cluster1 =
      calculate_sums_of_squares_for_map_values(cluster_size_map1);
  sum_squares_cluster2 =
      calculate_sums_of_squares_for_map_values(cluster_size_map2);
  sum_squares_overlap = calculate_sum_squares_overlap_mi(cluster_overlap_map);

  assert(sum_squares_cluster1 == 14);
  assert(sum_squares_cluster2 == 14);
  assert(sum_squares_overlap == 8);

  // For both clusterings 1 and 2, get map of cluster size -> number of
  // clusters in clustering with that size
  fill_cluster_size_count_map(cluster_size_map1, cluster_size_count_map1);
  fill_cluster_size_count_map(cluster_size_map2, cluster_size_count_map2);

  /*
  Create a YGM map of pair(size1,size2) -> pair(number of instances,
  adjustment_contribution), where (size1, size2) are the sizes for pairs
  of clusters in different clusterings and we always pick size1 <= size2,
  number of instances is how many distinct pairs of clusters in different
  clusterings have sizes (size1, size2), and adjustment_contribution is the
  contribution to the innermost sum in the expected value of mutual
  information
  */
  ygm::container::map<std::pair<uint64_t, uint64_t>,
                      std::pair<uint64_t, double>>
      size_pair_map(world);

  get_size_pairs_and_counts_for_size_pair_map(
      size_pair_map, cluster_size_count_map1, cluster_size_count_map2);
  fill_size_pair_map_expected_mi_contributions(size_pair_map, num_points);

  // Calculate the expected mutual information
  double expected_mutual_info =
      calculate_expected_mutual_information(size_pair_map);

  assert(std::abs(expected_mutual_info - 0.4703093) < 1e-6);

  // Calculate clustering metrics that only need sums of squares
  clustering_metrics_no_mi_type clustering_metrics_no_mi;
  clustering_metrics_no_mi = calculate_clustering_metrics_no_mi(
      num_points, sum_squares_overlap, sum_squares_cluster1,
      sum_squares_cluster2);

  assert(std::abs(clustering_metrics_no_mi.adjusted_rand_index - -0.0227273) <
         1e-6);
  assert(std::abs(clustering_metrics_no_mi.fowlkes_mallows - 0.25) < 1e-6);

  // Calculate purity
  double purity = calculate_purity(cluster_overlap_map, num_points);

  assert(std::abs(purity - 0.6666667) < 1e-6);

  return 0;
}