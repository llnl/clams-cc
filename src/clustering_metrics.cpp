#include <clustering_metrics_mi.hpp>
#include <clustering_metrics_utils.hpp>

using point_id_type   = uint64_t;
using cluster_id_type = int32_t;

/* Show usage */
void show_help() {
  std::cout
      << "Compare clusterings by calculating clustering metrics including "
         "information theoretic quantities."
      << std::endl
      << std::endl
      << "<<Usage>>" << std::endl
      << "Usage: clustering_metrics -g ground_truth_file [optional arguments] "
         "individual_clustering_file"
      << std::endl
      << std::endl
      << "Required arguments: " << std::endl
      << "  -g <path> Path to file(s) for clustering 1. This is the clustering "
         "that we compare everything else to, e.g., the `ground truth' "
         "clustering. Can provide a single file or a directory containing all "
         "files for the clustering."
      << std::endl
      << std::endl
      << "Optional arguments: " << std::endl
      << "  -l <path> Path to file containing a list of clusterings to compare "
         "to. In this case, we compare clustering 1 to each clustering in this "
         "list."
      << std::endl
      << "  -p Calculate and output purity. This assumes that clustering 1 is "
         "the ground truth and calculates the purity of clustering 2 with "
         "respect to clustering 1."
      << std::endl
      << "  -v Verbose. Print out more quantities leading up to calculation of "
         "clustering metrics."
      << std::endl
      << "  -h Show help." << std::endl
      << std::endl
      << "Note: Clustering files should be in the form `point_id cluster_id' "
         "separated by a white space on each line."
      << std::endl;
}

// --------------------------------

int main(int argc, char **argv) {
  ygm::comm world(&argc, &argv);
  {
    /*
    Print check for defining what cluster id's that represent noise look like.
    */

    if (std::is_signed<cluster_id_type>::value) {
      world.cout0(
          "Using signed cluster id type. Points are considered noise "
          "if cluster_id < 0");
      world.cout0("is_cluster_noise(-1) = ", is_cluster_noise(-1), "\n");
    } else {
      world.cout0(
          "Using unsigned cluster id type. Points are considered noise "
          "if cluster_id == max value for type");
      world.cout0(
          "is_cluster_noise(", std::numeric_limits<cluster_id_type>::max(),
          ") = ", is_cluster_noise(std::numeric_limits<cluster_id_type>::max()),
          "\n");
    }

    world.cout0("Calculating clustering metrics");
    world.cout0("The number of ranks is ", world.size(), "\n");
    ygm::utility::timer step_timer{};
    ygm::utility::timer per_clustering_timer{};

    /* Parse the user inputs into files we will read */
    clustering_metrics_option_t opt;
    bool                        unsuccessful_parse =
        parse_clustering_metrics_options(argc, argv, opt, world);
    if (unsuccessful_parse) {
      show_help();
      return 1;
    }

    /*
    Initial set up for ygm maps
    */

    // map of point id -> pair of cluster id in clustering 1, cluster id in
    // clustering 2
    ygm::container::map<point_id_type,
                        std::pair<cluster_id_type, cluster_id_type>>
        point_to_clusters_map(world);

    // map of cluster id in clustering 1 -> cluster size
    ygm::container::map<cluster_id_type, uint64_t> cluster_size_map1(world);

    // map of cluster size -> number of clusters in clustering 1 with that
    // cluster size
    ygm::container::map<uint64_t, uint64_t> cluster_size_count_map1(world);

    // map of cluster id in clustering 2 -> cluster size
    ygm::container::map<cluster_id_type, uint64_t> cluster_size_map2(world);

    // map of cluster size -> number of clusters in clustering 2 with that
    // cluster size
    ygm::container::map<uint64_t, uint64_t> cluster_size_count_map2(world);

    /*
    Read the first clustering file (the one we compare all the other clusterings
    to)
    */

    // number of noise points read in the clustering 1 file
    uint64_t num_noise_points_clustering1;

    // Read and process the first clustering file (see
    // clustering_metric_utils.hpp)
    num_noise_points_clustering1 = read_first_clustering_file(
        opt.clustering1_path, point_to_clusters_map, cluster_size_map1);

    world.barrier();
    world.cout0("Time to read first clustering file: ", step_timer.elapsed(),
                " seconds\n");

    uint64_t num_points_clustering1 = point_to_clusters_map.size();
    world.cout0("Number of points with valid cluster ids in clustering 1: ",
                num_points_clustering1);
    if (num_noise_points_clustering1 > 0) {
      world.cout0("Number of noise points in clustering 1: ",
                  num_noise_points_clustering1);
    }
    world.cout0("Number of clusters in clustering 1: ",
                cluster_size_map1.size());
    step_timer.reset();

    // // Later, we'll pool all the map entries on this rank together so each
    // // rank gets a full copy of the local map
    // // Vector contains pairs of cluster size, number of clusters with that
    // size static std::vector<std::pair<uint64_t, uint64_t>>
    // local_size_count_vector;

    /* Compare clusterings for each second clustering file */
    for (std::filesystem::path &clustering2_path :
         opt.vector_of_clustering2_paths) {
      per_clustering_timer.reset();
      step_timer.reset();

      world.cout0(
          "\n-------------------------------------------------------\n");

      /* Read the second clustering file */

      // Read clustering file 2 and fill points_to_cluster_map,
      // get info on numbers of points read
      auto [num_points, num_noise_points, num_points_in_file,
            num_noise_points_in_file] =
          read_second_clustering_file(clustering2_path, point_to_clusters_map);

      world.cout0("Time to read comparison clustering file: ",
                  step_timer.elapsed(), " seconds\n");

      world.cout0(
          "Number of points with valid cluster ids in clustering 2 that "
          "are also in clustering 1: ",
          num_points);
      world.cout0(
          "Number of noise points in clustering 2 (that are also points "
          "in clustering 1): ",
          num_noise_points);
      double percent_clustered =
          100.0 *
          static_cast<double>(num_points_in_file - num_noise_points_in_file) /
          static_cast<double>(num_points_in_file);
      world.cout0("Clustering 2 - Clustered %: ", percent_clustered);
      double percent_clustered_shared_with_clustering1 =
          100.0 * static_cast<double>(num_points) /
          static_cast<double>(num_points + num_noise_points);
      if (percent_clustered_shared_with_clustering1 != percent_clustered) {
        world.cout0(
            "Clustering 2 - Clustered % when only considering points that "
            "exist in clustering 1: ",
            percent_clustered_shared_with_clustering1);
      }

      if (opt.verbose) {
        world.cout0("\nNumber of points read from clustering 2 file: ",
                    num_points_in_file);
        world.cout0("Number of noise points read from clustering 2 file: ",
                    num_noise_points_in_file);
        world.cout0(
            "Number points that are noise in both clusterings 1 and 2: ",
            num_noise_points_in_file - num_noise_points);
      }

      world.barrier();
      step_timer.reset();

      /* Create a YGM map of hash(i,j) -> tuple(overlap size, size of cluster
      i, size of cluster j), where i is a cluster_id in the first clustering,
      j is a cluster_id in the second clustering, and overlap size is the size
      of the intersection between clusters i and j
      */
      ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                          std::tuple<uint64_t, uint64_t, uint64_t>>
          cluster_overlap_map(world);

      // Process point_to_clusters_map to fill cluster_overlap_map,
      // cluster_size_map1, and cluster_size_map2
      fill_cluster_overlap_and_size_maps_mi(
          point_to_clusters_map, cluster_overlap_map, cluster_size_map1,
          cluster_size_map2);

      world.cout0("\nNumber of overlapping cluster pairs: ",
                  cluster_overlap_map.size());
      world.cout0("Number of clusters in clustering 1: ",
                  cluster_size_map1.size());
      world.cout0("Number of clusters in clustering 2: ",
                  cluster_size_map2.size());
      if (opt.verbose) {
        world.cout0(
            "Time to create YGM maps of cluster overlaps and cluster sizes: ",
            step_timer.elapsed(), " seconds");
      }

      step_timer.reset();
      world.barrier();

      /* Calculate the entropies required to calculate mutual information */

      // Calculate the entropy for clustering 1
      double entropy1 =
          calculate_clustering_entropy(cluster_size_map1, num_points);

      // Calculate the entropy of clustering 2
      double entropy2 =
          calculate_clustering_entropy(cluster_size_map2, num_points);

      // Calculate the joint entropy
      double joint_entropy =
          calculate_joint_entropy(cluster_overlap_map, num_points);

      if (opt.verbose) {
        world.cout0("\nEntropy for clustering 1: ", entropy1);
        world.cout0("Entropy for clustering 2: ", entropy2);
        world.cout0("Joint entropy: ", joint_entropy);
      }

      /* Calculate the sums of sizes squared required for various clustering
       * metrics */

      uint64_t sum_squares_cluster1, sum_squares_cluster2, sum_squares_overlap;
      sum_squares_cluster1 =
          calculate_sums_of_squares_for_map_values(cluster_size_map1);
      sum_squares_cluster2 =
          calculate_sums_of_squares_for_map_values(cluster_size_map2);

      // Calculate the sums of overlap sizes squared
      sum_squares_overlap =
          calculate_sum_squares_overlap_mi(cluster_overlap_map);

      if (opt.verbose) {
        world.cout0("Sum of cluster sizes squared for clustering 1: ",
                    sum_squares_cluster1);
        world.cout0("Sum of cluster sizes squared for clustering 2: ",
                    sum_squares_cluster2);
        world.cout0("Sum of overlaps squared for the two clusterings: ",
                    sum_squares_overlap);
        world.cout0("Time to calculate entropies and sums of sizes squared: ",
                    step_timer.elapsed(), " seconds\n");
      }

      step_timer.reset();
      world.barrier();

      /*
      Get cluster size counts for calculations for the expected value of
      mutual information
      */

      cluster_size_count_map1.clear();
      cluster_size_count_map2.clear();
      world.barrier();

      // For both clusterings 1 and 2, get map of cluster size -> number of
      // clusters in clustering with that size
      fill_cluster_size_count_map(cluster_size_map1, cluster_size_count_map1);
      fill_cluster_size_count_map(cluster_size_map2, cluster_size_count_map2);

      if (opt.verbose) {
        world.cout0("Number of distinct cluster sizes for clustering 1: ",
                    cluster_size_count_map1.size());
        world.cout0("Number of distinct cluster sizes for clustering 2: ",
                    cluster_size_count_map2.size());
      }

      // Clear cluster size maps since we're done with them
      cluster_size_map1.clear();
      cluster_size_map2.clear();
      world.barrier();

      /* Calculate purity */
      double purity = 0.0;
      if (opt.calculate_purity) {
        purity = calculate_purity(cluster_overlap_map, num_points);
      }

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

      // Get the size pairs and how often they occur
      get_size_pairs_and_counts_for_size_pair_map(
          size_pair_map, cluster_size_count_map1, cluster_size_count_map2);

      world.cout0(
          "Number of distinct cluster size "
          "pairs: ",
          size_pair_map.size());
      if (opt.verbose) {
        world.cout0("Time to get the distinct cluster size pairs: ",
                    step_timer.elapsed(), " seconds");
      }
      step_timer.reset();
      world.barrier();

      // Clear some space by freeing up size count maps
      cluster_size_count_map1.clear();
      cluster_size_count_map2.clear();

      // Get the expected mutual information adjustment contribution for each
      // size pair
      fill_size_pair_map_expected_mi_contributions(size_pair_map, num_points);

      // Calculate the expected mutual information
      double expected_mutual_info =
          calculate_expected_mutual_information(size_pair_map);
      world.cout0("Expected mutual information: ", expected_mutual_info);
      if (opt.verbose) {
        world.cout0("Time to calculate expected mutual information: ",
                    step_timer.elapsed(), " seconds");
      }
      step_timer.reset();
      world.barrier();

      /* Calculate various clustering metrics */
      if (world.rank() == 0) {
        // Calculate normalized mutual information
        double mutual_info           = entropy1 + entropy2 - joint_entropy;
        double geometric_normalizer  = std::sqrt(entropy1 * entropy2);
        double arithmetic_normalizer = 0.5 * (entropy1 + entropy2);

        if (opt.verbose) {
          std::cout << "\nUnnormalized mutual information I(X,Y) = "
                    << mutual_info << std::endl;
          std::cout << "Geometric normalized mutual information "
                    << "I(X,Y) / sqrt(H(X)H(Y)) = "
                    << mutual_info / geometric_normalizer << std::endl;
          std::cout << "Arithmetic normalized mutual information "
                    << "I(X,Y) / 0.5(H(X) + H(Y)) = "
                    << mutual_info / arithmetic_normalizer << std::endl;
        }

        // Calculate adjusted mutual information
        double arithmetic_adjusted_mutual_info =
            (mutual_info - expected_mutual_info) /
            (arithmetic_normalizer - expected_mutual_info);
        double geometric_adjusted_mutual_info =
            (mutual_info - expected_mutual_info) /
            (geometric_normalizer - expected_mutual_info);
        std::cout
            << "\n(Arithmetic-normalized) Adjusted Mutual Information (AMI): "
            //   << "{I(X,Y) - E[I(X,Y)]} / {0.5(H(X) + H(Y)) - E[I(X,Y)]} = "
            << arithmetic_adjusted_mutual_info << std::endl;
        if (opt.verbose) {
          std::cout
              << "Geometric-normalized adjusted mutual information: "
              //   << "{I(X,Y) - E[I(X,Y)]} / {sqrt(H(X)H(Y)) - E[I(X,Y)]} = "
              << geometric_adjusted_mutual_info << std::endl;
        }

        // Calculate clustering metrics that only need sums of squares
        clustering_metrics_no_mi_type clustering_metrics_no_mi;
        clustering_metrics_no_mi = calculate_clustering_metrics_no_mi(
            num_points, sum_squares_overlap, sum_squares_cluster1,
            sum_squares_cluster2);

        std::cout << "\nAdjusted Rand Index (ARI): "
                  << clustering_metrics_no_mi.adjusted_rand_index << std::endl;
        std::cout << "Fowlkes Mallows Index: "
                  << clustering_metrics_no_mi.fowlkes_mallows << std::endl;
        if (opt.calculate_purity) {
          std::cout << "(Assumes clustering 1 is the ground truth) Purity: "
                    << purity << std::endl;
        }
        std::cout << "\nPair-confusion Balanced Accuracy: "
                  << clustering_metrics_no_mi.balanced_accuracy << std::endl;
        std::cout << "Pair-confusion Geometric Mean: "
                  << clustering_metrics_no_mi.geometric_mean << std::endl;

        if (opt.verbose) {
          // Calculate Rand Index (without adjustment)
          double numerator = static_cast<double>(
              num_points * (num_points - 1) + 2 * sum_squares_overlap -
              sum_squares_cluster1 - sum_squares_cluster2);
          double denominator =
              static_cast<double>(num_points * (num_points - 1));

          double rand_index = numerator / denominator;

          std::cout << "\nUnadjusted Rand Index: " << rand_index << std::endl;
        }
      }

      world.cout0("\nTime to calculate clustering metrics: ",
                  step_timer.elapsed(), " seconds");
      world.cout0("Total time to calculate metrics for this clustering: ",
                  per_clustering_timer.elapsed(), " seconds");

    }  // For each clustering_file2
  }

  return 0;
}
