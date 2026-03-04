#include <clustering_metrics_utils.hpp>

/* Show usage */
void show_help() {
  std::cout
      << "Compare clusterings by calculating clustering metrics without any "
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

// ------------------------------

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

    world.cout0(
        "Calculating clustering metrics without information theoretic "
        "measures");
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

    // map of cluster id in clustering 2 -> cluster size
    ygm::container::map<cluster_id_type, uint64_t> cluster_size_map2(world);

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

    /* Compare clusterings for each second clustering file */
    for (std::filesystem::path &clustering2_path :
         opt.vector_of_clustering2_paths) {
      per_clustering_timer.reset();
      step_timer.reset();

      // Make sure the cluster size maps are cleared
      cluster_size_map1.clear();
      cluster_size_map2.clear();
      world.barrier();

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

      /*
      Create a YGM map of cluster pair (i,j) -> overlap size
      where i is a cluster_id in the first clustering, j is a cluster_id
      in the second clustering, and overlap size is the size of the intersection
      between clusters i and j
      */
      ygm::container::map<std::pair<cluster_id_type, cluster_id_type>, uint64_t>
          cluster_overlap_map(world);

      fill_cluster_overlap_and_size_maps(point_to_clusters_map,
                                         cluster_overlap_map, cluster_size_map1,
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

      /* Calculate the sums of squares of cluster sizes */

      uint64_t sum_squares_cluster1, sum_squares_cluster2, sum_squares_overlap;
      sum_squares_cluster1 =
          calculate_sums_of_squares_for_map_values(cluster_size_map1);
      sum_squares_cluster2 =
          calculate_sums_of_squares_for_map_values(cluster_size_map2);
      sum_squares_overlap =
          calculate_sums_of_squares_for_map_values(cluster_overlap_map);

      if (opt.verbose) {
        world.cout0("\nSum of cluster sizes squared for clustering 1: ",
                    sum_squares_cluster1);
        world.cout0("Sum of cluster sizes squared for clustering 2: ",
                    sum_squares_cluster2);
        world.cout0("Sum of overlaps squared for the two clusterings: ",
                    sum_squares_overlap);
        world.cout0("Time to calculate sums of squares: ", step_timer.elapsed(),
                    " seconds");
      }

      step_timer.reset();
      world.barrier();

      // Clear cluster size maps before calculating purity
      cluster_size_map1.clear();
      cluster_size_map2.clear();
      world.barrier();

      /* Calculate purity */
      double purity = 0.0;
      if (opt.calculate_purity) {
        purity = calculate_purity(cluster_overlap_map, num_points);
      }

      /* Calculate various clustering metrics */
      if (world.rank() == 0) {
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
        std::cout << "\nPair-confusion balanced accuracy: "
                  << clustering_metrics_no_mi.balanced_accuracy << std::endl;
        std::cout << "Pair-confusion geometric mean: "
                  << clustering_metrics_no_mi.geometric_mean << std::endl;

        if (opt.verbose) {
          // Calculate Rand Index (not adjusted)
          double numerator = static_cast<double>(
              num_points * (num_points - 1) + 2 * sum_squares_overlap -
              sum_squares_cluster1 - sum_squares_cluster2);
          double denominator =
              static_cast<double>(num_points * (num_points - 1));

          double rand_index = numerator / denominator;

          std::cout << "\nUnadjusted Rand Index: " << rand_index << std::endl;
        }
      }

      if (opt.verbose) {
        world.cout0("\nTime to calculate clustering metrics: ",
                    step_timer.elapsed(), " seconds");
      }
      world.cout0("Total time to calculate metrics for this clustering: ",
                  per_clustering_timer.elapsed(), " seconds");

    }  // For each clustering_file2
  }

  return 0;
}
