#pragma once

#include <ygm/comm.hpp>
#include <ygm/container/map.hpp>
#include <ygm/detail/collective.hpp>
#include <ygm/io/line_parser.hpp>
#include <ygm/utility/timer.hpp>

#include <common_utils.hpp>

using point_id_type   = uint64_t;
using cluster_id_type = int32_t;

/**
 * @brief
 * Define what cluster id's that represent noise look like.
 */
// If cluster ids are signed, we say that negative cluster ids (usually -1)
// are noise points
template <typename T>
typename std::enable_if<std::is_signed<T>::value, bool>::type is_cluster_noise(
    const T &cluster_id) {
  return (cluster_id < 0) ? true : false;
}
// If the cluster ids are unsigned, we say that the max possible value of
// the cluster id type indicates that a point is noise
template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, bool>::type
is_cluster_noise(const T &cluster_id) {
  T max_value = std::numeric_limits<T>::max();
  return (cluster_id == max_value) ? true : false;
}

/**
 * @brief
 * Holds user options for running clustering metrics programs.
 */
struct clustering_metrics_option_t {
  std::filesystem::path              clustering1_path;
  std::vector<std::filesystem::path> vector_of_clustering2_paths;
  bool                               calculate_purity = false;
  bool                               verbose          = false;
};

/**
 * @brief
 * Parse user options for computing clustering metrics using getopt
 */
bool parse_clustering_metrics_options(int argc, char *argv[],
                                      clustering_metrics_option_t &opt,
                                      ygm::comm                   &world) {
  std::filesystem::path list_of_clustering2_files;

  // Read the optional arguments if available
  int opt_char;
  while ((opt_char = getopt(argc, argv, "pvhg:s:l:")) != -1) {
    switch (opt_char) {
      case 'g':
        opt.clustering1_path = optarg;
        break;
      case 'l':
        list_of_clustering2_files = optarg;
        break;
      case 'p':
        opt.calculate_purity = true;
        break;
      case 'v':
        opt.verbose = true;
        break;
      case 'h':
        return 1;
      default:
        if (world.rank() == 0) {
          std::cerr << "Unrecognized option: " << opt_char << ", ignore."
                    << std::endl;
        }
    }
  }

  // If we have any non-optional arguments, we'll assume they are more
  // clusterings to compare with clustering 1
  for (int index = optind; index < argc; index++) {
    if (std::filesystem::exists(argv[index])) {
      std::filesystem::path p = argv[index];
      opt.vector_of_clustering2_paths.push_back(p);
    } else {
      if (world.rank() == 0) {
        std::cout << "Invalid extra non-optional argument will be ignored: "
                  << argv[index] << std::endl;
      }
    }
  }

  // Throw an error and quit if no clustering 1 file is provided
  if (!std::filesystem::exists(opt.clustering1_path) ||
      std::filesystem::is_empty(opt.clustering1_path)) {
    if (world.rank() == 0) {
      std::cout << "Error: No valid clustering 1 file provided";
    }
    return 1;
  }

  // If we have file that contains a list (one file per line) of comparision
  // clusterings read in all the file names in that list
  if (std::filesystem::exists(list_of_clustering2_files) &&
      !std::filesystem::is_empty(list_of_clustering2_files)) {
    // Read the files in the list of clustering files into the vector
    world.cout0("Reading comparison files from list in: ",
                list_of_clustering2_files);
    std::ifstream clustering_ifs(list_of_clustering2_files);
    std::string   line;
    while (std::getline(clustering_ifs, line)) {
      opt.vector_of_clustering2_paths.push_back(line);
    }
    clustering_ifs.close();
    clustering_ifs.clear();
  }

  // If we weren't provided with a comparison file, compare clustering 1 to
  // itself
  if (opt.vector_of_clustering2_paths.size() == 0) {
    opt.vector_of_clustering2_paths.push_back(opt.clustering1_path);
  }

  return 0;
}

/**
 * Read and process the first clustering file (the one we compare all the other
 * clusterings to).
 *
 * @param clustering1_path String with path to first clustering file(s). This
 * can be a single file or a directory of multiple files.
 * @param point_to_clusters_map An empty YGM map for point_id -> (cluster_id in
 * clustering 1, cluster_id in clustering 2).
 * @param cluster_size_map1 An empty YGM map for cluster_id -> cluster size in
 * clustering 1.
 * @return The number of noise points read in the first clustering file.
 */

template <typename point_id_type, typename cluster_id_type>
uint64_t read_first_clustering_file(
    std::filesystem::path &clustering1_path,
    ygm::container::map<point_id_type,
                        std::pair<cluster_id_type, cluster_id_type>>
                                                   &point_to_clusters_map,
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map1) {
  ygm::comm &mpi_comm = point_to_clusters_map.comm();

  mpi_comm.cout0() << "Reading provided clustering 1 file/directory: "
                   << clustering1_path << std::endl;

  std::vector<std::filesystem::path> paths = find_file_paths(clustering1_path);

  static uint64_t local_num_noise_points_clustering1;
  local_num_noise_points_clustering1 = 0;

  // YGM line parser
  std::vector<std::string> file_vector;
  for (const std::filesystem::path &p : paths) {
    file_vector.push_back(p.c_str());
  }
  ygm::io::line_parser line_parser_clustering1(mpi_comm, file_vector);

  auto increment_cluster_size_lambda =
      []([[maybe_unused]] const cluster_id_type &cluster_id,
         uint64_t &cluster_size) { cluster_size += 1; };

  // Read each line in the input file
  auto read_input_line_lambda1 =
      [&point_to_clusters_map, &cluster_size_map1,
       &increment_cluster_size_lambda](const std::string &line) {
        if (std::isdigit(line[0])) {
          try {
            point_id_type                               point_id;
            cluster_id_type                             cluster_id;
            std::pair<cluster_id_type, cluster_id_type> cluster_pair(0, 0);
            std::stringstream                           ss(line);

            ss >> point_id >> cluster_id;

            // Only consider clusters that are not labeled as noise
            if (is_cluster_noise(cluster_id)) {
              ++local_num_noise_points_clustering1;
            } else {
              cluster_pair.first = cluster_id;
              point_to_clusters_map.async_insert(point_id, cluster_pair);
              cluster_size_map1.async_visit(cluster_id,
                                            increment_cluster_size_lambda);
            }
          } catch (...) {
            std::cout << "Error reading line: " << line << std::endl;
          }
        } else {
          std::cout << "Read comment line: " << line << std::endl;
        }
      };
  line_parser_clustering1.for_all(read_input_line_lambda1);

  mpi_comm.barrier();
  uint64_t num_noise_points_clustering1 =
      ygm::sum(local_num_noise_points_clustering1, mpi_comm);

  return num_noise_points_clustering1;
}

/**
 * Read and process the comparison clustering file.
 *
 * @param clustering2_path String with path to first clustering file(s). This
 * can be a single file or a directory of multiple files.
 * @param point_to_clusters_map An empty YGM map for point_id -> (cluster_id in
 * clustering 1, cluster_id in clustering 2).
 * @return Pair of (number of valid points, number of noise points) read in
 * comparision clustering.
 */
template <typename point_id_type, typename cluster_id_type>
std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> read_second_clustering_file(
    std::filesystem::path &clustering2_path,
    ygm::container::map<point_id_type,
                        std::pair<cluster_id_type, cluster_id_type>>
        &point_to_clusters_map) {
  ygm::comm &mpi_comm = point_to_clusters_map.comm();

  mpi_comm.cout0() << "Reading provided clustering 2 file/directory: "
                   << clustering2_path << std::endl;

  std::vector<std::filesystem::path> paths = find_file_paths(clustering2_path);

  // Number of points and number of noise points in clustering 2
  // that are also in clustering 1
  static uint64_t local_num_points;
  static uint64_t local_num_noise_points;
  local_num_points       = 0;
  local_num_noise_points = 0;

  // Total number of points and number of noise points read from file
  // (need not be in clustering 1)
  static uint64_t local_num_points_in_file;
  static uint64_t local_num_noise_points_in_file;
  local_num_points_in_file       = 0;
  local_num_noise_points_in_file = 0;

  std::vector<std::string> file_vector{};
  for (const std::filesystem::path &p : paths) {
    file_vector.push_back(p.c_str());
  }
  ygm::io::line_parser line_parser_clustering2(mpi_comm, file_vector);

  auto read_input_line_lambda2 =
      [&point_to_clusters_map](const std::string &line) {
        if (std::isdigit(line[0])) {
          try {
            point_id_type     point_id;
            cluster_id_type   cluster_id;
            std::stringstream ss(line);
            ss >> point_id >> cluster_id;

            // If we read a numeric point in the clustering 2 file,
            // track the total even if the point isn't in clustering 1
            ++local_num_points_in_file;

            // Add points in clustering 2 that exist in clustering 1 to the
            // point_to_clusters_map

            // If the point is a noise point in clustering 2, also note that
            if (is_cluster_noise(cluster_id)) {
              ++local_num_noise_points_in_file;

              auto add_second_cluster_thats_noise_lambda =
                  []([[maybe_unused]] const point_id_type        &point_id,
                     std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
                     const cluster_id_type &new_cluster_id) {
                    cluster_pair.second = new_cluster_id;
                    ++local_num_noise_points;
                  };
              point_to_clusters_map.async_visit_if_contains(
                  point_id, add_second_cluster_thats_noise_lambda, cluster_id);

            } else {
              auto add_second_cluster_lambda =
                  []([[maybe_unused]] const point_id_type        &point_id,
                     std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
                     const cluster_id_type &new_cluster_id) {
                    cluster_pair.second = new_cluster_id;
                    ++local_num_points;
                  };
              point_to_clusters_map.async_visit_if_contains(
                  point_id, add_second_cluster_lambda, cluster_id);
            }
          } catch (...) {
            std::cout << "Error reading line: " << line << std::endl;
          }
        } else {
          std::cout << "Read comment line: " << line << std::endl;
        }
      };

  line_parser_clustering2.for_all(read_input_line_lambda2);

  mpi_comm.barrier();

  // number of points read in clustering 2 whose id exists in clustering 1
  uint64_t num_points = ygm::sum(local_num_points, mpi_comm);
  // number of noise points read in clustering 2 whose id exists in clustering 1
  uint64_t num_noise_points = ygm::sum(local_num_noise_points, mpi_comm);
  // total number of points read from file for clustering 2
  uint64_t num_points_in_file = ygm::sum(local_num_points_in_file, mpi_comm);
  // total number of noise points read from file for clustering 2
  uint64_t num_noise_points_in_file =
      ygm::sum(local_num_noise_points_in_file, mpi_comm);
  mpi_comm.barrier();

  return std::make_tuple(num_points, num_noise_points, num_points_in_file,
                         num_noise_points_in_file);
}

/**
 * Process point_to_clusters_map to fill cluster_overlap_map for calculating
 * non-information theoretic clustering metrics
 *
 * @param point_to_clusters_map YGM map for point_id -> (cluster_id in
 * clustering 1, cluster_id in clustering 2).
 * @param cluster_overlap_map An empty YGM map of cluster pair (i,j) -> overlap
 * size where i is a cluster_id in the first clustering, j is a cluster_id in
 * the second clustering, overlap size is the size of the intersection between
 * clusters i and j. This is the cluster overlap map we use when calculating
 * clustering metrics without information theoretic metrics.
 * @param cluster_size_map1 An empty YGM map for cluster_id -> cluster size in
 * clustering 1.
 * @param cluster_size_map2 An empty YGM map for cluster_id -> cluster size in
 * clustering 2.
 *
 */
void fill_cluster_overlap_and_size_maps(
    ygm::container::map<point_id_type,
                        std::pair<cluster_id_type, cluster_id_type>>
        &point_to_clusters_map,
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>, uint64_t>
                                                   &cluster_overlap_map,
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map1,
    ygm::container::map<cluster_id_type, uint64_t> &cluster_size_map2) {
  ygm::comm &mpi_comm = point_to_clusters_map.comm();

  // Make sure maps are empty to start
  cluster_overlap_map.clear();
  cluster_size_map1.clear();
  cluster_size_map2.clear();

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
               uint64_t &overlap) { overlap = overlap + 1; };

        // Lambda to increase the cluster size in cluster_size_map 1 & 2
        auto increment_cluster_size_lambda =
            []([[maybe_unused]] const cluster_id_type &cluster_id,
               uint64_t &cluster_size) { cluster_size += 1; };

        // Only consider points that are not noise
        if (!is_cluster_noise(cluster_pair.second)) {
          cluster_overlap_map.async_visit(cluster_pair, update_overlap_lambda);
          cluster_size_map1.async_visit(cluster_pair.first,
                                        increment_cluster_size_lambda);
          cluster_size_map2.async_visit(cluster_pair.second,
                                        increment_cluster_size_lambda);
        }
      };

  point_to_clusters_map.for_all(process_point_lambda);
  mpi_comm.barrier();

  return;
}

/**
 * Calculate the sums of values squared in a ygm map of key -> value.
 *
 * @param ygm_map YGM map where the values are squarable.
 * @return Sum of squares of the values in the YGM map.
 */
template <typename map_key_type>
uint64_t calculate_sums_of_squares_for_map_values(
    ygm::container::map<map_key_type, uint64_t> &ygm_map) {
  ygm::comm &mpi_comm = ygm_map.comm();

  uint64_t local_sum_squares  = 0;
  auto     sum_squares_lambda = [&local_sum_squares](
                                [[maybe_unused]] const map_key_type &key,
                                const uint64_t                      &value) {
    local_sum_squares = local_sum_squares + value * value;
  };

  ygm_map.for_all(sum_squares_lambda);
  mpi_comm.barrier();

  uint64_t sum_squares = ygm::sum(local_sum_squares, mpi_comm);

  return sum_squares;
}

/**
 * @brief holds the non-information-theoretic clustering metrics from comparing
 * two clusterings
 */
struct clustering_metrics_no_mi_type {
  /** @brief Adjusted Rand index */
  double adjusted_rand_index;
  /** @brief Fowlkes--Mallows index */
  double fowlkes_mallows;
  /** @brief The balanced accuracy is the arithmetic mean of sensitivity and
   * specificity */
  double balanced_accuracy;
  /** @brief Geometric mean of sensitivity and specificity */
  double geometric_mean;
};

/**
 * Calculate non-information-theoretic clustering metrics.
 *
 * @param num_points The number of points.
 * @param sum_squares_overlap The sum of the squares of the sizes of the cluster
 * overlaps in the two clusterings.
 * @param sum_squares_cluster1 The sum of the squares of the sizes of the
 * clusters in clustering 1.
 * @param sum_squares_cluster2 The sum of the squares of the sizes of the
 * clusters in clustering 2.
 * @return Struct containing the calculated clustering metrics. The values are:
 *         adjusted_rand_index, fowlkes_mallows, balanced_accuracy,
 * geometric_mean.
 */
clustering_metrics_no_mi_type calculate_clustering_metrics_no_mi(
    uint64_t num_points, uint64_t sum_squares_overlap,
    uint64_t sum_squares_cluster1, uint64_t sum_squares_cluster2) {
  clustering_metrics_no_mi_type clustering_metrics_no_mi = {};

  // Calculate Adjusted Rand Index
  double temp = static_cast<double>((sum_squares_cluster1 - num_points) *
                                    (sum_squares_cluster2 - num_points)) /
                static_cast<double>(num_points * (num_points - 1));
  double numerator =
      static_cast<double>(sum_squares_overlap - num_points) - temp;
  double denominator =
      0.5 * static_cast<double>(sum_squares_cluster1 + sum_squares_cluster2 -
                                2 * num_points) -
      temp;

  clustering_metrics_no_mi.adjusted_rand_index = numerator / denominator;

  // Calculate Fowlkes Mallows Index
  double fowlkes_mallows =
      static_cast<double>(sum_squares_overlap - num_points) /
      static_cast<double>(sum_squares_cluster1 - num_points);
  fowlkes_mallows = fowlkes_mallows *
                    static_cast<double>(sum_squares_overlap - num_points) /
                    static_cast<double>(sum_squares_cluster2 - num_points);
  fowlkes_mallows = std::sqrt(fowlkes_mallows);

  clustering_metrics_no_mi.fowlkes_mallows = fowlkes_mallows;

  // Balanced accuracy and geometric mean version
  // Here, we assume that clustering 1 is the ``ground truth'' and
  // that the positive class is pairs of items that are clustered together in
  // clustering 1

  // Sensitivity = TP / P = TP / (TP + FN)
  numerator          = static_cast<double>(sum_squares_overlap - num_points);
  denominator        = static_cast<double>(sum_squares_cluster1 - num_points);
  double sensitivity = numerator / denominator;

  // Specificity = TN / N = TN / (TN + FP)
  numerator =
      static_cast<double>(num_points * num_points + sum_squares_overlap -
                          sum_squares_cluster1 - sum_squares_cluster2);
  denominator =
      static_cast<double>(num_points * num_points - sum_squares_cluster1);
  double specificity = numerator / denominator;

  clustering_metrics_no_mi.balanced_accuracy =
      0.5 * (sensitivity + specificity);
  clustering_metrics_no_mi.geometric_mean =
      std::sqrt(sensitivity * specificity);

  return clustering_metrics_no_mi;
}

/**
 * Processes a YGM cluster_overlap_map to fill a max_ground_truth_overlap_map.
 * Auxiliary function for calculating purity.
 *
 * @param cluster_overlap_map A YGM map of cluster pair (i,j) -> overlap size
 * where i is a cluster_id in the first clustering, j is a cluster_id in the
 * second clustering, overlap size is the size of the intersection between
 * clusters i and j. This is the cluster overlap map we use when calculating
 * clustering metrics without information theoretic metrics.
 * @param max_ground_truth_overlap_map An empty YGM map of cluster_id in
 * clustering 2 -> (ground truth cluster with largest overlap, overlap size)
 */
void fill_max_ground_truth_overlap_map(
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>, uint64_t>
        &cluster_overlap_map,
    ygm::container::map<cluster_id_type, std::pair<cluster_id_type, uint64_t>>
        &max_ground_truth_overlap_map) {
  ygm::comm &mpi_comm = cluster_overlap_map.comm();

  max_ground_truth_overlap_map.clear();
  auto max_ground_truth_overlap_map_ptr =
      max_ground_truth_overlap_map.get_ygm_ptr();

  auto ground_truth_overlap_lambda =
      [&max_ground_truth_overlap_map_ptr](
          const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
          const uint64_t                                    &overlap_size) {
        cluster_id_type ground_truth_id = cluster_pair.first;
        cluster_id_type cluster_id      = cluster_pair.second;

        auto visit_max_overlap_lambda =
            []([[maybe_unused]] const cluster_id_type cluster_id,
               std::pair<cluster_id_type, uint64_t>  &max_ground_truth_overlap,
               const cluster_id_type                 &ground_truth_id,
               const uint64_t                        &overlap_size) {
              if (overlap_size > max_ground_truth_overlap.second) {
                max_ground_truth_overlap.first  = ground_truth_id;
                max_ground_truth_overlap.second = overlap_size;
              }
            };

        max_ground_truth_overlap_map_ptr->async_visit(
            cluster_id, visit_max_overlap_lambda, ground_truth_id,
            overlap_size);
      };
  cluster_overlap_map.for_all(ground_truth_overlap_lambda);
  mpi_comm.barrier();

  return;
}

/**
 * Processes a YGM cluster_overlap_map to fill a max_ground_truth_overlap_map.
 * Auxiliary function for calculating purity.
 *
 * @param cluster_overlap_map A YGM map of cluster pair (i,j) -> tuple(overlap
 * size, size of cluster i, size of cluster j), where i is a cluster_id in the
 * first clustering, j is a cluster_id in the second clustering, overlap size is
 * the size of the intersection between clusters i and j. This is the cluster
 * overlap map we use when calculating clustering metrics with information
 * theoretic metrics.
 * @param max_ground_truth_overlap_map An empty YGM map of cluster_id in
 * clustering 2 -> (ground truth cluster with largest overlap, overlap size)
 */
void fill_max_ground_truth_overlap_map(
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                        std::tuple<uint64_t, uint64_t, uint64_t>>
        &cluster_overlap_map,
    ygm::container::map<cluster_id_type, std::pair<cluster_id_type, uint64_t>>
        &max_ground_truth_overlap_map) {
  ygm::comm &mpi_comm = cluster_overlap_map.comm();

  max_ground_truth_overlap_map.clear();
  auto max_ground_truth_overlap_map_ptr =
      max_ground_truth_overlap_map.get_ygm_ptr();

  auto ground_truth_overlap_lambda =
      [&max_ground_truth_overlap_map_ptr](
          const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
          const std::tuple<uint64_t, uint64_t, uint64_t>    &overlap_tuple) {
        cluster_id_type ground_truth_id = cluster_pair.first;
        cluster_id_type cluster_id      = cluster_pair.second;
        uint64_t        overlap_size    = std::get<0>(overlap_tuple);

        auto visit_max_overlap_lambda =
            []([[maybe_unused]] const cluster_id_type cluster_id,
               std::pair<cluster_id_type, uint64_t>  &max_ground_truth_overlap,
               const cluster_id_type                  ground_truth_id,
               const uint64_t                        &overlap_size) {
              if (overlap_size > max_ground_truth_overlap.second) {
                max_ground_truth_overlap.first  = ground_truth_id;
                max_ground_truth_overlap.second = overlap_size;
              }
            };

        max_ground_truth_overlap_map_ptr->async_visit(
            cluster_id, visit_max_overlap_lambda, ground_truth_id,
            overlap_size);
      };
  cluster_overlap_map.for_all(ground_truth_overlap_lambda);
  mpi_comm.barrier();

  return;
}

/**
 * Calculate purity
 *
 * @param cluster_overlap_map A ygm map of cluster pair (i,j) -> overlap size
 * where i is a cluster_id in the first clustering, j is a cluster_id in the
 * second clustering, overlap size is the size of the intersection between
 * clusters i and j
 * @param num_points The number of points.
 * @return The clustering purity where clustering 1 is taken as the `ground
 * truth`. This measures the extent to which clusters in clustering 2 contain a
 * single label from clustering 1. We take the sum over all clusters in
 * clustering 2 of the size of the majority clustering 1 label in that cluster,
 * and then divide by the number of points.
 */
template <typename overlap_info_type>
double calculate_purity(
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>,
                        overlap_info_type> &cluster_overlap_map,
    uint64_t                                num_points) {
  ygm::comm &mpi_comm = cluster_overlap_map.comm();

  // Map of cluster id -> (ground truth cluster with largest overlap,
  // overlap size)
  ygm::container::map<cluster_id_type, std::pair<cluster_id_type, uint64_t>>
      max_ground_truth_overlap_map(mpi_comm);

  fill_max_ground_truth_overlap_map(cluster_overlap_map,
                                    max_ground_truth_overlap_map);

  /* Calculate purity */
  uint64_t local_purity_sum = 0;
  auto     purity_sum_lambda =
      [&local_purity_sum]([[maybe_unused]] const cluster_id_type &cluster_id,
                          const std::pair<cluster_id_type, uint64_t>
                              &max_ground_truth_overlap) {
        local_purity_sum = local_purity_sum + max_ground_truth_overlap.second;
      };
  max_ground_truth_overlap_map.for_all(purity_sum_lambda);
  mpi_comm.barrier();
  uint64_t purity_sum = ygm::sum(local_purity_sum, mpi_comm);
  double   purity =
      static_cast<double>(purity_sum) / static_cast<double>(num_points);

  return purity;
}