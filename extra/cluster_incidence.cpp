
#include <clustering_metrics_utils.hpp>
#include <ygm/io/multi_output.hpp>

using point_id_type   = uint64_t;
using cluster_id_type = int32_t;

/* Show usage */
void show_help() {
  std::cout
      << "Outputs an incidence file where each row is \\t separated with "
      << "<cluster id in clustering2> <cluster id in clustering2> "
      << "<the number of cluster ids in clustering 1 present in both clusters>"
      << std::endl
      << std::endl
      << "<<Usage>>" << std::endl
      << "Usage: cluster_incidence -g clustering1_path -c clustering2_path -o "
         "output_file"
      << std::endl
      << std::endl
      << "Required arguments: " << std::endl
      << "  -g <path> Path to file(s) for clustering 1. This is the clustering "
         "that we use to determine the overlap between clusters, e.g., the "
         "'ground truth' clustering. Can provide a single file or a directory "
         "containing all files for the clustering."
      << std::endl
      << " -c <path> Path to file(s) for clustering 2. This is the clustering "
         "we "
         "will determine the incidence of. Can provide a single "
         "file or a directory containing all files for the clustering."
      << std::endl
      << "  -o <path> Output file." << std::endl
      << std::endl
      << "Optional arguments:" << std::endl
      << "  -v Verbose printout." << std::endl
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

    world.cout0("Calculating cluster incidence");
    world.cout0("The number of ranks is ", world.size(), "\n");
    ygm::utility::timer step_timer{};
    ygm::utility::timer per_clustering_timer{};

    /* Parse the user inputs into files we will read */

    // Getopt example:
    // https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html

    /* Read the file names for comparison*/
    std::filesystem::path clustering1_path;
    std::filesystem::path clustering2_path;
    std::filesystem::path output_path;

    static bool verbose_printout;
    verbose_printout = false;

    // Read the optional arguments if available
    int opt_char;
    while ((opt_char = getopt(argc, argv, "vhg:c:o:")) != -1) {
      switch (opt_char) {
        case 'g':
          clustering1_path = optarg;
          break;
        case 'c':
          clustering2_path = optarg;
          break;
        case 'o':
          output_path = optarg;
          break;
        case 'v':
          verbose_printout = true;
          break;
        case 'h':
          if (world.rank() == 0) {
            show_help();
          }
          return 1;
        default:
          if (world.rank() == 0) {
            std::cerr << "Unrecognized option: " << opt_char << ", ignore."
                      << std::endl;
            show_help();
          }
      }
    }

    // Throw an error and quit if no clustering 1 file is provided
    if (std::filesystem::is_empty(clustering1_path) ||
        !std::filesystem::exists(clustering1_path)) {
      if (world.rank() == 0) {
        std::cout << "Error: No valid clustering 1 file provided";
        show_help();
      }
      return 1;
    }

    // Throw an error and quit if no clustering 2 file is provided
    if (std::filesystem::is_empty(clustering2_path) ||
        !std::filesystem::exists(clustering2_path)) {
      if (world.rank() == 0) {
        std::cout << "Error: No valid clustering 2 file provided";
        show_help();
      }
      return 1;
    }

    /* Read the clustering 1 file(s) */

    // map of cluster id in clustering 1 -> vector of points in that cluster
    ygm::container::map<cluster_id_type, std::vector<point_id_type>>
        cluster_to_members_map1(world);

    {
      step_timer.reset();
      world.cout0() << "Reading provided clustering 1 file/directory: "
                    << clustering1_path << std::endl;

      std::vector<std::filesystem::path> paths =
          find_file_paths(clustering1_path);

      static uint64_t local_num_noise_points_clustering1;
      local_num_noise_points_clustering1 = 0;
      static uint64_t local_num_points_clustering1;
      local_num_points_clustering1 = 0;

      // YGM line parser
      std::vector<std::string> file_vector;
      for (const std::filesystem::path &p : paths) {
        file_vector.push_back(p.c_str());
      }
      ygm::io::line_parser line_parser_clustering1(world, file_vector);

      // Read each line in the input file
      auto read_input_line_lambda =
          [&cluster_to_members_map1](const std::string &line) {
            if (std::isdigit(line[0])) {
              try {
                point_id_type     point_id;
                cluster_id_type   cluster_id;
                std::stringstream ss(line);

                ss >> point_id >> cluster_id;

                // Only consider clusters that are not labeled as noise
                if (is_cluster_noise(cluster_id)) {
                  ++local_num_noise_points_clustering1;
                } else {
                  ++local_num_points_clustering1;
                  cluster_to_members_map1.async_visit(
                      cluster_id,
                      [](const cluster_id_type      &cluster_id,
                         std::vector<point_id_type> &members,
                         const point_id_type        &new_member) {
                        members.push_back(new_member);
                      },
                      point_id);
                }
              } catch (...) {
                std::cout << "Error reading line: " << line << std::endl;
              }
            } else {
              std::cout << "Read comment line: " << line << std::endl;
            }
          };
      line_parser_clustering1.for_all(read_input_line_lambda);

      world.barrier();
      world.cout0() << "Number of clusters in clustering 1: "
                    << cluster_to_members_map1.size() << std::endl;
      world.cout0()
          << "Number of points with valid cluster ids in clustering 1: "
          << ygm::sum(local_num_points_clustering1, world) << std::endl;
      uint64_t num_noise_points_clustering1 =
          ygm::sum(local_num_noise_points_clustering1, world);
      if (num_noise_points_clustering1 > 0) {
        world.cout0() << "Number of noise points in clustering 1: "
                      << num_noise_points_clustering1 << std::endl;
      }
      world.cout0() << "Time to read first clustering file (s): "
                    << step_timer.elapsed() << std::endl
                    << std::endl;
    }

    /* Read clustering 2 file(s) */

    // map of point id -> cluster id in clustering 2
    ygm::container::map<point_id_type, cluster_id_type> point_to_cluster_map2(
        world);

    {
      step_timer.reset();
      world.cout0() << "Reading provided clustering 2 file/directory: "
                    << clustering2_path << std::endl;

      std::vector<std::filesystem::path> paths =
          find_file_paths(clustering2_path);

      static uint64_t local_num_noise_points_clustering2;
      local_num_noise_points_clustering2 = 0;

      // YGM line parser
      std::vector<std::string> file_vector;
      for (const std::filesystem::path &p : paths) {
        file_vector.push_back(p.c_str());
      }
      ygm::io::line_parser line_parser_clustering2(world, file_vector);

      // Read each line in the input file
      auto read_input_line_lambda =
          [&point_to_cluster_map2](const std::string &line) {
            if (std::isdigit(line[0])) {
              try {
                point_id_type     point_id;
                cluster_id_type   cluster_id;
                std::stringstream ss(line);

                ss >> point_id >> cluster_id;

                // Only consider clusters that are not labeled as noise
                if (is_cluster_noise(cluster_id)) {
                  ++local_num_noise_points_clustering2;
                } else {
                  point_to_cluster_map2.async_insert(point_id, cluster_id);
                }
              } catch (...) {
                std::cout << "Error reading line: " << line << std::endl;
              }
            } else {
              std::cout << "Read comment line: " << line << std::endl;
            }
          };
      line_parser_clustering2.for_all(read_input_line_lambda);

      world.barrier();
      world.cout0()
          << "Number of points with valid cluster ids in clustering 2: "
          << point_to_cluster_map2.size() << std::endl;
      uint64_t num_noise_points_clustering2 =
          ygm::sum(local_num_noise_points_clustering2, world);
      if (num_noise_points_clustering2 > 0) {
        world.cout0() << "Number of noise points in clustering 2: "
                      << num_noise_points_clustering2 << std::endl;
      }
      world.cout0() << "Time to read second clustering file (s): "
                    << step_timer.elapsed() << std::endl
                    << std::endl;
    }

    /*
    Use the previous maps point_to_cluster_map2 and cluster_to_members_map1
    to get a map of
    cluster i in clustering 1
      --> set of (clusters in clustering 2 represented in cluster i)
    */

    ygm::container::map<cluster_id_type, std::set<cluster_id_type>>
         cluster_to_represented_clusters_map(world);
    auto cluster_to_represented_clusters_map_ptr =
        cluster_to_represented_clusters_map.get_ygm_ptr();

    {
      step_timer.reset();
      world.cout0() << "Getting map of clustering 1 clusters to clustering 2 "
                       "representation counts"
                    << std::endl;

      auto process_clusters_lambda =
          [&cluster_to_represented_clusters_map_ptr, &point_to_cluster_map2](
              [[maybe_unused]] const cluster_id_type &clustering1_cluster_id,
              std::vector<point_id_type>             &members) {
            auto add_cluster_rep_lambda =
                []([[maybe_unused]] const point_id_type &point,
                   const cluster_id_type                &clustering2_cluster_id,
                   const cluster_id_type                &clustering1_cluster_id,
                   auto cluster_to_represented_clusters_map_ptr) {
                  cluster_to_represented_clusters_map_ptr->async_visit(
                      clustering1_cluster_id,
                      []([[maybe_unused]] const cluster_id_type
                             &clustering1_cluster_id,
                         std::set<cluster_id_type>
                             &represented_clustering2_clusters,
                         const cluster_id_type &clustering2_cluster_id) {
                        represented_clustering2_clusters.insert(
                            clustering2_cluster_id);
                      },
                      clustering2_cluster_id);
                };
            for (const point_id_type &member : members) {
              point_to_cluster_map2.async_visit(
                  member, add_cluster_rep_lambda, clustering1_cluster_id,
                  cluster_to_represented_clusters_map_ptr);
            }
          };
      cluster_to_members_map1.for_all(process_clusters_lambda);
      world.barrier();

      world.cout0() << "Time elapsed (s): " << step_timer.elapsed()
                    << std::endl;
    }

    /* Calculate cluster incidence*/

    // Map of (cluster i in clustering2, cluster j in clustering 2) -> number of
    // clusters in clustering 1 represented in both i and j
    // Only store pairs i <= j to avoid double counting
    ygm::container::map<std::pair<cluster_id_type, cluster_id_type>, uint32_t>
        incidence_map(world);

    {
      step_timer.reset();
      world.cout0() << "Getting clustering 2 cluster incidence" << std::endl;

      auto get_incidence_lambda =
          [&incidence_map](
              [[maybe_unused]] const cluster_id_type &clustering1_cluster_id,
              const std::set<cluster_id_type>
                  &represented_clustering2_clusters) {
            for (const cluster_id_type &i : represented_clustering2_clusters) {
              for (const cluster_id_type &j :
                   represented_clustering2_clusters) {
                if (i <= j) {
                  std::pair<cluster_id_type, cluster_id_type> cluster_pair = {
                      i, j};
                  incidence_map.async_visit(
                      cluster_pair,
                      [](const std::pair<cluster_id_type, cluster_id_type>
                                  &cluster_pair,
                         uint32_t &incidence) { ++incidence; });
                }
              }
            }
          };
      cluster_to_represented_clusters_map.for_all(get_incidence_lambda);
      world.barrier();

      world.cout0() << "Time elapsed (s): " << step_timer.elapsed()
                    << std::endl;
    }

    /* Write to file */

    std::string        output_dir = output_path.parent_path().c_str();
    static std::string filename;
    filename = output_path.filename().c_str();

    // YGM multi output for writing to file, append = false
    ygm::io::multi_output mo(world, output_dir, 1024, false);

    {
      step_timer.reset();
      world.cout0() << "Writing incidence file" << std::endl;

      auto write_incidence_lambda =
          [&mo](const std::pair<cluster_id_type, cluster_id_type> &cluster_pair,
                const uint32_t                                    &incidence) {
            std::stringstream ss;
            ss << cluster_pair.first << "\t" << cluster_pair.second << "\t"
               << incidence;
            mo.async_write_line(filename, ss.str());
          };
      incidence_map.for_all(write_incidence_lambda);
      world.barrier();

      world.cout0() << "Time elapsed (s): " << step_timer.elapsed()
                    << std::endl;
    }
  }
  return 0;
}