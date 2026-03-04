#include <boost/functional/hash.hpp>
#include <ygm/comm.hpp>
#include <ygm/container/disjoint_set.hpp>
#include <ygm/container/map.hpp>
#include <ygm/io/line_parser.hpp>
#include <ygm/io/multi_output.hpp>
#include <ygm/utility/timer.hpp>

#include <common_utils.hpp>

/* Show usage */
template <typename cout_type>
void usage(std::string_view exe_name, cout_type &cout) {
  cout
      << "Outputs the counts of the number of inter- or intra- cluster "
         "edges for each combination of cluster ids. Each file line contains"
      << "<cluster id i><cluster id j><number edges between clusters i "
         "and j>"
      << std::endl
      << std::endl
      << "<<Usage>>" << std::endl
      << "Required arguments: " << std::endl
      << "  -p <path>   Path to clustering label file(s). Each file line "
         "should be of the form <point_id><cluster_id> separated by a white "
         "space. Can provide a single file or a directory containing all "
         "files for the clustering."
      << std::endl
      << "  -g <path>   Path to white-space separated edge list file(s) for "
         "the graph. Can provide a single file or a directory containing all "
         "files for the graph."
      << std::endl
      << "  -o <path>   Path to directory to write output file(s)." << std::endl
      << "Optional arguments:" << std::endl
      << "  -r <string> Path to file of nodes to ignore. For example, if "
         "we want to only consider the largest connected component (LCC), the "
         "file should have nodes not in the LCC. Each line in the file should "
         "contain a node. There can be a whitespace and other things following."
      << std::endl
      << "  -e          Write edge probability file. Outputs a file where "
         "each line is <cluster id><cluster id><fraction of possible edges "
         "between clusters>."
      << std::endl
      << "  -w          Write neighborhood purity to file. The neighborhood "
         "purity for a node is the fraction of its neighbors that belong to "
         "the majority cluster of its neighbors. Outputs a file where each "
         "line is <node><neighborhood purity>."
      << std::endl
      << "  -c          Write cluster size file. Outputs a file where each "
         "line is <cluster_id><size>."
      << std::endl
      << "  -d          Write degree file. Outputs a file where each line "
         "is <node><degree>."
      << std::endl
      << "  -s          Sort outputs when writing to file." << std::endl
      << std::endl;
}

using id_t         = uint32_t;
using cluster_id_t = uint32_t;

struct option_t {
  std::filesystem::path clustering_path;
  std::filesystem::path edge_list_path;
  std::string           output_dir;
  bool                  write_neighborhood_purity_file = false;
  bool                  write_edge_prob_file           = false;
  bool                  write_cluster_size_file        = false;
  bool                  write_degree_file              = false;
  bool                  sort_outputs                   = false;
  bool                  verbose                        = false;
  std::filesystem::path file_of_nodes_to_ignore;
};

bool parse_options(int, char **, option_t &, bool &);
template <typename cout_type>
void usage(std::string_view, cout_type &);
template <typename cout_type>
void show_options(const option_t &, cout_type &);

template <typename key_type, typename value_type>
void consolidate_key_value_vector(
    ygm::container::map<key_type, value_type> &ygm_map, ygm::comm &comm,
    std::vector<std::pair<key_type, value_type>> &consolidated_vector) {
  static std::vector<std::pair<key_type, value_type>> local_vector;
  static std::vector<std::pair<key_type, value_type>> global_vector;
  auto consolidate_lambda = [](const key_type &key, auto &value) {
    local_vector.push_back(std::make_pair(key, value));
  };
  ygm_map.for_all(consolidate_lambda);
  comm.barrier();

  // Broadcast the vector to all ranks
  comm.async_bcast(
      [](const auto &local_vector) {
        global_vector.insert(global_vector.end(), local_vector.begin(),
                             local_vector.end());
      },
      local_vector);
  comm.barrier();

  consolidated_vector.clear();
  std::copy(global_vector.begin(), global_vector.end(),
            std::back_inserter(consolidated_vector));
}

// --------------------------------

int main(int argc, char **argv) {
  ygm::comm world(&argc, &argv);
  {
    // Parse arguments
    option_t opt;
    {
      bool help        = false;
      bool parser_bool = parse_options(argc, argv, opt, help);
      if (help) {
        usage(argv[0], world.cout0());
        return 0;
      }
      show_options(opt, world.cout0());
      world.cout0() << "Number of ranks: " << world.size() << std::endl;
    }

    // timer
    ygm::utility::timer step_timer{};

    /*
      Read clustering file and get the labels for each point
    */
    ygm::container::map<id_t, cluster_id_t> point_cluster_map(world);
    {
      step_timer.reset();

      world.cout0() << "\n<<Read clustering label file>>" << std::endl;
      world.cout0() << "Provided clustering path: " << opt.clustering_path
                    << std::endl;

      std::vector<std::filesystem::path> paths =
          find_file_paths(opt.clustering_path);
      std::vector<std::string> file_vector;
      for (const std::filesystem::path &p : paths) {
        file_vector.push_back(p.c_str());
      }
      ygm::io::line_parser line_parser_clustering(world, file_vector);

      auto line_parser_lambda = [&point_cluster_map](const std::string &line) {
        if (std::isdigit(line[0])) {
          try {
            std::stringstream ss(line);
            id_t              point_id;
            cluster_id_t      cluster_id;
            ss >> point_id >> cluster_id;
            point_cluster_map.async_insert(point_id, cluster_id);
          } catch (...) {
            std::cout << "Error reading line in clustering file: " << line
                      << std::endl;
          }
        } else {
          std::cout << "Read comment line in clustering file: " << line
                    << std::endl;
        }
      };
      line_parser_clustering.for_all(line_parser_lambda);
      world.barrier();

      world.cout0() << "Number of points: " << point_cluster_map.size()
                    << std::endl;
      world.cout0() << "Time to read clustering file (s): "
                    << step_timer.elapsed() << std::endl;
    }

    /*
      Read file of nodes to ignore if provided
    */
    static bool ignore_some_nodes;
    ignore_some_nodes = false;
    static std::vector<id_t> nodes_to_ignore;
    {
      if (std::filesystem::exists(opt.file_of_nodes_to_ignore)) {
        world.cout0() << "\n<<Read nodes to ignore>>" << std::endl;
        step_timer.reset();

        static std::vector<id_t> local_nodes_to_ignore;
        if (world.rank() == 0) {
          std::ifstream ifs(opt.file_of_nodes_to_ignore);
          if (!ifs) {
            std::cerr << "Failed to open the file: "
                      << opt.file_of_nodes_to_ignore << std::endl;
          }
          std::string line;
          id_t        node;
          while (std::getline(ifs, line)) {
            std::istringstream iss(line);
            iss >> node;
            local_nodes_to_ignore.push_back(node);
          }
        }

        // Broadcast the nodes to ignore vector to all ranks
        world.async_bcast(
            [](const std::vector<id_t> &local_nodes_to_ignore) {
              nodes_to_ignore.insert(nodes_to_ignore.end(),
                                     local_nodes_to_ignore.begin(),
                                     local_nodes_to_ignore.end());
            },
            local_nodes_to_ignore);
      }
      world.barrier();

      if (nodes_to_ignore.size() > 0) {
        world.cout0() << "Number of nodes to ignore: " << nodes_to_ignore.size()
                      << std::endl;
        world.cout0() << "Time to read file: " << step_timer.elapsed()
                      << std::endl;
      }
    }

    // YGM adjacency map: node -> set of neighbors
    ygm::container::map<id_t, std::set<id_t>> adjacency_map(world);
    static auto                               insert_adjacent_node_lambda =
        []([[maybe_unused]] const id_t &id, auto &neighbor_set,
           auto &new_entry) { neighbor_set.insert(new_entry); };

    /*
       Read graph into adjacency map
    */
    {
      step_timer.reset();
      world.cout0() << "\n<<Read graph>>" << std::endl;
      world.cout0() << "Provided edge_list path: " << opt.edge_list_path
                    << std::endl;

      std::vector<std::filesystem::path> paths =
          find_file_paths(opt.edge_list_path);
      std::vector<std::string> file_vector;
      for (const std::filesystem::path &p : paths) {
        file_vector.push_back(p.c_str());
      }
      ygm::io::line_parser edge_list_line_parser(world, file_vector);

      // Lambda to edge list file and put neighbors in adjacency map
      auto read_edge_list_lambda = [&adjacency_map](const std::string &line) {
        std::stringstream ss(line);
        id_t              node1;
        id_t              node2;
        ss >> node1 >> node2;
        if (ignore_some_nodes) {
          if (std::find(nodes_to_ignore.begin(), nodes_to_ignore.end(),
                        node1) == nodes_to_ignore.end() &&
              std::find(nodes_to_ignore.begin(), nodes_to_ignore.end(),
                        node2) == nodes_to_ignore.end()) {
            adjacency_map.async_visit(node1, insert_adjacent_node_lambda,
                                      node2);
            adjacency_map.async_visit(node2, insert_adjacent_node_lambda,
                                      node1);
          }
        } else {
          adjacency_map.async_visit(node1, insert_adjacent_node_lambda, node2);
          adjacency_map.async_visit(node2, insert_adjacent_node_lambda, node1);
        }
      };
      edge_list_line_parser.for_all(read_edge_list_lambda);
      world.barrier();

      // Calculate the mean and max degree
      static size_t sum_num_neighbors;
      static size_t local_max_degree;
      sum_num_neighbors             = 0;
      local_max_degree              = 0;
      auto add_num_neighbors_lambda = []([[maybe_unused]] const id_t &id,
                                         const std::set<id_t> &neighbors) {
        sum_num_neighbors += neighbors.size();
        if (neighbors.size() > local_max_degree) {
          local_max_degree = neighbors.size();
        }
      };
      adjacency_map.for_all(add_num_neighbors_lambda);
      world.barrier();

      world.cout0() << "Number of nodes: " << adjacency_map.size() << std::endl;
      world.cout0() << "Average degree: "
                    << static_cast<double>(ygm::sum(sum_num_neighbors, world)) /
                           static_cast<double>(adjacency_map.size())
                    << std::endl;
      world.cout0() << "Max degree: " << ygm::max(local_max_degree, world)
                    << std::endl;
      world.cout0() << "Time to read graph (s): " << step_timer.elapsed()
                    << std::endl;
    }

    /* Get connected components */

    // Create a disjoint set
    ygm::container::disjoint_set<id_t> connected_components(world);

    // Map of component representative -> component size
    ygm::container::map<id_t, size_t> component_size_map(world);

    // Pair of (component representative, component size)
    static std::pair<id_t, size_t> largest_cc_pair;
    largest_cc_pair = std::make_pair(0, 0);

    // Number of connected components
    static size_t num_components;
    num_components = 0;

    {
      step_timer.reset();
      world.cout0() << "\n<<Get the connected components using a disjoint set>>"
                    << std::endl;

      // Add edges to disjoint set
      auto add_edges_lambda = [&connected_components](
                                  const id_t           &node,
                                  const std::set<id_t> &neighbors) {
        for (const id_t &neighbor : neighbors) {
          connected_components.async_union(node, neighbor);
        }
      };
      adjacency_map.for_all(add_edges_lambda);
      world.barrier();
      world.cout0() << "All edges added to disjoint set" << std::endl;

      connected_components.all_compress();
      num_components = connected_components.num_sets();
      world.cout0() << "Number of CCs:\t" << num_components << std::endl;

      // Get the size of the connected components
      auto increase_component_size_lambda =
          [&component_size_map]([[maybe_unused]] const id_t &node,
                                const id_t                  &representative) {
            component_size_map.async_visit(
                representative, []([[maybe_unused]] const id_t &key,
                                   size_t &value) { value = value + 1; });
          };

      connected_components.for_all(increase_component_size_lambda);
      world.barrier();

      // Find largest connected component and print out component sizes
      static std::pair<id_t, size_t> local_largest_cc_pair;
      local_largest_cc_pair = std::make_pair(0, 0);

      auto update_local_largest_cc_lambda = [](const id_t   &representative,
                                               const size_t &size) {
        if (size > local_largest_cc_pair.second) {
          local_largest_cc_pair.first  = representative;
          local_largest_cc_pair.second = size;
        }
        std::cout << "Component with representative " << representative
                  << " has size: " << size << std::endl;
      };
      component_size_map.for_all(update_local_largest_cc_lambda);

      // Go through all ranks to get the global largest connected
      // component
      // Pool all the local info to find the potential seed with
      // min local conductance with smallest local conductance
      world.async_bcast(
          [](const std::pair<id_t, size_t> &local_largest_cc_pair) {
            if (local_largest_cc_pair.second > largest_cc_pair.second) {
              largest_cc_pair.first  = local_largest_cc_pair.first;
              largest_cc_pair.second = local_largest_cc_pair.second;
            }
          },
          local_largest_cc_pair);
      world.barrier();

      world.cout0() << "Largest CC size: " << largest_cc_pair.second
                    << std::endl;
      world.cout0() << "Largest CC ID: " << largest_cc_pair.first << std::endl;
      world.cout0() << "Time to get connected components: \t"
                    << step_timer.elapsed() << std::endl;
    }

    // YGM map of node -> (its own cluster id, vector of cluster ids of
    // neighbors)
    ygm::container::map<id_t,
                        std::pair<cluster_id_t, std::vector<cluster_id_t>>>
                cluster_adjacency_map(world);
    auto        cluster_adjacency_map_ptr = cluster_adjacency_map.get_ygm_ptr();
    static auto update_own_cluster_id_lambda =
        []([[maybe_unused]] const id_t &node, auto &cluster_info_pair,
           auto &cluster_id) { cluster_info_pair.first = cluster_id; };
    static auto insert_neighbor_cluster_id_lambda =
        []([[maybe_unused]] const id_t &node, auto &cluster_info_pair,
           auto &neighbor_cluster_id) {
          cluster_info_pair.second.push_back(neighbor_cluster_id);
        };

    /* Process the graph so we link each node to its cluster id and get the
     * cluster ids of each node's set of neighboring nodes*/
    {
      world.cout0() << "\n<<Process graph to associate nodes with clusters>>"
                    << std::endl;
      step_timer.reset();

      // For each node in the graph, update its vector of cluster ids of
      // neighbors
      auto process_graph_lambda = [&cluster_adjacency_map_ptr,
                                   &point_cluster_map](
                                      const id_t           &node,
                                      const std::set<id_t> &neighbors) {
        auto visit_point_cluster_map_lambda =
            []([[maybe_unused]] const id_t &neighbor,
               const cluster_id_t &neighbor_cluster_id, const id_t &node,
               auto &cluster_adjacency_map_ptr) {
              cluster_adjacency_map_ptr->async_visit(
                  node, insert_neighbor_cluster_id_lambda, neighbor_cluster_id);
            };

        for (const auto &neighbor : neighbors) {
          point_cluster_map.async_visit(neighbor,
                                        visit_point_cluster_map_lambda, node,
                                        cluster_adjacency_map_ptr);
        }
      };
      adjacency_map.for_all(process_graph_lambda);
      world.barrier();

      adjacency_map.clear();

      // For each node in the graph, update its own cluster id
      auto get_own_node_cluster_ids_lambda =
          [&point_cluster_map, &cluster_adjacency_map_ptr](
              const id_t &node, [[maybe_unused]] auto &cluster_info_pair) {
            auto visit_point_cluster_map_lambda =
                [](const id_t &node, const cluster_id_t &cluster_id,
                   auto &cluster_adjacency_map_ptr) {
                  cluster_adjacency_map_ptr->async_visit(
                      node, update_own_cluster_id_lambda, cluster_id);
                };
            point_cluster_map.async_visit(node, visit_point_cluster_map_lambda,
                                          cluster_adjacency_map_ptr);
          };
      cluster_adjacency_map.for_all(get_own_node_cluster_ids_lambda);
      world.barrier();

      world.cout0() << "cluster_adjacency_map.size() = "
                    << cluster_adjacency_map.size() << std::endl;
      world.cout0() << "Time to go through graph and link cluster ids (s):"
                    << step_timer.elapsed() << std::endl;
    }

    // Map of number of inter-cluster edges:
    // (cluster 1, cluster 2) -> number of edges between the clusters
    ygm::container::map<std::pair<cluster_id_t, cluster_id_t>, uint32_t>
                    intercluster_edge_map(world);
    static uint64_t num_within_cluster_edges;
    num_within_cluster_edges = 0;
    static uint64_t num_between_cluster_edges;
    num_between_cluster_edges = 0;

    // Map of the neighborhood purity for each node
    ygm::container::map<id_t, float> neighborhood_purity_map(world);

    {
      world.cout0()
          << "\n<<Get neighborhood purities and count intercluster edges>>"
          << std::endl;
      step_timer.reset();

      auto process_neighbors_lambda = [&intercluster_edge_map,
                                       &neighborhood_purity_map](
                                          [[maybe_unused]] const id_t &node,
                                          const auto &cluster_info_pair) {
        cluster_id_t cluster_id                    = cluster_info_pair.first;
        uint32_t     num_neighbors                 = 0;
        uint32_t     num_neighbors_in_same_cluster = 0;
        for (const auto &neighbor_cluster_id : cluster_info_pair.second) {
          // Update number of neighbors and number of neighbors in the same
          // cluster
          ++num_neighbors;
          if (neighbor_cluster_id == cluster_id) {
            ++num_neighbors_in_same_cluster;
          }
          // Always have the first entry of the intercluster edge be <= the
          // second entry to avoid double counting
          // Each edge will show up twice - once with each entry as key
          // (assume no self loops)
          if (cluster_id <= neighbor_cluster_id) {
            if (cluster_id == neighbor_cluster_id) {
              ++num_within_cluster_edges;
            } else {
              ++num_between_cluster_edges;
            }

            std::pair<cluster_id_t, cluster_id_t> cluster_pair = {
                cluster_id, neighbor_cluster_id};
            intercluster_edge_map.async_visit(
                cluster_pair,
                [](const std::pair<cluster_id_t, cluster_id_t> &cluster_pair,
                   auto &count) { ++count; });
          }
        }

        // Calculate the neighborhood purity
        float neighborhood_purity =
            static_cast<float>(num_neighbors_in_same_cluster) /
            static_cast<float>(num_neighbors);
        neighborhood_purity_map.async_insert(node, neighborhood_purity);
      };

      cluster_adjacency_map.for_all(process_neighbors_lambda);
      world.barrier();

      // Calculate the average neighborhood purity
      step_timer.reset();
      static float purity_sum;
      purity_sum                 = 0.0;
      auto add_purity_sum_lambda = []([[maybe_unused]] const id_t &node,
                                      const float &neighborhood_purity) {
        purity_sum += neighborhood_purity;
      };
      neighborhood_purity_map.for_all(add_purity_sum_lambda);

      float    total_purity_sum = ygm::sum(purity_sum, world);
      uint64_t num_nodes        = neighborhood_purity_map.size();
      world.cout0() << "Sum of neighborhood purities: " << purity_sum
                    << std::endl;
      world.cout0() << "Number of nodes: " << num_nodes << std::endl;
      world.cout0() << "Average neighborhood purity: "
                    << total_purity_sum / static_cast<float>(num_nodes)
                    << std::endl;

      world.cout0() << "intercluster_edge_map.size() = "
                    << intercluster_edge_map.size() << std::endl;
      world.cout0() << "Number of intracluster edges: "
                    << ygm::sum(num_within_cluster_edges, world) << std::endl;
      world.cout0() << "Number of intercluster edges: "
                    << ygm::sum(num_between_cluster_edges, world) << std::endl;
      world.cout0()
          << "Time to get neighborhood purities and count intercluster "
             "edges (s):"
          << step_timer.elapsed() << std::endl;
    }

    // Map of the edge probabilities for each cluster pair
    // (cluster i, cluster j)
    //      -> fraction of possible edges from node in i to node
    ygm::container::map<std::pair<cluster_id_t, cluster_id_t>, double>
        edge_prob_map(world);

    // Map of cluster -> size (number of nodes)
    ygm::container::map<cluster_id_t, uint32_t> cluster_size_map(world);
    auto cluster_size_map_ptr = cluster_size_map.get_ygm_ptr();

    if (opt.write_cluster_size_file || opt.write_edge_prob_file) {
      world.cout0() << "\n<<Get cluster sizes and maybe edge probabilities>>"
                    << std::endl;
      step_timer.reset();

      // Get the cluster sizes
      auto increase_cluster_size_lambda =
          [&cluster_size_map]([[maybe_unused]] const id_t &point_id,
                              const cluster_id_t          &cluster_id) {
            cluster_size_map.async_visit(
                cluster_id, []([[maybe_unused]] const id_t &cluster_id,
                               auto                        &size) { ++size; });
          };
      point_cluster_map.for_all(increase_cluster_size_lambda);
      world.barrier();

      world.cout0() << "Number of clusters: " << cluster_size_map.size()
                    << std::endl;
    }

    if (opt.write_edge_prob_file) {
      // Get edge probabilities
      auto edge_prob_map_ptr = edge_prob_map.get_ygm_ptr();

      auto get_edge_prob_lambda =
          [&cluster_size_map_ptr, &edge_prob_map_ptr](
              const std::pair<cluster_id_t, cluster_id_t> &cluster_id_pair,
              const auto                                  &number_of_edges) {
            // Within cluster edge probability is #edges / (cluster size
            // choose 2)
            if (cluster_id_pair.first == cluster_id_pair.second) {
              auto intracluster_edge_lambda =
                  [](const cluster_id_t &cluster_id, const auto &cluster_size,
                     const auto &number_of_edges, auto &edge_prob_map_ptr) {
                    double edge_prob =
                        static_cast<double>(number_of_edges) /
                        static_cast<double>(cluster_size * (cluster_size - 1));
                    std::pair<cluster_id_t, cluster_id_t> cluster_id_pair = {
                        cluster_id, cluster_id};
                    edge_prob_map_ptr->async_insert(cluster_id_pair, edge_prob);
                  };
              cluster_size_map_ptr->async_visit(
                  cluster_id_pair.first, intracluster_edge_lambda,
                  number_of_edges, edge_prob_map_ptr);
            } else {
              auto intercluster_edge_lambda =
                  []([[maybe_unused]] const cluster_id_t &cluster_id,
                     const auto                          &cluster_size,
                     const std::pair<cluster_id_t, cluster_id_t>
                                &cluster_id_pair,
                     const auto &number_of_edges, auto &edge_prob_map_ptr,
                     auto &cluster_size_map_ptr) {
                    double partial_edge_prob =
                        static_cast<double>(number_of_edges) /
                        static_cast<double>(cluster_size);
                    auto visit_other_cluster_lambda =
                        []([[maybe_unused]] const cluster_id_t &cluster_id,
                           auto                                &cluster_size,
                           const std::pair<cluster_id_t, cluster_id_t>
                                        &cluster_id_pair,
                           const double &partial_edge_prob,
                           auto         &edge_prob_map_ptr) {
                          double edge_prob = partial_edge_prob /
                                             static_cast<double>(cluster_size);
                          edge_prob_map_ptr->async_insert(cluster_id_pair,
                                                          edge_prob);
                        };
                    cluster_size_map_ptr->async_visit(
                        cluster_id_pair.second, visit_other_cluster_lambda,
                        cluster_id_pair, partial_edge_prob, edge_prob_map_ptr);
                  };
              cluster_size_map_ptr->async_visit(
                  cluster_id_pair.first, intercluster_edge_lambda,
                  cluster_id_pair, number_of_edges, edge_prob_map_ptr,
                  cluster_size_map_ptr);
            }
          };
      intercluster_edge_map.for_all(get_edge_prob_lambda);
      world.barrier();
    }

    /* Write to files */

    // YGM multi output for writing to file, append = false
    // TODO: modify this for multiple ranks writing to one file
    ygm::io::multi_output mo(world, opt.output_dir, 1024, false);

    {
      world.cout0() << "\n<<Write files>>" << std::endl;
      step_timer.reset();

      if (opt.write_neighborhood_purity_file) {
        if (!opt.sort_outputs) {
          if (world.rank() == 0) {
            mo.async_write_line("neighborhood_purity.txt",
                                "node\tneighborhood_purity");
          }
          world.barrier();
          auto write_neighborhood_purity_lambda =
              [&mo](const id_t &node, const float &neighborhood_purity) {
                std::stringstream ss;
                ss << node << "\t" << neighborhood_purity;
                mo.async_write_line("neighborhood_purity.txt", ss.str());
              };
          neighborhood_purity_map.for_all(write_neighborhood_purity_lambda);
        } else {
          using key_type   = id_t;
          using value_type = float;
          static std::vector<std::pair<key_type, value_type>> local_vector;
          static std::vector<std::pair<key_type, value_type>> global_vector;
          auto consolidate_lambda = [](const id_t &key, float &value) {
            local_vector.push_back(std::make_pair(key, value));
          };
          neighborhood_purity_map.for_all(consolidate_lambda);
          world.barrier();

          // Broadcast the vector to all ranks
          world.async_bcast(
              [](const auto &local_vector) {
                global_vector.insert(global_vector.end(), local_vector.begin(),
                                     local_vector.end());
              },
              local_vector);
          world.barrier();

          if (world.rank() == 0) {
            std::sort(global_vector.begin(), global_vector.end(),
                      [](const auto &a, const auto &b) {
                        return a.second > b.second;
                      });

            std::string out_file = opt.output_dir + "/neighborhood_purity.txt";
            std::ofstream ofs(out_file);
            ofs << "node\tneighborhood_purity" << std::endl;
            for (auto &entry : global_vector) {
              ofs << entry.first << "\t" << entry.second << std::endl;
            }
          }
        }
        //  std::string output_filename =
        //      opt.output_dir + "/neighborhood_purity.txt";
        //  print_ygm_map_tsv_sorted_by_value(neighborhood_purity_map, world,
        //                                    output_filename);
      }

      if (opt.write_cluster_size_file) {
        if (!opt.sort_outputs) {
          if (world.rank() == 0) {
            mo.async_write_line("cluster_sizes.txt", "cluster_id\tnum_nodes");
          }
          world.barrier();
          auto write_cluster_size_lambda = [&mo](const cluster_id_t &cluster_id,
                                                 auto               &size) {
            std::stringstream ss;
            ss << cluster_id << "\t" << size;
            mo.async_write_line("cluster_sizes.txt", ss.str());
          };
          cluster_size_map.for_all(write_cluster_size_lambda);
        } else {
          using key_type   = cluster_id_t;
          using value_type = uint32_t;
          static std::vector<std::pair<key_type, value_type>> local_vector;
          static std::vector<std::pair<key_type, value_type>> global_vector;
          auto consolidate_lambda = [](const key_type   &key,
                                       const value_type &value) {
            local_vector.push_back(std::make_pair(key, value));
          };
          cluster_size_map.for_all(consolidate_lambda);
          world.barrier();

          // Broadcast the vector to all ranks
          world.async_bcast(
              [](const auto &local_vector) {
                global_vector.insert(global_vector.end(), local_vector.begin(),
                                     local_vector.end());
              },
              local_vector);
          world.barrier();

          if (world.rank() == 0) {
            std::sort(global_vector.begin(), global_vector.end(),
                      [](const auto &a, const auto &b) {
                        return a.second > b.second;
                      });

            std::string   out_file = opt.output_dir + "/cluster_sizes.txt";
            std::ofstream ofs(out_file);
            ofs << "cluster_id\tnum_nodes" << std::endl;
            for (const auto &entry : global_vector) {
              ofs << entry.first << "\t" << entry.second << std::endl;
            }
          }
        }
      }

      if (opt.write_degree_file) {
        if (!opt.sort_outputs) {
          if (world.rank() == 0) {
            mo.async_write_line("degrees.txt", "node\tdegree");
          }
          world.barrier();
          auto write_degree_lambda =
              [&mo](const id_t &node,
                    std::pair<cluster_id_t, std::vector<cluster_id_t>>
                        &cluster_info_pair) {
                std::stringstream ss;
                ss << node << "\t" << cluster_info_pair.second.size();
                mo.async_write_line("degrees.txt", ss.str());
              };
          cluster_adjacency_map.for_all(write_degree_lambda);
        } else {
          using key_type   = id_t;
          using value_type = std::pair<cluster_id_t, std::vector<cluster_id_t>>;
          static std::vector<std::pair<key_type, value_type>> local_vector;
          static std::vector<std::pair<key_type, value_type>> global_vector;
          auto consolidate_lambda = [](const key_type   &key,
                                       const value_type &value) {
            local_vector.push_back(std::make_pair(key, value));
          };
          cluster_adjacency_map.for_all(consolidate_lambda);
          world.barrier();

          // Broadcast the vector to all ranks
          world.async_bcast(
              [](const auto &local_vector) {
                global_vector.insert(global_vector.end(), local_vector.begin(),
                                     local_vector.end());
              },
              local_vector);
          world.barrier();

          if (world.rank() == 0) {
            std::sort(global_vector.begin(), global_vector.end(),
                      [](const auto &a, const auto &b) {
                        return a.second.second.size() > b.second.second.size();
                      });

            std::string   out_file = opt.output_dir + "/degrees.txt";
            std::ofstream ofs(out_file);
            ofs << "node\tdegree" << std::endl;
            for (auto &entry : global_vector) {
              ofs << entry.first << "\t" << entry.second.second.size()
                  << std::endl;
            }
          }
        }
      }

      if (world.rank() == 0) {
        mo.async_write_line("intercluster_edge_counts.txt",
                            "cluster_id1\tcluster_id2\tnum_edges");
      }
      world.barrier();

      auto write_intercluster_edge_counts_lambda =
          [&mo](const std::pair<cluster_id_t, cluster_id_t> &cluster_id_pair,
                const auto                                  &count) {
            std::stringstream ss;
            ss << cluster_id_pair.first << "\t" << cluster_id_pair.second
               << "\t" << count;
            mo.async_write_line("intercluster_edge_counts.txt", ss.str());
          };
      intercluster_edge_map.for_all(write_intercluster_edge_counts_lambda);

      if (opt.write_edge_prob_file) {
        static std::vector<std::tuple<cluster_id_t, cluster_id_t, double>>
            local_edge_prob_vector;
        static std::vector<std::tuple<cluster_id_t, cluster_id_t, double>>
            edge_prob_vector;

        auto consolidate_edge_prob_lambda =
            [](const std::pair<cluster_id_t, cluster_id_t> &cluster_id_pair,
               const double                                &edge_prob) {
              std::tuple<cluster_id_t, cluster_id_t, double> edge_prob_entry = {
                  cluster_id_pair.first, cluster_id_pair.second, edge_prob};
              local_edge_prob_vector.push_back(edge_prob_entry);
            };
        edge_prob_map.for_all(consolidate_edge_prob_lambda);
        world.barrier();

        // Broadcast the edge_probability vector to all ranks
        world.async_bcast(
            [](const auto &local_vector) {
              edge_prob_vector.insert(edge_prob_vector.end(),
                                      local_vector.begin(), local_vector.end());
            },
            local_edge_prob_vector);
        world.barrier();

        if (world.rank() == 0) {
          std::sort(edge_prob_vector.begin(), edge_prob_vector.end(),
                    [](const auto &a, const auto &b) {
                      return std::get<2>(a) > std::get<2>(b);
                    });

          std::string edge_prob_file =
              opt.output_dir + "/edge_probabilities.txt";
          std::ofstream ofs(edge_prob_file);
          ofs << "cluster_id1\tcluster_id2\tedge_prob" << std::endl;
          for (auto &entry : edge_prob_vector) {
            ofs << std::get<0>(entry) << "\t" << std::get<1>(entry) << "\t"
                << std::get<2>(entry) << std::endl;
          }
        }
      }

      world.barrier();
      world.cout0() << "Time to write files (s):\t" << step_timer.elapsed()
                    << std::endl;
    }
  }

  return 0;
}

bool parse_options(int argc, char **argv, option_t &opt, bool &help) {
  opt.clustering_path.clear();
  opt.edge_list_path.clear();
  help = false;

  int opt_char;
  while ((opt_char = ::getopt(argc, argv, "p:g:o:r:wecsdvh")) != -1) {
    switch (opt_char) {
      case 'p':
        opt.clustering_path = std::filesystem::path(optarg);
        break;

      case 'g':
        opt.edge_list_path = std::filesystem::path(optarg);
        break;

      case 'o':
        opt.output_dir = optarg;
        break;

      case 'r':
        opt.file_of_nodes_to_ignore = std::filesystem::path(optarg);
        break;

      case 'e':
        opt.write_edge_prob_file = true;
        break;

      case 'w':
        opt.write_neighborhood_purity_file = true;
        break;

      case 'c':
        opt.write_cluster_size_file = true;
        break;

      case 's':
        opt.sort_outputs = true;
        break;

      case 'd':
        opt.write_degree_file = true;
        break;

      case 'h':
        help = true;
        return true;

      default:
        help = true;
        std::cout << "Unrecognized option: " << opt_char << ", ignore."
                  << std::endl;
    }
  }

  return true;
}

template <typename cout_type>
void show_options(const option_t &opt, cout_type &cout) {
  cout << "Options:" << std::endl;
  cout << "  Clustering path: " << opt.clustering_path << std::endl;
  cout << "  Edge list path: " << opt.edge_list_path << std::endl;
  cout << "  Output directory:  " << opt.output_dir << std::endl;
  cout << "  write neighborhood purity to file: "
       << opt.write_neighborhood_purity_file << std::endl;
  cout << "  write edge probability file: " << opt.write_edge_prob_file
       << std::endl;
  cout << "  write cluster size file: " << opt.write_cluster_size_file
       << std::endl;
  cout << "  write degree file: " << opt.write_degree_file << std::endl;
  cout << "  sort outputs: " << opt.sort_outputs << std::endl;
  cout << "  File of nodes to ignore:  " << opt.file_of_nodes_to_ignore
       << std::endl;
}
