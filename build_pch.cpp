#include "build_pch.hpp"

#include "dijkstra.hpp"
#include "query.hpp"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr,
            "WeightedGraph only\n"
            "Usage: %s [-i input_graph]\n"
            "Options:\n"
            "\t-o,\toutput_ch_graph\n"
            "\t-p,\tmax_pop_count\n"
            "\t-s,\tselect fraction\n"
            "\t-t,\ts-t query verify numn"
            "\t-q,\tsssp query verify num\n"
            "\t-b,\tdegree bound\n"
            "\t-d,\tprint per round detail\n",
            argv[0]);
    return 0;
  }
  char c;
  int max_pop_count = 500;
  int bidirect_verify_num = 0;
  int sssp_verify_num = 0;
  int degree_bound = 0;
  bool degree_bounded = false, print_detail = false, write_ch = false;
  double sample_bound = 1;
  EdgeTy k_value = std::numeric_limits<EdgeTy>::max();
  std::vector<size_t> vertex_label_offsets = {};
  while ((c = getopt(argc, argv, "i:o:p:s:t:q:b:dk:v:")) != -1) {
    switch (c) {
      case 'i':
        INPUT_FILEPATH = optarg;
        break;
      case 'o':
        write_ch= true;
        OUTPUT_FILEPATH = optarg;
        break;
      case 'p':
        max_pop_count = atol(optarg);
        if (max_pop_count < 0) {
          fprintf(stderr, "Error: max_pop_count must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 's':
        sample_bound = atof(optarg);
        if (sample_bound <= 0 || sample_bound >1) {
          fprintf(stderr, "Error: selection_fraction must be larger than 0 and smaller than or equal to 1\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 't':
        bidirect_verify_num = atol(optarg);
        if (bidirect_verify_num < 0) {
          fprintf(stderr, "Error: bidirect_verify_num must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'q':
        sssp_verify_num = atol(optarg);
        if (sssp_verify_num < 0) {
          fprintf(stderr, "Error: sssp_verify_num must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'b':
        degree_bounded = true;
        degree_bound = atol(optarg);
        if (degree_bound < 0) {
          fprintf(stderr, "Error: degree_bound must be non negative\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'd':
        print_detail = true;
        break;
      case 'k':
        k_value = atof(optarg);  // Parse k as a floating point value (EdgeTy is assumed to be a float/double)
        if (k_value < 0) {
          fprintf(stderr, "Error: k (shortcut weight bound) must be non-negative\n");
          exit(EXIT_FAILURE);
        }
      case 'v':
        LABEL_OFFSETS_FILEPATH = optarg;
        break;
      break;
      default:
        fprintf(stderr, "Error: Unknown option %c\n", optopt);
        exit(EXIT_FAILURE);
    }
  }
  Graph origin_graph = read_graph(INPUT_FILEPATH);

  if (LABEL_OFFSETS_FILEPATH != nullptr) {
    std::ifstream label_offsets_file(LABEL_OFFSETS_FILEPATH);  // Open the file for reading
    if (!label_offsets_file) {
        std::cerr << "Error opening vertex labels file." << std::endl;
    } else {
        size_t label_offset;
        while (label_offsets_file >> label_offset) { // Read space-separated sequences as single strings
            vertex_label_offsets.push_back(label_offset);
        }
        label_offsets_file.close();
    }
  }
  std::vector<size_t> label_lengths;
  if (!vertex_label_offsets.empty()) {
      label_lengths.reserve(vertex_label_offsets.size() - 1); // Reserve for efficiency
      for (size_t i = 1; i < vertex_label_offsets.size(); ++i) {
          label_lengths.push_back(vertex_label_offsets[i] - vertex_label_offsets[i - 1]);
      }
  }


  PCH *solver =
      new PCH(origin_graph, max_pop_count, degree_bounded, degree_bound, sample_bound, k_value, print_detail, label_lengths);
  PchGraph contracted_graph = solver->createContractionHierarchy();
  delete (solver);
  if(!label_lengths.empty() && 1 < *std::max_element(label_lengths.begin(), label_lengths.end())) {
    printf("skipping queries because of string labels longer than one character\n");
    bidirect_verify_num=0;
    sssp_verify_num=0;
  }

  PchQuery query(contracted_graph, origin_graph);
  ofstream ofs("pch.tsv", ios::app);
  ofs << fixed << setprecision(6);
  printf("Start query\n");
  query.make_inverse();
  for (int i = 0; i < bidirect_verify_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    NodeId t = hash32(i + s) % origin_graph.n;
    internal::timer tm;
    auto [d, itr] = query.stQuery(s, t);
    tm.stop();
    internal::timer tm2;
    auto [exp_dist, itr2] = query.stVerifier(s, t);
    tm2.stop();
    if (exp_dist != d) {
      printf("Error: s: %u, t: %u, output_dist: %u, exp_dist: %u\n", s, t, d,
             exp_dist);
      abort();
    }
    printf("s: %u, t: %u, itr: %u, d: %u\n", s, t, itr, d);
    ofs << s << '\t' << t << '\t' << itr << '\t' << itr2 << '\t' << d << '\t'
        << tm.total_time() << '\t' << tm2.total_time() << '\n';
  }
  for (int i = 0; i < sssp_verify_num; i++) {
    NodeId s = hash32(i) % origin_graph.n;
    internal::timer tm;
    auto dist = query.ssspQuery(s);
    tm.stop();
    printf("%d SSSP query: s: %u\n", i, s);
    internal::timer tm2;
    auto exp_dist = query.ssspVerifier(s);
    tm2.stop();
    if (exp_dist != dist) {
      for (size_t i = 0; i < dist.size(); i++) {
        if (dist[i] != exp_dist[i]) {
          printf("output_dist[%zu]: %u, exp_dist[%zu]: %u\n", i, dist[i], i,
                 exp_dist[i]);
        }
      }
      abort();
    }
    ofs << s << '\t' << tm.total_time() << '\t' << tm2.total_time() << '\n';
  }
  ofs.close();
  if(write_ch)write_pbbs_format(contracted_graph, OUTPUT_FILEPATH);
  return 0;
}