#ifndef PARAMS_H
#define PARAMS_H
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>

static const int EXPECTED_NUMBER_OF_ARGUMENTS = 5;
static const int ARG_N_IDX = 1, ARG_M_IDX = 2, ARG_GENOME_IN_IDX = 3,
                 ARG_QUERIES_IN_IDX = 4, ARG_QUERIES_OUT_IDX = 5;

auto parse_args(int &argc, char **(&argv)) {
    if (argc < 1 + EXPECTED_NUMBER_OF_ARGUMENTS) {
        std::cerr << "Expected at least " << EXPECTED_NUMBER_OF_ARGUMENTS
                  << " parameters in the following format:"
                  << "\n"
                  << "\t"
                  << "mpiexec ./genome_index :n :m :genome_in :queries_in "
                     ":queries_out"
                  << std::endl;
        exit(1);
    }

    const auto n = std::stoi(argv[ARG_N_IDX]);
    const auto m = std::stoi(argv[ARG_M_IDX]);
    const auto genome_in = argv[ARG_GENOME_IN_IDX];
    const auto queries_in = std::string(argv[ARG_QUERIES_IN_IDX]);
    const auto queries_out = std::string(argv[ARG_QUERIES_OUT_IDX]);

    argv[EXPECTED_NUMBER_OF_ARGUMENTS] = argv[0];
    argv += EXPECTED_NUMBER_OF_ARGUMENTS;
    argc -= EXPECTED_NUMBER_OF_ARGUMENTS;

    return std::make_tuple(n, m, genome_in, queries_in, queries_out);
}

#endif /* PARAMS_H */
