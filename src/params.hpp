#ifndef PARAMS_H
#define PARAMS_H
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>

static const int EXPECTED_NUMBER_OF_ARGUMENTS = 5;
static const int ARG_N_IDX = 1, ARG_M_IDX = 2, ARG_GENOME_IN_IDX = 3,
                 ARG_QUERIES_IN_IDX = 4, ARG_QUERIES_OUT_IDX = 5;

std::tuple<const int, const int, char *, std::string, std::string> parse_args(
    int &argc, char *argv[]) {
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

    char **original_argv = argv;
    argv += EXPECTED_NUMBER_OF_ARGUMENTS;
    argv[0] = original_argv[0];
    argc -= EXPECTED_NUMBER_OF_ARGUMENTS;

    return std::make_tuple(
        std::stoi(original_argv[ARG_N_IDX]),
        std::stoi(original_argv[ARG_M_IDX]), original_argv[ARG_GENOME_IN_IDX],
        original_argv[ARG_QUERIES_IN_IDX], original_argv[ARG_QUERIES_OUT_IDX]);
}
#endif /* PARAMS_H */
