#include <mpi.h>

#include <iostream>

#include "params.hpp"

using std::cerr;

int main(int argc, char* argv[]) {
    auto [n, m, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    MPI_Init(&argc, &argv); /* intialize the library with parameters caught by
                               the runtime */
    MPI_Finalize();
    return 0;
}
