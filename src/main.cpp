#include <mpi.h>

#include <iostream>

#include "params.hpp"

using std::cerr;

int main(int argc, char* argv[]) {
    int number_of_processes;
    int myRank;
    auto [n, m, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    MPI_Init(&argc, &argv); /* intialize the library with parameters caught by
                               the runtime */

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    MPI_Finalize();
    return 0;
}
