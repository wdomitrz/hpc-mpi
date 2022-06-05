#include <mpi.h>

#include <iostream>

#include "data_source.h"
#include "params.hpp"
#include "sa.hpp"

using std::cerr;

int main(int argc, char* argv[]) {
    int number_of_processes;
    int my_rank;
    auto [n, m, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    DataSource data_source(genome_in);

    sa(my_rank, number_of_processes, n, m, data_source, queries_in,
       queries_out);

    MPI_Finalize();
    return 0;
}
