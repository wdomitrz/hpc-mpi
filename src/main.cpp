#include <mpi.h>

#include <iostream>

#include "data_source.h"
#include "params.hpp"

using std::cerr;

int main(int argc, char* argv[]) {
    int number_of_processes;
    int my_rank;
    auto [n, m, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    MPI_Init(&argc, &argv); /* intialize the library with parameters caught by
                               the runtime */
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    DataSource ds(genome_in);
    for (int i = 0; i < n; i++) {
        auto x = ds.getTotalGenomeSize(i), y = ds.getNodeGenomeSize(i),
             z = ds.getNodeGenomeOffset(i);
        for (int j = 0; j < number_of_processes; ++j) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (j == my_rank) {
                cerr << my_rank << "/" << number_of_processes << "\t" << i
                     << "\t" << x << "\t" << y << "\t" << z << "\n";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
