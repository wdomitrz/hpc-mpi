#include <iostream>

#include "params.hpp"

using std::cout;

int main(int argc, char* argv[]) {
    auto [n, m, genome_in, queries_in, queries_out] = parse_args(argc, argv);

    cout << n << "\t" << m << "\t" << genome_in << "\t" << queries_in << "\t"
         << queries_out << "\n";

    return 0;
}
