#include "sa_seq.hpp"

#include <mpi.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "data_source.h"
using std::vector, std::pair;

#define SA(i) (B[i].second)

const size_t char_size = 3;  // There are 4 characters -- A, C, T, G and one
                             // special character -- the end of the word
static const size_t k = (64 - char_size) / char_size;
static inline uint64_t char_to_word(char c) {
    switch (c) {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 3;
        case 'T':
            return 4;
    }
    return 0;
}

static inline void non_seq_fail() {
    std::cerr << "This is a sequential version, which requires to be "
                 "executed on exactly one processor with access to the "
                 "whole input."
              << std::endl;
    exit(1);
}

inline bool rebucket_and_check_all_singleton(
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B) {
    bool res = true;
    std::pair<uint64_t, uint64_t> last_val = B[0].first;
    uint64_t g = 1;
    B[0].first = std::make_pair(g, 0);
    for (size_t i = 1; i < B.size(); i++) {
        if (last_val == B[i].first) {
            res = false;
        } else {
            last_val = B[i].first;
            g++;
        }
        B[i].first = std::make_pair(g, 0);
    }
    return res;
}

inline void printB(
    const std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const char *buffer) {
    for (uint64_t i = 0; i < B.size(); i++) {
        std::cerr << B[i].first.first << " ";
        if (B[i].second >= 10 && B[i].first.first < 10) std::cerr << " ";
    }
    std::cerr << std::endl;
    for (uint64_t i = 0; i < B.size(); i++) {
        std::cerr << B[i].first.second << " ";
        if (B[i].second >= 10 && B[i].first.second < 10) std::cerr << " ";
    }
    std::cerr << std::endl;
    for (uint64_t i = 0; i < B.size(); i++) std::cerr << B[i].second << " ";
    std::cerr << std::endl << std::endl;

    for (uint64_t i = 0; i < B.size(); i++) {
        std::cerr << &buffer[B[i].second] << std::endl;
    }
    std::cerr << std::endl << std::endl;
}

const std::vector<uint64_t> sa_word_size_param(
    int my_rank, int number_of_processes, int which, DataSource &data_source,
    std::vector<std::string> queries) {
    const uint64_t genome_size = data_source.getTotalGenomeSize(which),
                   my_genome_part_size = data_source.getNodeGenomeSize(which),
                   my_genome_offset = data_source.getNodeGenomeOffset(which);

    if (number_of_processes != 1 || my_rank != 0 ||
        genome_size != my_genome_part_size || my_genome_offset != 0)
        non_seq_fail();

    std::string buffer(my_genome_part_size, 0);
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> B(
        my_genome_part_size);
    std::vector<uint64_t> B_prim(my_genome_part_size);

    data_source.getNodeGenomeValues(which, buffer.data());

    // Create B with k-mers
    uint64_t M = 1;
    uint64_t current_value = 0;
    for (uint64_t i = 0; i < k; i++) {
        current_value *= (1 << char_size);
        if (i < buffer.size()) current_value += char_to_word(buffer[i]);
        if (i > 0) M *= (1 << char_size);
        // std::cerr << M << std::endl;
    }
    for (size_t i = 0; i < std::min(my_genome_part_size, B.size()); i++) {
        B[i] = std::make_pair(std::make_pair(current_value, 0), i);
        current_value %= M;
        current_value *= (1 << char_size);
        if (i + k < buffer.size()) current_value += char_to_word(buffer[i + k]);
    }

    // Create SA
    std::sort(B.begin(), B.end());
    bool done = rebucket_and_check_all_singleton(B);
    for (uint64_t h = k;; h <<= 1) {
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            B_prim[B[i].second] = B[i].first.first;
        }
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            B[i].first.first = B_prim[i];
        }
        if (done) {
            break;
        }
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            if (i + h < B.size())
                B[i].first.second = B[i + h].first.first;
            else
                B[i].first.second = 0;
            B[i].second = i;
        }
        std::sort(B.begin(), B.end());
        done = rebucket_and_check_all_singleton(B);
    }

    // Answer the queries
    std::vector<uint64_t> res(queries.size());
    for (uint64_t i = 0; i < queries.size(); i++) {
        uint64_t first_occurrence, after_last_occurrence;
        // find first occurrence
        {
            uint64_t b = 0, e = B.size(), m;
            while (b < e) {
                m = (b + e) / 2;
                if (strcmp(queries[i].c_str(), &buffer.c_str()[B[m].second]) >
                    0) {  // queries[i] > &buffer.c_str()[B[m].second]
                    b = m + 1;
                } else {
                    e = m;
                }
            }
            first_occurrence = b;
        }
        queries[i][queries[i].size() - 1]++;
        {
            uint64_t b = 0, e = B.size(), m;
            while (b < e) {
                m = (b + e) / 2;
                if (strcmp(queries[i].c_str(), &buffer.c_str()[B[m].second]) >
                    0) {  // queries[i] > &buffer.c_str()[B[m].second]
                    b = m + 1;
                } else {
                    e = m;
                }
            }
            after_last_occurrence = b;
        }
        queries[i][queries[i].size() - 1]--;
        res[i] = after_last_occurrence - first_occurrence;
    }
    return res;
}

void sa(int my_rank, int number_of_processes, uint64_t n, uint64_t m,
        DataSource &data_source, const std::string &queries_in,
        const std::string &queries_out) {
    // Read queries
    std::vector<std::string> queries(m);
    {
        std::ifstream queries_in_file(queries_in);
        for (uint64_t i = 0; i < m; i++) queries_in_file >> queries[i];
    }

    // Compute SA and answer the queries
    std::vector<std::vector<uint64_t>> res(n);
    for (uint64_t i = 0; i < n; i++) {
        res[i] = sa_word_size_param(my_rank, number_of_processes, (int)i,
                                    data_source, queries);
    }

    // Write the results to a file
    {
        std::ofstream queries_out_file(queries_out);
        for (uint64_t j = 0; j < m; j++) {
            for (uint64_t i = 0; i < n; i++) {
                if (i > 0) queries_out_file << " ";
                queries_out_file << res[i][j];
            }
            queries_out_file << "\n";
        }
    }
}
