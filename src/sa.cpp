#include "sa.hpp"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "data_source.h"
using std::vector, std::pair;

#define K_VAL (std::min(k_max, genome_size / number_of_processes))

#define whose(rank) (rank * number_of_processes / genome_size)
#define how_much_x_has(rank) \
    (how_much_node_has(rank, number_of_processes, genome_size))
#define offset(rank) \
    (get_node_genome_offset(rank, number_of_processes, genome_size))
#define my_sort(B)                                             \
    (my_sort_params(my_rank, number_of_processes, genome_size, \
                    my_genome_part_size, B))
#define printB() printB_fun(B, buffer.c_str(), my_genome_part_size)
#define ok() \
    { std::cerr << "ok:\t" << my_rank << "\t" << __LINE__ << std::endl; }

MPI_Request global_request;
const size_t char_size = 3;  // There are 4 characters -- A, C, T, G and one
                             // special character -- the end of the word
static const size_t k_max = (64 - char_size) / char_size;
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

uint64_t how_much_node_has(int rank, int nprocs, uint64_t totalSize) {
    uint64_t nodeSize = (totalSize + nprocs - rank - 1) / nprocs;
    return nodeSize;
}

inline static uint64_t get_node_genome_offset(int rank, int nprocs,
                                              uint64_t totalSize) {
    MPI_Offset offset = (uint64_t)rank * (totalSize / (uint64_t)nprocs) +
                        std::min((uint64_t)totalSize % nprocs, (uint64_t)rank);
    return offset;
}

inline bool rebucket_and_check_all_singleton(
    int my_rank, int number_of_processes, const uint64_t,
    const uint64_t my_genome_part_size, const uint64_t,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B) {
    bool my_res = true;
    std::vector<uint64_t> my_partial_results(my_genome_part_size);

    if (my_rank > 0)
        MPI_Isend(&B[0], sizeof(B[0]), MPI_BYTE, my_rank - 1, 0, MPI_COMM_WORLD,
                  &global_request);  // Instead of creating a custom MPI type, I
                                     // use MPI_BYTE and send data as bytes.

    MPI_Request get_next_one_request;
    if (my_rank < number_of_processes - 1)
        MPI_Irecv(
            &B[my_genome_part_size], sizeof(B[0]), MPI_BYTE, my_rank + 1,
            MPI_ANY_TAG, MPI_COMM_WORLD,
            &get_next_one_request);  // Instead of creating a custom MPI type, I
                                     // use MPI_BYTE and send data as bytes.

    uint64_t my_count = 0;
    std::pair<uint64_t, uint64_t> prev_val = B[0].first;
    B[0].first = std::make_pair(0, 0);
    for (size_t i = 1; i < my_genome_part_size; i++) {
        if (prev_val == B[i].first) {
            B[i].first = B[i - 1].first;
            my_res = false;
        } else {
            B[i].first = B[i - 1].first;
            B[i].first.first++;
            my_count++;
        }
    }

    if (my_rank + 1 != number_of_processes) {
        MPI_Wait(&get_next_one_request, nullptr);
        if (prev_val == B[my_genome_part_size].first) {
            my_res = false;
        } else {
            my_count++;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &my_res, 1, MPI_C_BOOL, MPI_LAND,
                  MPI_COMM_WORLD);
    MPI_Exscan(MPI_IN_PLACE, &my_count, 1, MPI_UINT64_T, MPI_SUM,
               MPI_COMM_WORLD);

    if (my_rank != 0)
        for (size_t i = 0; i < my_genome_part_size; i++)
            B[i].first.first += my_count;

    return my_res;
}

inline void my_sort_params(
    const int my_rank, const int number_of_processes,
    const uint64_t genome_size, const uint64_t my_genome_part_size,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B) {
    // std::sort(B.data(), &B.data()[my_genome_part_size]);

    std::sort(B.data(), &B.data()[number_of_processes]);

    if (my_rank == 0) {
        std::vector<int> recv_size(number_of_processes);
        std::vector<int> recv_offset(number_of_processes);
        recv_size[0] = (int)how_much_x_has(0) * sizeof(B[0]);
        recv_offset[0] = 0;
        for (int i = 1; i < number_of_processes; i++) {
            recv_offset[i] = recv_offset[i - 1] + recv_size[i - 1];
            recv_size[i] = (int)how_much_x_has(i) * sizeof(B[0]);
        }

        static std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>>
            B_all(genome_size);

        MPI_Gatherv(B.data(), (int)sizeof(B[0]) * (int)my_genome_part_size,
                    MPI_BYTE, B_all.data(), recv_size.data(),
                    recv_offset.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        std::sort(B_all.begin(), B_all.end());
        MPI_Scatterv(B_all.data(), recv_size.data(), recv_offset.data(),
                     MPI_BYTE, B.data(),
                     (int)sizeof(B[0]) * (int)my_genome_part_size, MPI_BYTE, 0,
                     MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(B.data(), (int)sizeof(B[0]) * (int)my_genome_part_size,
                    MPI_BYTE, nullptr, nullptr, nullptr, MPI_BYTE, 0,
                    MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_BYTE, B.data(),
                     (int)sizeof(B[0]) * (int)my_genome_part_size, MPI_BYTE, 0,
                     MPI_COMM_WORLD);
    }
}

inline void printB_fun(
    const std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const char *, uint64_t genome_size) {
    for (uint64_t i = 0; i < genome_size; i++) {
        std::cerr << B[i].first.first << " ";
        if (B[i].second >= 10 && B[i].first.first < 10) std::cerr << " ";
    }
    std::cerr << std::endl;
    for (uint64_t i = 0; i < genome_size; i++) {
        std::cerr << B[i].first.second << " ";
        if (B[i].second >= 10 && B[i].first.second < 10) std::cerr << " ";
    }
    std::cerr << std::endl;
    for (uint64_t i = 0; i < genome_size; i++) std::cerr << B[i].second << " ";
    // std::cerr << std::endl << std::endl;

    // for (uint64_t i = 0; i < genome_size; i++) {
    //     std::cerr << &buffer[B[i].second] << std::endl;
    // }
}

const std::vector<uint64_t> sa_word_size_param(
    int my_rank, int number_of_processes, int which, DataSource &data_source,
    std::vector<std::string> queries) {
    const uint64_t genome_size = data_source.getTotalGenomeSize(which),
                   my_genome_part_size = data_source.getNodeGenomeSize(which),
                   my_genome_offset = data_source.getNodeGenomeOffset(which);

    std::vector<int> send_counts(number_of_processes),
        recv_counts(number_of_processes), send_offsets(number_of_processes),
        recv_offsets(number_of_processes);
    std::string buffer(my_genome_part_size + K_VAL, 0);
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> B(
        my_genome_part_size + 1);
    std::vector<std::pair<uint64_t, uint64_t>> B_prim_source(
        my_genome_part_size),
        B_prim(my_genome_part_size);
    std::vector<uint64_t> B_plus_h(my_genome_part_size),
        B_helper(my_genome_part_size);

    data_source.getNodeGenomeValues(which, buffer.data());

    // Create B with K_VAL-mers

    // Send my K_VAL first elements to the previous node and get K_VAL next
    // from the next node. We choose K_VAL, such that it fits onto every
    // node
    if (my_rank > 0)
        MPI_Isend(buffer.c_str(), (int)K_VAL, MPI_CHAR, my_rank - 1, 0,
                  MPI_COMM_WORLD, &global_request);
    MPI_Request get_k_request;
    if (my_rank < number_of_processes - 1)
        MPI_Irecv(&buffer.data()[my_genome_part_size], (int)K_VAL, MPI_CHAR,
                  my_rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &get_k_request);

    uint64_t M = 1;
    uint64_t current_value = 0;

    for (uint64_t i = 0; i < K_VAL; i++) {
        current_value *= (1 << char_size);
        if (i < buffer.size()) {
            if (i == my_genome_part_size && my_genome_offset + i < genome_size)
                MPI_Wait(&get_k_request, nullptr);
            current_value += char_to_word(buffer[i]);
        }
        if (i > 0) M *= (1 << char_size);
    }

    for (size_t i = 0; i < my_genome_part_size; i++) {
        B[i] = std::make_pair(std::make_pair(current_value, 0),
                              my_genome_offset + i);
        current_value %= M;
        current_value *= (1 << char_size);
        const size_t j = i + K_VAL;
        if (j < buffer.size()) {
            if (j == my_genome_part_size && my_genome_offset + j < genome_size)
                MPI_Wait(&get_k_request, nullptr);
            current_value += char_to_word(buffer[j]);
        }
    }

    // Create SA
    my_sort(B);  // sorting

    bool done = rebucket_and_check_all_singleton(
        my_rank, number_of_processes, genome_size, my_genome_part_size,
        my_genome_offset, B);

    for (uint64_t h = K_VAL;; h <<= 1) {
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> to_send_to(
            number_of_processes);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            const int to = (int)whose(B[i].second);
            to_send_to[to].emplace_back(B[i].second - offset(to),
                                        B[i].first.first);
        }

        // Communicate processes sending values at desired positions
        size_t pos_in_B_prim = 0;
        send_offsets[0] = 0;
        for (int i = 0; i < number_of_processes; i++) {
            send_counts[i] = (int)to_send_to[i].size() *
                             sizeof(std::pair<uint64_t, uint64_t>);
            if (i > 0)
                send_offsets[i] = send_offsets[i - 1] + send_counts[i - 1];
            for (size_t j = 0; j < to_send_to[i].size(); j++, pos_in_B_prim++) {
                B_prim_source[pos_in_B_prim] = to_send_to[i][j];
            }
        }

        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1,
                     MPI_INT, MPI_COMM_WORLD);

        recv_offsets[0] = 0;
        for (int i = 1; i < number_of_processes; i++)
            recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];

        MPI_Alltoallv(B_prim_source.data(), send_counts.data(),
                      send_offsets.data(), MPI_BYTE, B_prim.data(),
                      recv_counts.data(), recv_offsets.data(), MPI_BYTE,
                      MPI_COMM_WORLD);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            const size_t j = B_prim[i].first;
            B[j].first.first = B_prim[i].second;
            B_helper[j] = B[j].first.first;
        }
        if (done) {
            break;
        }

        // communicate (send to two, receive from two)
        MPI_Request got_B_plus_h_request[2];
        int first_receiver_rank;
        if (my_genome_offset >= h && (first_receiver_rank = (int)whose(
                                          my_genome_offset - h)) != my_rank) {
            uint64_t first_receiver_relative_offset =
                my_genome_offset - h - offset(first_receiver_rank);
            uint64_t first_receiver_size = how_much_x_has(first_receiver_rank) -
                                           first_receiver_relative_offset;
            MPI_Isend(B_helper.data(), (int)first_receiver_size, MPI_UINT64_T,
                      first_receiver_rank, 0, MPI_COMM_WORLD, &global_request);
        }
        int second_receiver_rank;
        if (my_genome_offset + my_genome_part_size > h + 1 &&
            (second_receiver_rank = (int)whose(
                 my_genome_offset + my_genome_part_size - h - 1)) != my_rank) {
            uint64_t second_receiver_size = my_genome_offset +
                                            my_genome_part_size - h -
                                            offset(second_receiver_rank);
            uint64_t second_receiver_relative_offset =
                my_genome_part_size - second_receiver_size;
            ;
            if (second_receiver_relative_offset < my_genome_part_size)
                MPI_Isend(&B_helper.data()[second_receiver_relative_offset],
                          (int)second_receiver_size, MPI_UINT64_T,
                          second_receiver_rank, 0, MPI_COMM_WORLD, nullptr);
        }
        int first_sender_rank = -1;
        if (h >= my_genome_part_size && my_genome_offset + h < genome_size &&
            (first_sender_rank = (int)whose(my_genome_offset + h)) != my_rank) {
            uint64_t first_sender_relative_offset =
                my_genome_offset + h - offset(first_sender_rank);
            uint64_t first_sender_size = how_much_x_has(first_sender_rank) -
                                         first_sender_relative_offset;

            MPI_Irecv(B_plus_h.data(), (int)first_sender_size, MPI_UINT64_T,
                      first_sender_rank, MPI_ANY_TAG, MPI_COMM_WORLD,
                      &got_B_plus_h_request[0]);
        }
        int second_sender_rank;
        uint64_t second_sender_relative_offset = 0;
        if (my_genome_offset + my_genome_part_size + h - 1 < genome_size &&
            (second_sender_rank = (int)whose(
                 my_genome_offset + my_genome_part_size + h - 1)) != my_rank) {
            uint64_t second_sender_size = my_genome_offset +
                                          my_genome_part_size + h -
                                          offset(second_sender_rank);
            second_sender_relative_offset =
                my_genome_part_size - second_sender_size;
            MPI_Irecv(&B_plus_h.data()[second_sender_relative_offset],
                      (int)second_sender_size, MPI_UINT64_T, second_sender_rank,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &got_B_plus_h_request[1]);
        }

        if (first_sender_rank != -1)
            MPI_Wait(&got_B_plus_h_request[0], nullptr);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            if (i + h < my_genome_part_size)
                B[i].first.second = B[i + h].first.first;
            else if (my_genome_offset + i + h < genome_size) {
                if (i == second_sender_relative_offset)
                    MPI_Wait(&got_B_plus_h_request[1], nullptr);
                B[i].first.second = B_plus_h[i];
            } else
                B[i].first.second = 0;
        }
        my_sort(B);
        done = rebucket_and_check_all_singleton(
            my_rank, number_of_processes, genome_size, my_genome_part_size,
            my_genome_offset, B);
    }

    // TODO send borders
    // Answer the queries
    std::vector<uint64_t> res(queries.size());
    for (uint64_t i = 0; i < queries.size(); i++) {
        uint64_t first_occurrence, after_last_occurrence;
        // find first occurrence
        {
            uint64_t b = 0, e = B.size(), m;
            while (b < e) {
                m = (b + e) / 2;
                if (strcmp(queries[i].c_str(),
                           &buffer.c_str()[B[m].second]) >
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
                if (strcmp(queries[i].c_str(),
                           &buffer.c_str()[B[m].second]) >
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
