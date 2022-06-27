#include "sa.hpp"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include "data_source.h"
using std::vector, std::pair;

#define K_VAL (std::min(k_max, genome_size / number_of_processes))

#define whose(x) (whose_param(x, number_of_processes, genome_size))
#define how_much_x_has(rank) \
    (how_much_node_has(rank, number_of_processes, genome_size))
#define offset(rank) \
    (get_node_genome_offset(rank, number_of_processes, genome_size))
#define my_sort_full(B)                                        \
    (my_sort_params(my_rank, number_of_processes, genome_size, \
                    my_genome_part_size, B, 0, genome_size))
#define printB(x, y) printB_fun(B, buffer.c_str(), my_genome_part_size, x, y)
#define ok() \
    { std::cerr << "ok:\t" << my_rank << "\t" << __LINE__ << std::endl; }
#define assertm(exp, msg) assert(((void)msg, exp))

// std::random_device random_device;
std::mt19937 random_generator(1);

uint64_t max_query_size = 0;
const int ROOT = 0;

MPI_Request global_request;
MPI_Status global_status;
MPI_Datatype my_MPI_UINT64_Triplet;
MPI_Datatype my_MPI_UINT64_Pair;

const size_t char_size = 3;  // There are 4 characters -- A, C, T, G and one
                             // special character -- the end of the word

static inline uint64_t whose_param(uint64_t x, uint64_t number_of_processes,
                                   uint64_t total_size) {
    const uint64_t M = total_size % number_of_processes,
                   step = total_size / number_of_processes;
    const uint64_t threshold = M * (step + 1);

    if (x < threshold) return x / (step + 1);

    return (x - M) / (step);
}

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

void printBs(
    int number_of_processes, uint64_t my_genome_part_size, int my_rank,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B) {
    for (int r = 0; r < number_of_processes; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (r == my_rank) {
            for (uint64_t i = 0; i < my_genome_part_size; i++) {
                std::cerr << B[i].first.first << " ";
            }
            std::cerr << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (r == my_rank) {
            for (uint64_t i = 0; i < my_genome_part_size; i++) {
                std::cerr << B[i].first.second << " ";
            }
            std::cerr << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (r == my_rank) {
            for (uint64_t i = 0; i < my_genome_part_size; i++) {
                std::cerr << B[i].second << " ";
            }
            if (my_rank == number_of_processes - 1) std::cerr << std::endl;
            std::cerr << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

inline static uint64_t how_much_node_has(int rank, int nprocs,
                                         uint64_t totalSize) {
    uint64_t nodeSize = (totalSize + nprocs - rank - 1) / nprocs;
    return nodeSize;
}

inline static uint64_t get_node_genome_offset(int rank, int nprocs,
                                              uint64_t totalSize) {
    MPI_Offset offset =
        static_cast<uint64_t>(rank) * (totalSize / (uint64_t)nprocs) +
        std::min((uint64_t)totalSize % nprocs, (uint64_t)rank);
    return offset;
}

inline bool rebucket_and_check_all_singleton(
    int my_rank, int number_of_processes, const uint64_t,
    const uint64_t my_genome_part_size, const uint64_t,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B) {
    std::pair<std::pair<uint64_t, uint64_t>, uint64_t> first_val = B[0];
    bool my_res = true;
    std::vector<uint64_t> my_partial_results(my_genome_part_size);

    if (my_rank > 0)
        MPI_Isend(&first_val, 1, my_MPI_UINT64_Triplet, my_rank - 1, 0,
                  MPI_COMM_WORLD,
                  &global_request);  // Instead of creating a custom MPI type, I
                                     // use MPI_UINT64_T and send data as bytes.

    MPI_Request get_next_one_request;
    if (my_rank < number_of_processes - 1)
        MPI_Irecv(
            &B[my_genome_part_size], 1, my_MPI_UINT64_Triplet, my_rank + 1,
            MPI_ANY_TAG, MPI_COMM_WORLD,
            &get_next_one_request);  // Instead of creating a custom MPI type, I
                                     // use MPI_UINT64_T and send data as bytes.

    uint64_t my_count = 0;
    std::pair<uint64_t, uint64_t> prev_val = B[0].first;
    B[0].first = std::make_pair(1, 0);
    for (size_t i = 1; i < my_genome_part_size; i++) {
        if (prev_val == B[i].first) {
            B[i].first = B[i - 1].first;
            my_res = false;
        } else {
            prev_val = B[i].first;
            B[i].first = B[i - 1].first;
            B[i].first.first++;
            my_count++;
        }
    }

    if (my_rank + 1 != number_of_processes) {
        MPI_Wait(&get_next_one_request, &global_status);
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

inline void my_sort_params_old(
    const int my_rank, const int number_of_processes,
    const uint64_t genome_size, const uint64_t my_genome_part_size,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const uint64_t, const uint64_t) {
    if (my_rank == ROOT) {
        std::vector<int> recv_size(number_of_processes);
        std::vector<int> recv_offset(number_of_processes);
        recv_size[0] = static_cast<int>(how_much_x_has(0));
        recv_offset[0] = 0;
        for (int i = 1; i < number_of_processes; i++) {
            recv_offset[i] = recv_offset[i - 1] + recv_size[i - 1];
            recv_size[i] = static_cast<int>(how_much_x_has(i));
        }

        std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> B_all(
            genome_size);

        MPI_Gatherv(B.data(), static_cast<int>(my_genome_part_size),
                    my_MPI_UINT64_Triplet, B_all.data(), recv_size.data(),
                    recv_offset.data(), my_MPI_UINT64_Triplet, 0,
                    MPI_COMM_WORLD);
        std::sort(B_all.begin(), B_all.end());
        MPI_Scatterv(B_all.data(), recv_size.data(), recv_offset.data(),
                     my_MPI_UINT64_Triplet, B.data(),
                     static_cast<int>(my_genome_part_size),
                     my_MPI_UINT64_Triplet, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(B.data(), static_cast<int>(my_genome_part_size),
                    my_MPI_UINT64_Triplet, nullptr, nullptr, nullptr,
                    my_MPI_UINT64_Triplet, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, my_MPI_UINT64_Triplet, B.data(),
                     static_cast<int>(my_genome_part_size),
                     my_MPI_UINT64_Triplet, 0, MPI_COMM_WORLD);
    }
}

inline void my_sort_params(
    const int my_rank, const int number_of_processes,
    const uint64_t genome_size, const uint64_t my_genome_part_size,
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const uint64_t begin, const uint64_t end) {
    if (end <= begin) return;
    assert(begin == 0 && end == genome_size);

    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> pivots(
        number_of_processes);
    std::vector<int> send_offsets(number_of_processes),
        send_counts(number_of_processes, 0), recv_offsets(number_of_processes),
        recv_counts(number_of_processes);
    std::uniform_int_distribution<uint64_t> uniform_distribution(
        0, my_genome_part_size - 1);

    pivots[my_rank] = B[uniform_distribution(random_generator)];
    MPI_Request pivots_request;
    MPI_Iallgather(MPI_IN_PLACE, 1, my_MPI_UINT64_Triplet, pivots.data(), 1,
                   my_MPI_UINT64_Triplet, MPI_COMM_WORLD, &pivots_request);
    std::sort(B.data(), &B.data()[my_genome_part_size]);
    MPI_Wait(&pivots_request, &global_status);
    std::sort(pivots.begin(), pivots.end() - 1);

    int current_receiver = 0;
    for (uint64_t i = 0; i < my_genome_part_size; i++) {
        while (current_receiver + 1 < number_of_processes &&
               pivots[current_receiver] < B[i])
            current_receiver++;
        send_counts[current_receiver]++;
    }

    send_offsets[0] = 0;
    for (int i = 1; i < number_of_processes; i++) {
        send_offsets[i] = send_offsets[i - 1] + send_counts[i - 1];
    }

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);

    recv_offsets[0] = 0;
    for (int i = 1; i < number_of_processes; i++) {
        recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
    }

    uint64_t my_temp_size = recv_offsets[number_of_processes - 1] +
                            recv_counts[number_of_processes - 1],
             my_temp_offset;

    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> temp_B(
        my_temp_size + 1);

    MPI_Alltoallv(B.data(), send_counts.data(), send_offsets.data(),
                  my_MPI_UINT64_Triplet, (uint64_t *)temp_B.data(),
                  recv_counts.data(), recv_offsets.data(),
                  my_MPI_UINT64_Triplet, MPI_COMM_WORLD);

    MPI_Request temp_offset_request;

    MPI_Iexscan(&my_temp_size, &my_temp_offset, 1, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD, &temp_offset_request);

    std::sort(temp_B.begin(), temp_B.end() - 1);
    std::fill(send_counts.begin(), send_counts.end(), 0);
    std::fill(send_offsets.begin(), send_offsets.end(), 0);

    MPI_Wait(&temp_offset_request, &global_status);
    if (my_rank == 0) {
        my_temp_offset = 0;
    }

    send_offsets[0] = 0;
    for (uint64_t i = (whose(my_temp_offset));
         my_temp_offset < genome_size && my_temp_size > 0 &&
         i <= whose(my_temp_offset + my_temp_size - 1);
         i++) {
        send_counts[i] = static_cast<int>(std::min(
            my_temp_size, how_much_x_has(static_cast<int>(i)) +
                              offset(static_cast<int>(i)) - my_temp_offset));
        my_temp_size -= send_counts[i];
        my_temp_offset += send_counts[i];
        if (i > 0) send_offsets[i] = send_offsets[i - 1] + send_counts[i - 1];
    }

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);

    recv_offsets[0] = 0;
    for (int i = 1; i < number_of_processes; i++) {
        recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
    }

    MPI_Alltoallv(temp_B.data(), send_counts.data(), send_offsets.data(),
                  my_MPI_UINT64_Triplet, B.data(), recv_counts.data(),
                  recv_offsets.data(), my_MPI_UINT64_Triplet, MPI_COMM_WORLD);
}

inline void printB_fun(
    const std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const char *buffer, uint64_t genome_size, const uint64_t x1,
    const uint64_t x2) {
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
    std::cerr << std::endl << std::endl;

    for (uint64_t i = 0; i < genome_size; i++) {
        if (i == x1 || i == x2) std::cerr << "\t";
        std::cerr << &buffer[B[i].second] << std::endl;
    }
}

uint64_t my_lower_bound_old(  // TODO: compare with current implementation
    const int my_rank, const int number_of_processes,
    const uint64_t genome_size, const uint64_t my_genome_part_size,
    const uint64_t my_genome_offset, const std::string &query,
    const std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const std::string &buffer) {
    std::vector<uint64_t> last_values(number_of_processes);
    std::vector<int8_t> is_greater_or_equal(number_of_processes, false);
    MPI_Allgather(&B[my_genome_part_size - 1].second, 1, MPI_UINT64_T,
                  last_values.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

    for (int i = 0; i < number_of_processes; i++) {
        if (my_genome_offset <= last_values[i] &&
            last_values[i] < my_genome_offset + my_genome_part_size)
            is_greater_or_equal[i] =
                std::strcmp(
                    query.c_str(),
                    &buffer.c_str()[last_values[i] - my_genome_offset]) <= 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, is_greater_or_equal.data(), number_of_processes,
                  MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    int found_rank = 0;
    while (found_rank < number_of_processes &&
           !is_greater_or_equal[found_rank]) {
        found_rank++;
    }
    if (found_rank == number_of_processes) return genome_size;

    uint64_t size_of_found = how_much_x_has(found_rank);
    std::vector<uint64_t> SA(size_of_found);
    std::vector<int8_t> is_here(size_of_found, false);
    if (my_rank == found_rank) {
        for (uint64_t i = 0; i < SA.size(); i++) {
            SA[i] = B[i].second;
        }
        MPI_Bcast(SA.data(), static_cast<int>(SA.size()), MPI_UINT64_T,
                  found_rank, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(SA.data(), static_cast<int>(SA.size()), MPI_UINT64_T,
                  found_rank, MPI_COMM_WORLD);
    }

    for (size_t i = 0; i < SA.size(); i++) {
        if (my_genome_offset <= SA[i] &&
            SA[i] < my_genome_offset + my_genome_part_size)
            is_here[i] = strcmp(query.c_str(),
                                &buffer.c_str()[SA[i] - my_genome_offset]) <= 0;
    }

    MPI_Allreduce(MPI_IN_PLACE, is_here.data(),
                  static_cast<int>(is_here.size()), MPI_C_BOOL, MPI_LOR,
                  MPI_COMM_WORLD);

    uint64_t pos = 0;
    while (!is_here[pos])  // We know that is_here[is_here.size()-1] is true.
        pos++;

    return offset(found_rank) + pos;
}

inline uint64_t my_lower_bound(
    const int my_rank, const int number_of_processes,
    const uint64_t genome_size, const uint64_t, const uint64_t my_genome_offset,
    const std::string &query,
    const std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> &B,
    const std::string &buffer) {
    uint64_t b = 0, e = genome_size, m;
    uint64_t pos;
    int who_now;
    bool comp;
    while (b < e) {
        m = (b + e) / 2;

        if ((who_now = static_cast<int>(whose(m))) == my_rank)
            pos = B[m - my_genome_offset].second;
        MPI_Bcast(&pos, 1, MPI_UINT64_T, who_now, MPI_COMM_WORLD);

        if ((who_now = static_cast<int>(whose(pos))) == my_rank)
            comp = strcmp(query.c_str(),
                          &buffer.c_str()[pos - my_genome_offset]) <= 0;
        MPI_Bcast(&comp, 1, MPI_C_BOOL, who_now, MPI_COMM_WORLD);

        if (!comp) {
            b = m + 1;
        } else {
            e = m;
        }
    }
    return b;
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
    const uint64_t extension_size = std::max(max_query_size, K_VAL);
    std::string buffer(my_genome_part_size + extension_size, 0);
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
    // from the next node.
    uint64_t to_get = extension_size;
    recv_counts[my_rank] = static_cast<int>(my_genome_part_size);
    for (int i = my_rank + 1; i < number_of_processes; i++) {
        recv_counts[i] = static_cast<int>(std::min(to_get, how_much_x_has(i)));
        to_get -= recv_counts[i];
        recv_offsets[i] =
            recv_offsets[i - 1] +
            recv_counts[i - 1];  // I can use i - 1, as i = my_rank
                                 // + 1, so always i > 0
    }
    MPI_Alltoall(recv_counts.data(), 1, MPI_INT, send_counts.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);
    MPI_Alltoallv(buffer.data(), send_counts.data(), send_offsets.data(),
                  MPI_CHAR, buffer.data(), recv_counts.data(),
                  recv_offsets.data(), MPI_CHAR, MPI_COMM_WORLD);

    uint64_t M = 1;
    uint64_t current_value = 0;

    for (uint64_t i = 0; i < K_VAL; i++) {
        current_value *= (1 << char_size);
        if (i < buffer.size()) {
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
            current_value += char_to_word(buffer[j]);
        }
    }

    // Create SA

    my_sort_full(B);  // sorting
    bool done = rebucket_and_check_all_singleton(
        my_rank, number_of_processes, genome_size, my_genome_part_size,
        my_genome_offset, B);

    for (uint64_t h = K_VAL; h < genome_size; h *= 2) {
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> to_send_to(
            number_of_processes);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            const uint64_t to = whose(B[i].second);
            to_send_to[to].emplace_back(B[i].second, B[i].first.first);
        }

        // Communicate processes sending values at desired positions
        size_t pos_in_B_prim = 0;
        send_offsets[0] = 0;
        for (int i = 0; i < number_of_processes; i++) {
            send_counts[i] = static_cast<int>(to_send_to[i].size());
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
                      send_offsets.data(), my_MPI_UINT64_Pair, B_prim.data(),
                      recv_counts.data(), recv_offsets.data(),
                      my_MPI_UINT64_Pair, MPI_COMM_WORLD);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            const size_t j = B_prim[i].first - my_genome_offset;
            B[j].first.first = B_prim[i].second;
        }
        if (done) {
            break;
        }

        for (uint64_t i = 0; i < my_genome_part_size; i++)
            B_helper[i] = B[i].first.first;
        // communicate (send to two, receive from two)
        MPI_Request got_B_plus_h_request[2];
        int first_receiver_rank;
        if (my_genome_offset >= h &&
            (first_receiver_rank =
                 static_cast<int>(whose(my_genome_offset - h))) != my_rank) {
            uint64_t first_receiver_relative_offset =
                my_genome_offset - h - offset(first_receiver_rank);
            uint64_t first_receiver_size = how_much_x_has(first_receiver_rank) -
                                           first_receiver_relative_offset;
            MPI_Isend(B_helper.data(), static_cast<int>(first_receiver_size),
                      MPI_UINT64_T, first_receiver_rank, 0, MPI_COMM_WORLD,
                      &global_request);
        }
        int second_receiver_rank;
        if (my_genome_offset + my_genome_part_size > h + 1 &&
            (second_receiver_rank = static_cast<int>(whose(
                 my_genome_offset + my_genome_part_size - h - 1))) != my_rank) {
            uint64_t second_receiver_size = my_genome_offset +
                                            my_genome_part_size - h -
                                            offset(second_receiver_rank);
            uint64_t second_receiver_relative_offset =
                my_genome_part_size - second_receiver_size;
            ;
            if (second_receiver_relative_offset < my_genome_part_size)
                MPI_Isend(&B_helper.data()[second_receiver_relative_offset],
                          static_cast<int>(second_receiver_size), MPI_UINT64_T,
                          second_receiver_rank, 0, MPI_COMM_WORLD,
                          &global_request);
        }
        bool first_recv = false;
        int first_sender_rank;
        if (h >= my_genome_part_size && my_genome_offset + h < genome_size &&
            (first_sender_rank =
                 static_cast<int>(whose(my_genome_offset + h))) != my_rank) {
            uint64_t first_sender_relative_offset =
                my_genome_offset + h - offset(first_sender_rank);
            uint64_t first_sender_size = how_much_x_has(first_sender_rank) -
                                         first_sender_relative_offset;

            first_recv = true;
            MPI_Irecv(B_plus_h.data(), static_cast<int>(first_sender_size),
                      MPI_UINT64_T, first_sender_rank, MPI_ANY_TAG,
                      MPI_COMM_WORLD, &got_B_plus_h_request[0]);
        }
        int second_sender_rank;
        uint64_t second_sender_relative_offset = 0;
        if (my_genome_offset + my_genome_part_size + h - 1 < genome_size &&
            (second_sender_rank = static_cast<int>(whose(
                 my_genome_offset + my_genome_part_size + h - 1))) != my_rank) {
            uint64_t second_sender_size = my_genome_offset +
                                          my_genome_part_size + h -
                                          offset(second_sender_rank);
            uint64_t second_sender_relative_offset =
                my_genome_part_size - second_sender_size;
            MPI_Irecv(&B_plus_h.data()[second_sender_relative_offset],
                      static_cast<int>(second_sender_size), MPI_UINT64_T,
                      second_sender_rank, MPI_ANY_TAG, MPI_COMM_WORLD,
                      &got_B_plus_h_request[1]);
        }

        if (first_recv) MPI_Wait(&got_B_plus_h_request[0], &global_status);
        for (uint64_t i = 0; i < my_genome_part_size; i++) {
            B[i].second = my_genome_offset + i;
            if (i + h < my_genome_part_size)
                B[i].first.second = B[i + h].first.first;
            else if (my_genome_offset + i + h < genome_size) {
                if (i == second_sender_relative_offset)
                    MPI_Wait(&got_B_plus_h_request[1], &global_status);
                B[i].first.second = B_plus_h[i];
            } else
                B[i].first.second = 0;
        }
        my_sort_full(B);
        done = rebucket_and_check_all_singleton(
            my_rank, number_of_processes, genome_size, my_genome_part_size,
            my_genome_offset, B);
    }

    // Answer the queries
    std::vector<uint64_t> res(queries.size());

    for (uint64_t i = 0; i < queries.size(); i++) {
        // find first occurrence
        const uint64_t first_occurrence = my_lower_bound(
            my_rank, number_of_processes, genome_size, my_genome_part_size,
            my_genome_offset, queries[i], B, buffer);

        // find position after last occurrence
        queries[i][queries[i].size() - 1]++;
        const uint64_t after_last_occurrence = my_lower_bound(
            my_rank, number_of_processes, genome_size, my_genome_part_size,
            my_genome_offset, queries[i], B, buffer);
        queries[i][queries[i].size() - 1]--;

        res[i] = after_last_occurrence - first_occurrence;
    }
    return res;
}

void sa(int my_rank, int number_of_processes, uint64_t n, uint64_t m,
        DataSource &data_source, const std::string &queries_in,
        const std::string &queries_out) {
    MPI_Type_contiguous(2, MPI_UINT64_T, &my_MPI_UINT64_Pair);
    MPI_Type_contiguous(3, MPI_UINT64_T, &my_MPI_UINT64_Triplet);
    MPI_Type_commit(&my_MPI_UINT64_Pair);
    MPI_Type_commit(&my_MPI_UINT64_Triplet);
    // Read queries
    std::vector<std::string> queries(m);
    {
        std::ifstream queries_in_file(queries_in);
        for (uint64_t i = 0; i < m; i++) queries_in_file >> queries[i];
    }

    max_query_size = 0;
    for (size_t i = 0; i < queries.size(); i++)
        max_query_size = std::max(max_query_size, queries[i].size());

    // Compute SA and answer the queries
    std::vector<std::vector<uint64_t>> res(n);
    for (uint64_t i = 0; i < n; i++) {
        res[i] = sa_word_size_param(my_rank, number_of_processes,
                                    static_cast<int>(i), data_source, queries);
    }

    // Write the results to a file
    if (my_rank == ROOT) {
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
