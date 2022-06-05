#ifndef SA_SEQ_H
#define SA_SEQ_H
#include "data_source.h"

void sa(int my_rank, int number_of_processes, uint64_t n, uint64_t m,
        DataSource &data_source, const std::string &queries_in,
        const std::string &queries_out);

#endif /* SA_SEQ_H */
