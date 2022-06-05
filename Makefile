.PHONY: all clean test-clean test

all: genome_index

CC := mpic++
WARNING_OPTIONS := -Wextra -Wconversion -Wall -Wpedantic -Wno-cast-function-type
WEAK_LINKER_OPTIONS := -O3 -fopenmp -std=c++17
LINKER_OPTIONS := $(WARNING_OPTIONS) $(WEAK_LINKER_OPTIONS)
WEAK_COMPILTER_OPTIONS := $(WEAK_LINKER_OPTIONS) -c
COMPILTER_OPTIONS := $(LINKER_OPTIONS) -c

genome_index: src/main.cpp src/params.hpp data_source.o src/sa.cpp src/sa.hpp src/data_source.h
	$(CC) $(LINKER_OPTIONS) -o $@ $^

genome_index_seq: src/main.cpp src/params.hpp data_source.o src/sa_seq.cpp src/sa.hpp src/data_source.h
	$(CC) $(LINKER_OPTIONS) -o $@ $^

data_source.o : src/data_source.cpp src/data_source.h Makefile
	$(CC) $(WEAK_COMPILTER_OPTIONS) $<

clean:
	rm -f src/*.o *.o genome_index genome_index_seq

report.pdf: report.md
	pandoc $^ -o $@
