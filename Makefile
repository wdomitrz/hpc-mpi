.PHONY: all clean test-clean test

all: genome_index

CC := CC
WARNING_OPTIONS := -Wextra -Wconversion -Wall -Wpedantic -Wno-sign-conversion
WEAK_LINKER_OPTIONS := -O3 -fopenmp -std=c++17
LINKER_OPTIONS := $(WARNING_OPTIONS) $(WEAK_LINKER_OPTIONS)
WEAK_COMPILTER_OPTIONS := $(WEAK_LINKER_OPTIONS) -c
COMPILTER_OPTIONS := $(LINKER_OPTIONS) -c

genome_index: genome_index_v1
	cp genome_index_v1 genome_index

genome_index_v2: src/main.cpp data_source.o sa_v2.o src/params.hpp data_source.o src/sa.hpp src/data_source.h
	$(CC) $(LINKER_OPTIONS) -o $@ $< data_source.o sa_v2.o

genome_index_v1: src/main.cpp data_source.o sa.o src/params.hpp data_source.o src/sa.hpp src/data_source.h
	$(CC) $(LINKER_OPTIONS) -o $@ $< data_source.o sa.o

genome_index_seq: src/main.cpp src/params.hpp data_source.o sa_seq.o src/sa.hpp src/data_source.h
	$(CC) $(LINKER_OPTIONS) -o $@ $< data_source.o sa_seq.o

data_source.o : src/data_source.cpp src/data_source.h Makefile
	$(CC) $(WEAK_COMPILTER_OPTIONS) $<

sa_v2.o: src/sa_v2.cpp src/sa.hpp src/data_source.h Makefile
	$(CC) $(COMPILTER_OPTIONS) $<

sa.o: src/sa.cpp src/sa.hpp src/data_source.h Makefile
	$(CC) $(COMPILTER_OPTIONS) $<

sa_seq.o: src/sa_seq.cpp src/sa.hpp src/data_source.h Makefile
	$(CC) $(COMPILTER_OPTIONS) $<

clean:
	rm -f src/*.o *.o genome_index genome_index_seq genome_index_v1 genome_index_v2

report.pdf: README.md
	pandoc $^ -o $@
