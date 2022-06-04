.PHONY: all clean distclean test-clean test

all: genome_index

CC := mpic++
WARNING_OPTIONS := -Wextra -Wconversion -Wall -Wpedantic
WEAK_COMPILER_OPTIONS := -O3 -fopenmp
COMPILER_OPTIONS := $(WARNING_OPTIONS) $(WEAK_COMPILER_OPTIONS)
WEAK_LINKER_OPTIONS := $(WEAK_COMPILER_OPTIONS) -c
LINKER_OPTIONS := $(COMPILER_OPTIONS) -c

genome_index: src/main.cpp data_source.o
	$(CC) $(COMPILER_OPTIONS) -o $@ $^

data_source.o : src/data_source.cpp src/data_source.h Makefile
	$(CC) $(WEAK_LINKER_OPTIONS) $<

%.o : %.cpp src/data_source.h Makefile
	$(CC) $(LINKER_OPTIONS) $<

clean:
	rm -f src/*.o *.o

distclean: clean
	rm -f genome_index
