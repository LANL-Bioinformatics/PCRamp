OBJS = update.o options.o \
	seq_overlap.o sequence.o word.o optimize.o \
	background_match.o parse_fasta.o \
	optimize_pcr.o pcr_assay.o \
	sample.o \
	nuc_cruc.o nuc_cruc_santa_lucia.o \
	valid_pcr.o \
	select_words.o mpi_util.o \
	JSON.o
	
CC = mpic++
cc = mpicc

# 1) Define NDEBUG to turn off assert() checking
# 2) Unordered map and set require -std=c++0x
# 3) For reasons that I do not understand, the -mpopcnt
#    flag is needed to enable the builtin popcnt function to use
#    the hardware acceleration present in the x86 instruction set!
#    **update** It appears that -msse4.2 provides even more acceleration
#    (but I'm not sure if it is turning on -mpopcnt or not). Combining
#    -msse4.2 and -mpopcnt seems to slightly slow things down (in a single
#    benchmark experiment).
# 4) Define POPCNT to use the __builtin_popcountl command (for machines that support it). 
#    Otherwise, the code defaults to a slightly slower version.
# 5) Code now requires -DATA16 (and not -DATA32) in order to support efficient TaqMan PCR
#    assay searching.
# 6) Getting OpenMP to with on OS X is thanks to:
# 		https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
#
# 	 Please note that DYLD_LIBRARY_PATH must be set to the directory that contains
# 	 the OpenMP libomp.dylib file:
# 	 export DYLD_LIBRARY_PATH=/Users/jgans/llvm-project/build-openmp/runtime/src

#FLAGS = -O3 -Wall -std=c++0x -DATA32 -DPOPCNT -msse4.2
FLAGS = -O3 -Wall -std=c++0x -DATA16 -DPOPCNT -msse4.2

PROFILE = #-pg
OPENMP = -Xpreprocessor -fopenmp
	
INC = -I. -I/Users/jgans/llvm-project/build-openmp/runtime/src
LIBS = -lm -lz -L/Users/jgans/llvm-project/build-openmp/runtime/src -lomp

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<
 
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<

all: pcramp

pcramp: $(OBJS) main.o 
	$(CC) $(PROFILE) -o pcramp $(OBJS) $(OPENMP) main.o $(LIBS)
 	
clean:
	-rm -f *.o


