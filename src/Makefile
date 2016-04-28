
# list of executable files to produce

DRAW = draw
THEORETICAL_TEST = theoretical
EMPIRICAL_TEST = empirical
CONTROL_TEST = control

# .o files necessary to build the executables
OBJ_DRAW = AsymptoticSpeed.o Text.o Bernoulli.o Utils.o Lattice.o  PowerSet.o MatchingMachine.o PatternMatchingAlgorithms.o Explore.o Draw.o
OBJ_THEORETICAL_TEST = AsymptoticSpeed.o Text.o Bernoulli.o Utils.o Lattice.o TableUtils.o PowerSet.o MatchingMachine.o PatternMatchingAlgorithms.o Explore.o TheoreticalTest.o 
OBJ_EMPIRICAL_TEST = AsymptoticSpeed.o Text.o Bernoulli.o Utils.o Lattice.o TableUtils.o  PowerSet.o MatchingMachine.o PatternMatchingAlgorithms.o Explore.o ReferenceAlgorithmsCount.o EmpiricalTest.o 
OBJ_CONTROL_TEST = AsymptoticSpeed.o Text.o Bernoulli.o Utils.o Lattice.o TableUtils.o PowerSet.o MatchingMachine.o PatternMatchingAlgorithms.o ReferenceAlgorithmsCount.o Explore.o ControlTest.o 
########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
#CFLAGS = -Wall -Wno-char-subscripts  -std=c99 -Werror -pedantic -O3-march=native
#CFLAGS = -Wall -Wno-char-subscripts  -std=c99 -Wall -pedantic -O3 -march=native 
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -O3 -g -march=native



# linker and its options
LD = $(CC)
LDFLAGS = -lsuperlu -lm -O3


############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DRAW) $(THEORETICAL_TEST) $(EMPIRICAL_TEST) $(EMPIRICAL_TEST_BIN) $(CONTROL_TEST)

# build the executable
$(DRAW): $(OBJ_DRAW)
	$(LD) $^ -o $@ $(LDFLAGS)
# build the executable
$(THEORETICAL_TEST): $(OBJ_THEORETICAL_TEST)
	$(LD) $^ -o $@ $(LDFLAGS)
# build the executable
$(EMPIRICAL_TEST): $(OBJ_EMPIRICAL_TEST)
	$(LD) $^ -o $@ $(LDFLAGS)
# build the executable
$(CONTROL_TEST): $(OBJ_CONTROL_TEST)
	$(LD) $^ -o $@ $(LDFLAGS)

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release

-include Makefile.d
