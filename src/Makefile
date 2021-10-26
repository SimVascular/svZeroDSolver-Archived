# g++ -o unit_tests unit_tests.cpp blocks.cpp

# reference: https://web.stanford.edu/class/archive/cs/cs107/cs107.1194/resources/make

# The default C++ compiler
CC = g++

# The CFLAGS variable sets compile flags for g++:
#  -std=gnu99  use the GNU99 standard language definition
CFLAGS =

unit_tests: unit_tests.cpp blocks.cpp
	$(CC) $(CFLAGS) -o unit_tests unit_tests.cpp blocks.cpp

.PHONY: clean

clean:
	rm -f unit_tests *.o
