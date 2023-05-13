# adapted from https://stackoverflow.com/questions/2394609/makefile-header-dependencies sophie and maximumfps
ifeq ($(strip $(SYLDEV)),1)
$(info Using sylvan-develop)
GXX = g++ -DSYLDEV
else
GXX = g++
SYLDEV =
$(info Using plain sylvan (not develop))
endif
CXXFLAGS = #-Wfatal-errors -Wall -Wextra -Wpedantic -Wconversion -Wshadow
ASSEMBLY_OPTIONS = -g -std=c++2a -Wno-register -O0
LINK_OPTIONS = -lgmpxx -lgmp -lpthread -static # order matters

################################################################################

CXXOPTS_TARGET = libraries/cxxopts/include/cxxopts.hpp

CUDD_DIR = libraries/cudd
CUDD_LIB_DIR = $(CUDD_DIR)/cudd/.libs
CUDD_TARGET = $(CUDD_LIB_DIR)/libcudd.a
CUDD_INCLUSIONS = -I$(CUDD_DIR) -I$(CUDD_DIR)/cudd -I$(CUDD_DIR)/epd -I$(CUDD_DIR)/mtr -I$(CUDD_DIR)/st
CUDD_LINKS = -L$(CUDD_LIB_DIR) -lcudd

ifdef SYLDEV
SYLVAN_DIR = libraries/sylvan-develop/
else
SYLVAN_DIR = libraries/sylvan/
endif
SYLVAN_BUILD_DIR = $(SYLVAN_DIR)/build
SYLVAN_LIB_DIR = $(SYLVAN_BUILD_DIR)/src
SYLVAN_TARGET = $(SYLVAN_LIB_DIR)/libsylvan.a
SYLVAN_INCLUSIONS = -I$(SYLVAN_DIR)/src -I$(SYLVAN_BUILD_DIR)/_deps/lace-src/src/ -I$(SYLVAN_BUILD_DIR)/_deps/lace-build
SYLVAN_LINKS = -L$(SYLVAN_LIB_DIR) -lsylvan

LACE_LINKS = -L$(SYLVAN_BUILD_DIR)/_deps/lace-build -llace

CMSAT_DIR = libraries/cryptominisat
CMSAT_BUILD_DIR = $(CMSAT_DIR)/build
CMSAT_LIB_DIR = $(CMSAT_BUILD_DIR)/lib
CMSAT_TARGET = $(CMSAT_LIB_DIR)/libcryptominisat5.a
CMSAT_INCLUSIONS = -I$(CMSAT_BUILD_DIR)/include
CMSAT_LINKS = -L$(CMSAT_LIB_DIR) -lcryptominisat5

COLAMD_DIR = libraries/colamd
COLAMD_LINKS = -L$(COLAMD_DIR) -lcolamd
################################################################################


# Final binary
BIN = dmc
# Put all auto generated stuff to this build dir.
BUILD_DIR = ./build

# List of all .cpp source files.
CPP = $(wildcard src/*.cpp)

# All .o files go to build dir.
OBJ = $(CPP:%.cpp=$(BUILD_DIR)/%.o)
# Gcc/Clang will create these .d files containing dependencies.
DEP = $(OBJ:%.o=%.d)

# Default target named after the binary.
$(BIN) : $(BUILD_DIR)/$(BIN)

# Actual target of the binary - depends on all .o files.
$(BUILD_DIR)/$(BIN) : $(OBJ)
	# Create build directories - same structure as sources.
	mkdir -p $(@D)
	# Just link all the object files.
	$(GXX) $(ASSEMBLY_OPTIONS) $(CXXFLAGS) $^ -o $@ $(CUDD_LINKS) $(SYLVAN_LINKS) $(LACE_LINKS) $(CMSAT_LINKS) $(COLAMD_LINKS) $(LINK_OPTIONS)

# Include all .d files
-include $(DEP)

# Build target for every single object file.
# The potential dependency on header files is covered
# by calling `-include $(DEP)`.
$(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(@D)
	# The -MMD flags additionaly creates a .d file with
	# the same name as the .o file.
	$(GXX) -MMD $(ASSEMBLY_OPTIONS) $(CXXFLAGS) $(CUDD_INCLUSIONS) $(SYLVAN_INCLUSIONS) $(CMSAT_INCLUSIONS) -c $< -o $@

.PHONY : clean
clean :
	# This should remove all generated files.
	-rm $(BUILD_DIR)/$(BIN) $(OBJ) $(DEP)
