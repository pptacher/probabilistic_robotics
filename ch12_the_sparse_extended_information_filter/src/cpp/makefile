help:
	@echo
	@echo "make lib      <-- Compile lib."
	@echo "make build    <-- Compile all."
	@echo "make clean    <-- Clean up."
	@echo "make install  <-- Compile and install"
	@echo "make dist     <-- Create .tar.gz project file"
	@echo

PROJECT_NAME = seif
VERSION      = 001

ARCHIVE_DIR_NAME = $(PROJECT_NAME)-$(VERSION)
ARCHIVE_NAME     = $(ARCHIVE_DIR_NAME).tar.gz

FILES_TO_BE_ARCHIVED = *.h *.cc Makefile

INCLUDE_PATH = -I. -I/opt/intel/compilers_and_libraries/linux/mkl/include/
#problem with linking using gcc
COMPILO     =  g++
CXXFLAGS      =  $(RELEASEFLAGS) -m64 -Wall -Wextra -Wshadow -Werror -pedantic   $(INCLUDE_PATH) -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused
CXXFLAGS_SHARED  = $(CXXFLAGS) -fPIC -shared
LIBS        = -L. -L/usr/lib/x86_64-linux-gnu/ -lseif -larmadillo -lhdf5_cpp

INSTALL_PREFIX      = ./LocalInstall

INSTALL_BIN_DIR     = $(INSTALL_PREFIX)/bin
INSTALL_LIB_DIR     = $(INSTALL_PREFIX)/lib
INSTALL_INCLUDE_DIR = $(INSTALL_PREFIX)/include/seif

INSTALL_DIRS =					\
	$(INSTALL_BIN_DIR)			\
	$(INSTALL_LIB_DIR)			\
	$(INSTALL_INCLUDE_DIR)

##########################################################
#patch for gcc: pre process to convert UTF-8 encoding to Unicode code point.
#ref: https://gcc.gnu.org/wiki/FAQ#What_is_the_status_of_adding_the_UTF-8_support_for_identifier_names_in_GCC.3F
PROCESS_UTF8 = ./to_UCN.sh

#SHELL := /bin/sh
#CC := g++
#CXX := g++

export LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64
#export MKL_VERBOSE=1

COMPILER_OPTIONS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic -Iinclude -Wno-unused-parameter -Wno-unused-variable -Wno-deprecated-declarations
#INCLUDE_PATH = -I.

CXX_DEFINES = -DSHADERS_DIR=\"/Users/TACHER/sparklepp/src/shaders\" -DUSE_GLEW

#CFLAGS := -std=c11 $(COMPILER_OPTIONS)
#CXXFLAGS := -std=c++14 -Weffc++ $(COMPILER_OPTIONS) $(INCLUDE_PATH)

DEBUGFLAGS := -g -O0 -D _DEBUG
RELEASEFLAGS := -O3 -D NDEBUG
##########################################################

# Fabrication des .o (lib)

seif.o : seif.cc
			$(PROCESS_UTF8) seif.cc
			$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/seif.cc

motion.o : motion.cc
			$(PROCESS_UTF8) motion.cc
			$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/motion.cc

estimate.o : estimate.cc
				$(PROCESS_UTF8) estimate.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/estimate.cc

correspondence.o : correspondence.cc
				$(PROCESS_UTF8) correspondence.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/correspondence.cc

markov_blanket.o : markov_blanket.cc
				$(PROCESS_UTF8) markov_blanket.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/markov_blanket.cc

measurement.o : measurement.cc
				$(PROCESS_UTF8) measurement.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/measurement.cc

sparsification.o : sparsification.cc
				$(PROCESS_UTF8) sparsification.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/sparsification.cc

linear_model.o : linear_model.cc
				$(PROCESS_UTF8) linear_model.cc
				$(COMPILO)  -c -std=c++14  $(CXXFLAGS_SHARED) /tmp/linear_model.cc

# Fabrication des .o (hors lib)

main.o : main.cc
			$(COMPILO) -c -std=c++14 $(CXXFLAGS) main.cc

# Fabrication de la lib

libseif.so : seif.o motion.o estimate.o correspondence.o measurement.o sparsification.o linear_model.o markov_blanket.o
	$(COMPILO) -o libseif.so -shared -larmadillo seif.o motion.o estimate.o correspondence.o measurement.o sparsification.o linear_model.o markov_blanket.o

# Fabrication de l'exécutable

seif.bin : main.o libseif.so
		$(COMPILO) -o   seif.bin main.o  $(LIBS) -lpthread -lmkl_core -lmkl_def -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl
#-DARMA_DONT_USE_WRAPPER -lmkl_intel_thread -lmkl_intel_ilp64
phony :
	echo $(LIBRARY_PATH)

# Les cibles de l'aide

lib : libseif.so

build : seif.bin

clean :
	@rm -f *~ *.o *.so *.bin
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@echo
	@echo "Cleaned up."

install : seif.bin libseif.so
	@mkdir -p $(INSTALL_DIRS)
	@cp seif.bin  $(INSTALL_BIN_DIR)
	@cp libseif.so $(INSTALL_LIB_DIR)
	@cp *.h  $(INSTALL_INCLUDE_DIR)
	@echo "File installed in " $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR) $(INSTALL_INCLUDE_DIR)

dist :
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@mkdir $(ARCHIVE_DIR_NAME)
	@cp $(FILES_TO_BE_ARCHIVED) $(ARCHIVE_DIR_NAME)
	@tar zcf $(ARCHIVE_NAME) $(ARCHIVE_DIR_NAME)
	@rm -rf $(ARCHIVE_DIR_NAME)
	@echo "File" $(ARCHIVE_NAME) " generated."
