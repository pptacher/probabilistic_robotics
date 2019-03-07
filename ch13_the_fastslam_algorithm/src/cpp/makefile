help:
	@echo
	@echo "make lib      <-- Compile lib."
	@echo "make build    <-- Compile all."
	@echo "make clean    <-- Clean up."
	@echo "make install  <-- Compile and install"
	@echo "make dist     <-- Create .tar.gz project file"
	@echo

PROJECT_NAME = fastslam
VERSION      = 001

ARCHIVE_DIR_NAME = $(PROJECT_NAME)-$(VERSION)
ARCHIVE_NAME     = $(ARCHIVE_DIR_NAME).tar.gz

FILES_TO_BE_ARCHIVED = *.h *.cc Makefile

INCLUDE_PATH = -I. -I/opt/intel/compilers_and_libraries/linux/mkl/include/  -I/usr/include/hdf5/serial
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64:.
COMPILO     =  g++
CXXFLAGS      =  $(RELEASEFLAGS) -m64 -Wall -Wextra -Wshadow -Werror -pedantic   $(INCLUDE_PATH) -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused
CXXFLAGS_SHARED  = $(CXXFLAGS) -fPIC -shared
LIBS        = -L. -L/usr/local/hdf5/lib/ -lfslam -larmadillo -lhdf5_cpp

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

#memos
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64:/usr/local/lib64:/usr/local/lib:.
export LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64
#export MKL_VERBOSE=1
#gnuplot
#plot "poses.dat" using 91:92 pt 7 ps 0.1 lc rgb "black", "ldmarks.dat" using 1:2 pt 7 ps 0.6 lc rgb "green"

COMPILER_OPTIONS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic -Iinclude -Wno-unused-parameter -Wno-unused-variable -Wno-deprecated-declarations
#INCLUDE_PATH = -I.

#CXX_DEFINES = -DSHADERS_DIR=\"/Users/TACHER/sparklepp/src/shaders\" -DUSE_GLEW

#CFLAGS := -std=c11 $(COMPILER_OPTIONS)
#CXXFLAGS := -std=c++14 -Weffc++ $(COMPILER_OPTIONS) $(INCLUDE_PATH)

DEBUGFLAGS := -g -O0 -D _DEBUG
RELEASEFLAGS := -O2 -D NDEBUG
##########################################################

# Fabrication des .o (lib)

particle.o : particle.cc
			$(PROCESS_UTF8) particle.h
			$(PROCESS_UTF8) btree.h
			$(PROCESS_UTF8) particle.cc
			$(COMPILO)  -c -std=c++17 $(CXXFLAGS_SHARED) /tmp/particle.cc

fastslam.o : fastslam.cc
			$(PROCESS_UTF8) particle.h
			$(PROCESS_UTF8) btree.h
			$(PROCESS_UTF8) fastslam.cc
			$(COMPILO)  -c -std=c++17 $(CXXFLAGS_SHARED) /tmp/fastslam.cc

linear_model.o : linear_model.cc
						$(PROCESS_UTF8) linear_model.cc
						$(COMPILO)  -c -std=c++17  $(CXXFLAGS_SHARED) /tmp/linear_model.cc

# Fabrication des .o (hors lib)

main.o : main.cc
			$(COMPILO) -c -std=c++17 $(CXXFLAGS) main.cc

main0.o : main0.cc
			$(PROCESS_UTF8) particle.h
			$(PROCESS_UTF8) main0.cc
			$(COMPILO) -c -std=c++17 $(CXXFLAGS) /tmp/main0.cc

# Fabrication de la lib

libfslam.so : particle.o fastslam.o linear_model.o
	$(COMPILO) -o libfslam.so -shared -larmadillo particle.o fastslam.o linear_model.o

# Fabrication de l'exécutable

fslam.bin : main.o libfslam.so
		$(COMPILO) -o   fslam.bin main.o  $(LIBS) -lpthread -lmkl_core -lmkl_def -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl
#-DARMA_DONT_USE_WRAPPER -lmkl_intel_thread -lmkl_intel_ilp64
test.bin : main0.o libfslam.so
		$(COMPILO) -o   test.bin main0.o  $(LIBS) -lpthread -lmkl_core -lmkl_def -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl
# Les cibles de l'aide

lib : libfslam.so

build : fslam.bin

clean :
	@rm data/poses.dat data/ldmarks.dat
	@rm -f *~ *.o *.so *.bin
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@echo
	@echo "Cleaned up."

install : fslam.bin libfslam.so
	@mkdir -p $(INSTALL_DIRS)
	@cp fslam.bin  $(INSTALL_BIN_DIR)
	@cp data -r  $(INSTALL_BIN_DIR)
	@cp libfslam.so $(INSTALL_LIB_DIR)
	@cp *.h  $(INSTALL_INCLUDE_DIR)
	@echo "File installed in " $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR) $(INSTALL_INCLUDE_DIR)

dist :
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@mkdir $(ARCHIVE_DIR_NAME)
	@cp $(FILES_TO_BE_ARCHIVED) $(ARCHIVE_DIR_NAME)
	@tar zcf $(ARCHIVE_NAME) $(ARCHIVE_DIR_NAME)
	@rm -rf $(ARCHIVE_DIR_NAME)
	@echo "File" $(ARCHIVE_NAME) " generated."
