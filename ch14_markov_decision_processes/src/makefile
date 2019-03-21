help:
	@echo
	@echo "make lib      <-- Compile lib."
	@echo "make build    <-- Compile all."
	@echo "make clean    <-- Clean up."
	@echo "make install  <-- Compile and install"
	@echo "make dist     <-- Create .tar.gz project file"
	@echo

PROJECT_NAME = markovdp
VERSION      = 001

ARCHIVE_DIR_NAME = $(PROJECT_NAME)-$(VERSION)
ARCHIVE_NAME     = $(ARCHIVE_DIR_NAME).tar.gz

FILES_TO_BE_ARCHIVED = *.h *.cc Makefile

INCLUDE_PATH = -I. -I/opt/intel/compilers_and_libraries/linux/mkl/include/ -I/usr/include/hdf5/serial
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64:.
COMPILO     =  g++
CXXFLAGS      =  $(DEBUGFLAGS) -m64 -Wall -Wextra -Wshadow -Werror -pedantic   $(INCLUDE_PATH) -Wno-deprecated-declarations -Wno-unused-variable -Wno-unused
CXXFLAGS_SHARED  = $(CXXFLAGS) -fPIC -shared
LIBS        = -L.  -L/usr/local/hdf5/lib -lmarkovdp -larmadillo -lhdf5_serial -lhdf5_cpp

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
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/serial:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64:/usr/local/lib64:/usr/local/lib:.
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


# Fabrication des .o (hors lib)

markovdp.o : markovdp.cc
			$(PROCESS_UTF8) markovdp.h
			$(PROCESS_UTF8) markovdp.cc
			$(COMPILO)  -c -std=c++17 $(CXXFLAGS_SHARED) /tmp/markovdp.cc

main.o : main.cc
			$(COMPILO) -c -std=c++17 $(CXXFLAGS) main.cc

# Fabrication de la lib

libmarkovdp.so : markovdp.o
	$(COMPILO) -o libmarkovdp.so  -shared -larmadillo markovdp.o

# Fabrication de l'exécutable

vi.bin : main.o libmarkovdp.so
		$(COMPILO)  -o vi.bin main.o   $(LIBS)  -lpthread -lmkl_core -lmkl_def -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl

# Les cibles de l'aide

lib :

build : vi.bin

clean :
	@rm -f *~ *.o *.so *.bin
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@echo
	@echo "Cleaned up."

install : fslam.bin libfslam.so
	@mkdir -p $(INSTALL_DIRS)
	@cp vi.bin  $(INSTALL_BIN_DIR)
	@cp *.h  $(INSTALL_INCLUDE_DIR)
	@echo "File installed in " $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR) $(INSTALL_INCLUDE_DIR)

dist :
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@mkdir $(ARCHIVE_DIR_NAME)
	@cp $(FILES_TO_BE_ARCHIVED) $(ARCHIVE_DIR_NAME)
	@tar zcf $(ARCHIVE_NAME) $(ARCHIVE_DIR_NAME)
	@rm -rf $(ARCHIVE_DIR_NAME)
	@echo "File" $(ARCHIVE_NAME) " generated."
