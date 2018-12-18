help:
	@echo
	@echo "make lib      <-- Compile la lib."
	@echo "make build    <-- Compile tout."
	@echo "make clean    <-- Fais le m�nage."
	@echo "make install  <-- Compile et installe"
	@echo "make dist     <-- Fais un .tar.gz du projet"
	@echo

PROJECT_NAME = dgemm
VERSION      = 001

ARCHIVE_DIR_NAME = $(PROJECT_NAME)-$(VERSION)
ARCHIVE_NAME     = $(ARCHIVE_DIR_NAME).tar.gz

FILES_TO_BE_ARCHIVED = *.hpp *.cpp Makefile

INCLUDE_PATH = -I. -I/opt/intel/compilers_and_libraries/linux/mkl/include/

#doesn't support -std=c++14 when using "g++"
COMPILO     =  g++
CXXFLAGS      =  -Wall   $(INCLUDE_PATH) -Wno-deprecated-declarations
CXXFLAGS_COOK  = $(CXXFLAGS) -fPIC -shared
LIBS        = -L. -ldgemm

INSTALL_PREFIX      = ./LocalInstall

INSTALL_BIN_DIR     = $(INSTALL_PREFIX)/bin
INSTALL_LIB_DIR     = $(INSTALL_PREFIX)/lib
INSTALL_INCLUDE_DIR = $(INSTALL_PREFIX)/include/dgemm

INSTALL_DIRS =					\
	$(INSTALL_BIN_DIR)			\
	$(INSTALL_LIB_DIR)			\
	$(INSTALL_INCLUDE_DIR)

# Pour l'instant, aucun .hpp n'est mis dans les d�pendances.... on
# devrait, mais c'est difficile (car les .hpp peuvent s'inclure les
# uns les autres).
########################################################
#SHELL := /bin/sh
#CC := g++
#CXX := g++

export LIBRARY_PATH=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries/linux/lib/intel64
export MKL_VERBOSE=1

COMPILER_OPTIONS := -m64 -Wall -Wextra -Wshadow -Werror -pedantic -Iinclude     -fdiagnostics-color=auto -lsqlite3 -fno-builtin-qsort -lglfw -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo -Wno-deprecated-declarations
#INCLUDE_PATH = -I.

CXX_DEFINES = -DSHADERS_DIR=\"/Users/TACHER/sparklepp/src/shaders\" -DUSE_GLEW

#CFLAGS := -std=c11 $(COMPILER_OPTIONS)
#CXXFLAGS := -std=c++14 -Weffc++ $(COMPILER_OPTIONS) $(INCLUDE_PATH)

DEBUGFLAGS := -g -O0 -D _DEBUG
RELEASEFLAGS := -O2 -D NDEBUG
##########################################################

# Fabrication des .o (lib)

app.o :	app.cc
			$(COMPILO)  -c -std=c++14  $(CXXFLAGS_COOK) app.cc

# Fabrication des .o (hors lib)

main.o : main.cc
			$(COMPILO) -c -std=c++14 $(CXXFLAGS) main.cc

main1.o : main1.cc
			./to_UCN.sh main1.cc
			$(COMPILO) -c -std=c++14 $(CXXFLAGS) /tmp/main1.cc

main0.o : main0.cc
			$(COMPILO) -c -std=c++14 $(CXXFLAGS) main0.cc

main10.o : main10.cc
				$(COMPILO) -c -std=c++14 $(CXXFLAGS) main10.cc

example1.o : example1.cc
				$(COMPILO) -c -std=c++14 $(CXXFLAGS) example1.cc

# Fabrication de la lib

libsparkle.so : app.o events.o opengl.o scene.o append_consume_buffer.o gpu_particle.o random_buffer.o vector_field.o noise.o
	$(COMPILO) -o libsparkle.so -shared -lglfw3  -lFreetype -lGlew -framework Cocoa -framework OpenGL -framework Glut -framework IOKit -framework CoreVideo  app.o events.o opengl.o scene.o append_consume_buffer.o gpu_particle.o random_buffer.o vector_field.o noise.o

# Fabrication de l'ex�cutable

test2.bin : main.o libsparkle.so
	$(COMPILO) -o   test1.bin main.o  -lsoil -lglfw -lglfw3  -lGlew -framework Cocoa -framework OpenGL -framework Glut -framework IOKit -framework CoreVideo $(LIBS)

test.bin : main.o
		$(COMPILO) -o   test.bin main.o  -larmadillo -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin -lmkl_intel_ilp64  -lmkl_core -lmkl_def -lmkl_sequential

test1.bin : main1.o
		$(COMPILO) -o   test1.bin main1.o   -lmkl_intel_ilp64  -lmkl_core -lmkl_def -lmkl_intel_thread -lpthread -ldl -lm -liomp5

test10.bin : main10.o
		$(COMPILO) -o   test10.bin main10.o   -lmkl_intel_ilp64  -lmkl_core -lmkl_def -lmkl_intel_thread -lpthread -ldl -lm -liomp5 -larmadillo

example1.bin : example1.o
				$(COMPILO) -o   example1.bin example1.o   -lmkl_intel_ilp64  -lmkl_core -lmkl_def -lmkl_intel_thread -lpthread -ldl -lm -liomp5 -larmadillo

test0.bin : main0.o
		$(COMPILO) -o   test0.bin main0.o -larmadillo

#-DARMA_DONT_USE_WRAPPER -lmkl_intel_ilp64  -lmkl_core -lmkl_def -lmkl_intel_thread -lpthread -ldl -lm -liomp5

phony :
	echo $(LIBRARY_PATH)
	@(echo $(MKL_VERBOSE))

# Les cibles de l'aide

lib : libsparkle.so

build : test.bin

clean :
	@rm -f *~ *.o *.so *.bin
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@echo
	@echo "Cleaned up."

install : test.bin libcook.so
	@mkdir -p $(INSTALL_DIRS)
	@cp test.bin  $(INSTALL_BIN_DIR)
	@cp libcook.so $(INSTALL_LIB_DIR)
	@cp cook*.hpp  $(INSTALL_INCLUDE_DIR)
	@echo "File installed in " $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR) $(INSTALL_INCLUDE_DIR)

dist :
	@rm -rf $(ARCHIVE_DIR_NAME) $(ARCHIVE_NAME)
	@mkdir $(ARCHIVE_DIR_NAME)
	@cp $(FILES_TO_BE_ARCHIVED) $(ARCHIVE_DIR_NAME)
	@tar zcf $(ARCHIVE_NAME) $(ARCHIVE_DIR_NAME)
	@rm -rf $(ARCHIVE_DIR_NAME)
	@echo "File" $(ARCHIVE_NAME) " generated."
