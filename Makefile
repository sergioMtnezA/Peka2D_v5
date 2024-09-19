# Set no debug  by default
# Specify in command line: make DEBUG=yes to set debug
DEBUG=no
MUDFLAP=no
DLEVEL=2
WARN=no
OMP=no
FRICTION=0
COMPILER_OPT=yes

##########################################################################################################################
#	valgrind --tool=memcheck --leak-check=full --track-origins=yes --show-reachable=yes
##########################################################################################################################

# Define C compiler
CC = g++

# Define C flags
C_PROFILE_FLAGS = -pg
C_REQ_FLAGS =
C_DEBUG_FLAGS = -g
#C_DEBUG_FLAGS = -static-libgcc -static-libstdc++ -ggdb
C_OPT_FLAGS = -O3 -Wno-unused-result -Wno-format
C_WARN_FLAGS = -Wall
C_OMP_FLAGS = -fopenmp
OPENMPACT = 0


CFLAGS=$(C_REQ_FLAGS)

ifeq ($(COMPILER_OPT),yes)
	CFLAGS += $(C_OPT_FLAGS)
endif

ifeq ($(WARN),yes)
	CFLAGS += $(C_WARN_FLAGS)
endif

ifeq ($(PROFILE),yes)
	CFLAGS += $(C_PROFILE_FLAGS)
endif

ifeq ($(OMP),yes)
	CFLAGS += $(C_OMP_FLAGS)
endif

ifeq ($(DEBUG),yes)
	CDEBUG=$(C_DEBUG_FLAGS)$(D_LEVEL)

	CFLAGS+=$(CDEBUG)
endif

ifneq ($(FRICTION),0)
	CFLAGS+=-DHYDRO_FRICTION=$(FRICTION)
endif

#LDFLAGS=$(C_FLAGS) -lrt


#OBJ = src/mesh.o src/roe.o src/shallow_water.o src/read_io.o src/boundaries.o src/hydronia.o src/sediments.o src/roeSedim.o src/roeMud.o src/loader.o src/libSwapEndian.o src/simplify1D.o peka2d.o
OBJ = src/utilities.o src/water.o src/calculate.o src/writeInOut.o src/initialize.o src/mesh.o src/loadData.o src/mainKernel.o peka2d.o
BIN = peka2d

## $@ Es el token de la izquierda de los :
## $< Es el token de la derecha de los :

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@  $(OBJ) -lm

clean:
	$(RM) $(OBJ) $(BIN)	


# ------------------------------------------------------
peka2d.o: peka2d.cpp
	$(CC) $(CFLAGS) $< -c -o $@

src/mainKernel.o: src/mainKernel.cpp
	$(CC) $(CFLAGS) $< -c -o $@

src/loadData.o: src/loadData.cpp
	$(CC) $(CFLAGS) $< -c -o $@

src/mesh.o: src/mesh.cpp
	$(CC) $(CFLAGS) $< -c -o $@	

src/initialize.o: src/initialize.cpp
	$(CC) $(CFLAGS) $< -c -o $@		

# CPU-KERNEL compilation -----------------------------
src/writeInOut.o: src/writeInOut.cpp
	$(CC) $(CFLAGS) $< -c -o $@

src/calculate.o: src/calculate.cpp
	$(CC) $(CFLAGS) $< -c -o $@

src/water.o: src/water.cpp
	$(CC) $(CFLAGS) $< -c -o $@	

src/utilities.o: src/utilities.cpp
	$(CC) $(CFLAGS) $< -c -o $@	
