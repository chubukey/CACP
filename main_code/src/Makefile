
TARGET      = ./lmp_CACP

CPP_SRC     = 
C_SRC       = ./superlu_c2f_dwrap.c ./libfwrapper.c ./dcreate_dist_matrix.c
F_SRC       = ./fem_uel.f ./umat_nonlinear.f ./error.f ./ESCM_couple.f ./loadvc_C.f ./uel_disp_brick.f ./disinc_C.f ./DG_calc.f ./uel_disp_t.f ./init_slu.f
F90_SRC       =./inverse_slu_dist.f90
MOD_SRC		= ./superlupara.f90 ./superlu_mod.f90
SRC         = ${CPP_SRC} ${C_SRC} ${F_SRC} ${F90_SRC}

HEADER      = ./pardis.h ./ESCM_couple.h

OBJ         = ${CPP_SRC:.C=.o} ${C_SRC:.c=.o} ${MOD_SRC:.f90=.o} ${F_SRC:.f=.o} ${F90_SRC:.f90=.o}
#OBJ         = ${CPP_SRC:.C=.o} ${MOD_SRC:.f90=.o} ${F_SRC:.f=.o} ${F90_SRC:.f90=.o}

#--gcc
#DEBUG_FLAGS = -g -fbounds-check -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wuninitialized 
#RELEASE_FLAGS = -O3
#--intel fortran
#DEBUG_FLAGS = -g -debug inline-debug-info -debug extended -check all -ftrapuv -fpe-all=0 -traceback
#DEBUG_FLAGS = -g -debug inline-debug-info -check bounds -check uninit -ftrapuv -fpe-all=0 -traceback
#RELEASE_FLAGS = -O2 
RELEASE_FLAGS = -O2

#--gcc
#CC         = mpicc
#F77        = mpif90
#F90        = mpif90
#--intel fortran
CC         = mpiicc
F77        = mpiifort
F90        = mpiifort

SUPERLU_LIBS=-L/cm/shared/apps/SuperLU/4.1/lib -lsuperlu_dist_4.1
SUPERLU_INC=-I/cm/shared/apps/SuperLU/4.1/include
BLAS_LIBS=-L${MKL}
METIS_LIB=/cm/shared/apps/metis/5.1.0_32/lib/libmetis.a
PARMETIS_LIB=/cm/shared/apps/parmetis/4.0.3_32/lib/libparmetis.a
export METISPARMETIS_LIBS=${PARMETIS_LIB} ${METIS_LIB}




CGAL_INCLUDE=-I/home-2/jzhang74@jhu.edu/opt/local_boost/CGAL-4.9/include
EIGEN_INCLUDE=-I/home-2/jzhang74@jhu.edu/usr/local/include/eigen3
CRYSTALANALYSIS_INCLUDE=-I/home-2/jzhang74@jhu.edu/opt/local_boost/CrystalAnalysis150115/src
BOOST_INCLUDE= -I/home-2/jzhang74@jhu.edu/opt/local_boost/boost_1_55_0
GMP_INCLUDE=-I/home-2/jzhang74@jhu.edu/usr/local/include
DXA_INC=${CGAL_INCLUDE} ${EIGEN_INCLUDE} ${CRYSTALANALYSIS_INCLUDE} ${GMP_INCLUDE} ${BOOST_INCLUDE}

CGAL_LIB=-L/home-2/jzhang74@jhu.edu/opt/local_boost/CGAL-4.9/lib -lCGAL
CRYSTALANALYSIS_LIB=-L/home-2/jzhang74@jhu.edu/opt/local_boost/CrystalAnalysis150115/build/calib -lcrystalanalysis
BOOST_LIB=-L/home-2/jzhang74@jhu.edu/opt/local_boost/boost_1_55_0/stage/lib -lboost_iostreams -lboost_serialization -lboost_thread
GMP_LIB=-L/home-2/jzhang74@jhu.edu/usr/local/lib -lgmp
DXA_LIBS=${CGAL_LIB} ${CRYSTALANALYSIS_LIB} ${GMP_LIB} ${BOOST_LIB}


LAMMPS_INC=-I/home-2/jzhang74@jhu.edu/work/jiaxi/lammps/lammps_DXAtest/src ${DXA_INC}
LAMMPS_LIBS=-L/home-2/jzhang74@jhu.edu/work/jiaxi/lammps/lammps_DXAtest/src -llammps_cacp ${DXA_LIBS}





#--------------------------------------------------------------------------------#
INCLUDES = ${SUPERLU_INC} ${LAMMPS_INC}

C_64 = -m64 -D_LONGINT
F_64 = -i8
CF_lammps    = -g -O2 -DGZIP -DLINUX -fPIC 
other_libs =   -lstdc++
#-----------------------------------------------------------------------------#

CFLAGS   = -mcmodel=large ${C_64} ${CF_lammps}
FFLAGS	 = -mcmodel=large ${F_64}
F90FLAGS = -mcmodel=large -no-wrap-margin ${F_64}
#other_libs = -L/cm/shared/apps/Intel/openmpi/1.8.4/lib  -lmpi_cxx -lstdc++
LIBS  = ${SUPERLU_LIBS} ${BLAS_LIBS} ${METISPARMETIS_LIBS} ${LAMMPS_LIBS} ${other_libs}

# DEFAULT:
# adopt high-performance/release flags
release: FFLAGS += ${RELEASE_FLAGS}
release: F90FLAGS += ${RELEASE_FLAGS}
release: ${TARGET}

# adopt low-performance/debug flags
debug: FFLAGS += ${DEBUG_FLAGS}
debug: F90FLAGS += ${DEBUG_FLAGS}
debug: ${TARGET}

# linker
${TARGET}: ${OBJ}
	${F90} ${F90FLAGS} ${OBJ} ${LIBS} -o ${TARGET}

# recompile object files if header file changes
${OBJ}: ${HEADER}

runUmat: ${RUNUMAT_SRC:.f90=.o}
	gfortran ${RUNUMAT_SRC} ${BLAS_LIBS} ${F90FLAGS} ${DEBUG_FLAGS} -fno-align-commons -o runUmat
clean:
	rm -f ${OBJ}
	rm *.mod
tfile:
	cat ${SRC} ${HEADER} > d.d
	
%.o: %.mod
	
.SUFFIXES: .C .c .f .f90 .o

# rules for compiling source files into object files
#
.C.o:  
	${CC} ${CFLAGS} -c $<
.c.o:  
	${CC} ${INCLUDES} ${CFLAGS} -c $<
.f.o:
	${F77} ${INCLUDES} ${FFLAGS} -c $<
.f90.o:
	${F90} ${INCLUDES} ${F90FLAGS} -c $<
