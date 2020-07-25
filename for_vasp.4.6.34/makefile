.SUFFIXES: .inc .f .f90 .F
#-----------------------------------------------------------------------
# Makefile for Intel Fortran compiler for Operton systems
#
# The makefile was tested only under Linux on Intel platforms
# (Suse 5.3- Suse 9.2)
# the followin compiler versions have been tested
# 5.0, 6.0, 7.0 and 7.1 (some 8.0 versions seem to fail compiling the code,
#  8.1 is slower than 8.0)
# presently we recommend version 7.1 or 7.0, since these
# releases have been used to compile the present code versions
#
# it might be required to change some of library pathes, since
# LINUX installation vary a lot
# Hence check ***ALL**** options in this makefile very carefully
#-----------------------------------------------------------------------
#
# BLAS must be installed on the machine
# there are several options:
#   Kazushige Goto's BLAS is required
#   http://www.cs.utexas.edu/users/kgoto/signup_first.html
# (see below libgoto comments)
#
# FFT:
#   the fftw.3.0.1 must be available and installed, since
#   the ifc compiler creates crap code if the build in fft routines are used
# (see below fftw comments)
#-----------------------------------------------------------------------

# all CPP processed fortran files have the extension .f90
SUFFIX=.f90

#-----------------------------------------------------------------------
# fortran compiler and linker
#-----------------------------------------------------------------------
FC=ifort
# fortran linker
FCL=$(FC)


#-----------------------------------------------------------------------
# whereis CPP ?? (I need CPP, can't use gcc with proper options)
# that's the location of gcc for SUSE 5.3
#
#  CPP_   =  /usr/lib/gcc-lib/i486-linux/2.7.2/cpp -P -C 
#
# that's probably the right line for some Red Hat distribution:
#
#  CPP_   =  /usr/lib/gcc-lib/i386-redhat-linux/2.7.2.3/cpp -P -C
#
#  SUSE X.X, maybe some Red Hat distributions:

CPP_ =  ./preprocess <$*.F | /usr/bin/cpp -P -C -traditional >$*$(SUFFIX)

#-----------------------------------------------------------------------
# possible options for CPP:
# NGXhalf             charge density   reduced in X direction
# wNGXhalf            gamma point only reduced in X direction
# avoidalloc          avoid ALLOCATE if possible
# IFC                 work around some IFC bugs
# CACHE_SIZE          1000 for PII,PIII, 5000 for Athlon, 8000-12000 P4
# RPROMU_DGEMV        use DGEMV instead of DGEMM in RPRO (depends on used BLAS)
# RACCMU_DGEMV        use DGEMV instead of DGEMM in RACC (depends on used BLAS)
#-----------------------------------------------------------------------

CPP     = $(CPP_)  -DHOST=\"LinuxIFC\" \
          -Dkind8 -DNGXhalf  -DCACHE_SIZE=20000 -DPGF90 -Davoidalloc \
          -Duse_cray_ptr   -DRPROMU_DGEMV  -DRACCMU_DGEMV

#-----------------------------------------------------------------------
# general fortran flags  (there must a trailing blank on this line)
#-----------------------------------------------------------------------

FFLAGS =  -FR 
#-assume byterecl 

#-----------------------------------------------------------------------
# optimization
# we have tested whether higher optimisation improves performance
# -axK  SSE1 optimization,  but also generate code executable on all mach.
#       xK improves performance somewhat on XP, and a is required in order
#       to run the code on older Athlons as well
# -xW   SSE2 optimization
# -axW  SSE2 optimization,  but also generate code executable on all mach.
# -tpp6 P3 optimization
# -tpp7 P4 optimization
#-----------------------------------------------------------------------

OFLAG=-O3  -axSSE4.1
#-tpp7

OFLAG_HIGH = $(OFLAG)
OBJ_HIGH = 

OBJ_NOOPT = 
DEBUG  = -FR -O0
INLINE = $(OFLAG)


#-----------------------------------------------------------------------
# the following lines specify the position of BLAS  and LAPACK
# on Operton, you really need the  libgoto library
#-----------------------------------------------------------------------

#BLAS= $(HOME)/OPT/src/vasp.4.6/libgoto_opteron-32-r0.99.so -lpthread
#BLAS= $(HOME)/OPT/src/vasp.4.6/libgoto_opt32-r0.96.so -lpthread
#BLAS=  /opt/libs/libgoto/libgoto_p4_512-r0.6.so
#BLAS=/home/haiping/vasp/vasp.4.lib/libgoto_prescott32p-r1.00.so -lpthread
#build from src_pack
#BLAS=/home/haiping/vasp/vasp.4.lib/libgoto_core2p-r1.16.so -lpthread
#mkl 
#BLAS= -L/opt/intel/cmkl/8.0/lib/32 -lmkl_p4 -lguide -lpthread -lsvml
#BLAS=-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_em64t  -lguide -lpthread
BLAS=-L/opt/intel/mkl/lib/intel64/ -lirc -lpthread -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -openmp

# LAPACK, simplest use vasp.5.lib/lapack_double
LAPACK=../vasp.4.lib/lapack_double.o 

# use the mkl Intel lapack
#LAPACK= -lmkl_lapack

#-----------------------------------------------------------------------

LIB  = -L../vasp.4.lib -ldmy \
     ../vasp.4.lib/linpack_double.o $(LAPACK) \
     $(BLAS)

# options for linking, nothing is required (usually)
LINK    = /opt/intel/lib/intel64/libsvml.so /opt/intel/lib/intel64/libintlc.so.5 
# compiler version 7.0 generates some vector statments which are located
# in the svml library, add the LIBPATH and the library (just in case)
#LINK    =  -L/opt/intel/compiler70/ia32/lib/ -lsvml 

#-----------------------------------------------------------------------
# fft libraries:
# On Operton you really have to used fftw.3.0.X (http://www.fftw.org)
# the ifc compiler creates suboptimal performance on the Opteron for
# the build in fft routines
#
# fftw.3.0.1 was compiled using the following command lines:
# > export CC="gcc -m32"
# > export F77="f77 -m32"
# > ./configure  --enable-sse2 --prefix=/home/kresse/ifc_opt/fftw-3.0.1/
# > make 
# > make install
# PLEASE do not send querries related to fftw to the vasp site
#-----------------------------------------------------------------------

FFT3D   = fft3dfurth.o fft3dlib.o
#FFT3D   = fftw3d.o fft3dlib.o   /opt/libs/fftw-3.0.1/lib/libfftw3.a


#=======================================================================
# MPI section, uncomment the following lines
# 
# one comment for users of mpich or lam:
# You must *not* compile mpi with g77/f77, because f77/g77             
# appends *two* underscores to symbols that contain already an        
# underscore (i.e. MPI_SEND becomes mpi_send__).  The pgf90/ifc
# compilers however append only one underscore.
# Precompiled mpi version will also not work !!!
#
# We found that mpich.1.2.1 and lam-6.5.X to lam-7.0.4 are stable
# mpich.1.2.1 was configured with 
#  ./configure -prefix=/usr/local/mpich_nodvdbg -fc="pgf77 -Mx,119,0x200000"  \
# -f90="pgf90 " \
# --without-romio --without-mpe -opt=-O \
# 
# lam was configured with the line
#  ./configure  -prefix /opt/libs/lam-7.0.4 --with-cflags=-O -with-fc=ifc \
# --with-f77flags=-O --without-romio
# 
# please note that you might be able to use a lam or mpich version 
# compiled with f77/g77, but then you need to add the following
# options: -Msecond_underscore (compilation) and -g77libs (linking)
#
# !!! Please do not send me any queries on how to install MPI, I will
# certainly not answer them !!!!
#=======================================================================
#-----------------------------------------------------------------------
# fortran linker for mpi: if you use LAM and compiled it with the options
# suggested above,  you can use the following line
#-----------------------------------------------------------------------

FC=mpif90
FCL=$(FC)

#-----------------------------------------------------------------------
# additional options for CPP in parallel version (see also above):
# NGZhalf               charge density   reduced in Z direction
# wNGZhalf              gamma point only reduced in Z direction
# scaLAPACK             use scaLAPACK (usually slower on 100 Mbit Net)
#-----------------------------------------------------------------------

CPP    = $(CPP_) -DMPI  -DHOST=\"LinuxIFC\" -DIFC \
     -Dkind8 -DNGZhalf  -DCACHE_SIZE=10000 -DPGF90 -Davoidalloc \
     -DMPI_BLOCK=2000  \
     -Duse_cray_ptr
#    -DRPROMU_DGEMV  -DRACCMU_DGEMV

#-----------------------------------------------------------------------
# location of SCALAPACK
# if you do not use SCALAPACK simply uncomment the line SCA
#-----------------------------------------------------------------------

#BLACS=$(HOME)/archives/SCALAPACK/BLACS/
#SCA_=$(HOME)/archives/SCALAPACK/SCALAPACK

#SCA= $(SCA_)/libscalapack.a  \
# $(BLACS)/LIB/blacsF77init_MPI-LINUX-0.a $(BLACS)/LIB/blacs_MPI-LINUX-0.a $(BLACS)/LIB/blacsF77init_MPI-LINUX-0.a

#SCA=

#-----------------------------------------------------------------------
# libraries for mpi
#-----------------------------------------------------------------------

#LIB     = -L../vasp.4.lib -ldmy  \
#      ../vasp.4.lib/linpack_double.o $(LAPACK) \
#      $(SCA) $(BLAS)

# FFT: fftmpi.o with fft3dlib of Juergen Furthmueller
FFT3D   = fftmpi.o fftmpi_map.o fft3dlib.o 

# fftw.3.0.1 is much faster on Opteron
#FFT3D   = fftmpiw.o fftmpi_map.o fft3dlib.o   /opt/libs/fftw-3.0.1/lib/libfftw3.a

#-----------------------------------------------------------------------
# general rules and compile lines
#-----------------------------------------------------------------------
BASIC=   symmetry.o symlib.o   lattlib.o  random.o   

SOURCE=  base.o     mpi.o      smart_allocate.o      xml.o  \
         constant.o jacobi.o   main_mpi.o  scala.o   \
         asa.o      lattice.o  poscar.o   ini.o      setex.o     radial.o  \
         pseudo.o   mgrid.o    mkpoints.o wave.o      wave_mpi.o  $(BASIC) \
         nonl.o     nonlr.o    dfast.o    choleski2.o    \
         mix.o      charge.o   xcgrad.o   xcspin.o    potex1.o   potex2.o  \
         metagga.o  constrmag.o pot.o      cl_shift.o force.o    dos.o      elf.o      \
         tet.o      hamil.o    steep.o    \
         chain.o    dyna.o     relativistic.o LDApU.o sphpro.o  paw.o   us.o \
         ebs.o      wavpre.o   wavpre_noio.o broyden.o \
         dynbr.o    rmm-diis.o reader.o   writer.o   tutor.o xml_writer.o \
         brent.o    stufak.o   fileio.o   opergrid.o stepver.o  \
         dipol.o    xclib.o    chgloc.o   subrot.o   optreal.o   davidson.o \
         edtest.o   electron.o shm.o      pardens.o  paircorrection.o \
         optics.o   constr_cell_relax.o   stm.o    finite_diff.o \
         elpol.o    setlocalpp.o  aedens.o getnabij.o 
 
INC=

vasp: $(SOURCE) $(FFT3D) $(INC) main.o 
	rm -f vasp
	$(FCL) -o vasp $(LINK) main.o  $(SOURCE)   $(FFT3D) $(LIB) 
makeparam: $(SOURCE) $(FFT3D) makeparam.o main.F $(INC)
	$(FCL) -o makeparam  $(LINK) makeparam.o $(SOURCE) $(FFT3D) $(LIB)
zgemmtest: zgemmtest.o base.o random.o $(INC)
	$(FCL) -o zgemmtest $(LINK) zgemmtest.o random.o base.o $(LIB)
dgemmtest: dgemmtest.o base.o random.o $(INC)
	$(FCL) -o dgemmtest $(LINK) dgemmtest.o random.o base.o $(LIB) 
ffttest: base.o smart_allocate.o mpi.o mgrid.o random.o ffttest.o $(FFT3D) $(INC)
	$(FCL) -o ffttest $(LINK) ffttest.o mpi.o mgrid.o random.o smart_allocate.o base.o $(FFT3D) $(LIB)
kpoints: $(SOURCE) $(FFT3D) makekpoints.o main.F $(INC)
	$(FCL) -o kpoints $(LINK) makekpoints.o $(SOURCE) $(FFT3D) $(LIB)

clean:	
	-rm -f *.g *.f *.o *.L *.mod ; touch *.F

main.o: main$(SUFFIX)
	$(FC) $(FFLAGS)$(DEBUG)  $(INCS) -c main$(SUFFIX)
xcgrad.o: xcgrad$(SUFFIX)
	$(FC) $(FFLAGS) $(INLINE)  $(INCS) -c xcgrad$(SUFFIX)
xcspin.o: xcspin$(SUFFIX)
	$(FC) $(FFLAGS) $(INLINE)  $(INCS) -c xcspin$(SUFFIX)

makeparam.o: makeparam$(SUFFIX)
	$(FC) $(FFLAGS)$(DEBUG)  $(INCS) -c makeparam$(SUFFIX)

makeparam$(SUFFIX): makeparam.F main.F 
#
# MIND: I do not have a full dependency list for the include
# and MODULES: here are only the minimal basic dependencies
# if one strucuture is changed then touch_dep must be called
# with the corresponding name of the structure
#
base.o: base.inc base.F
mgrid.o: mgrid.inc mgrid.F
constant.o: constant.inc constant.F
lattice.o: lattice.inc lattice.F
setex.o: setexm.inc setex.F
pseudo.o: pseudo.inc pseudo.F
poscar.o: poscar.inc poscar.F
mkpoints.o: mkpoints.inc mkpoints.F
wave.o: wave.inc wave.F
nonl.o: nonl.inc nonl.F
nonlr.o: nonlr.inc nonlr.F

$(OBJ_HIGH):
	$(CPP)
	$(FC) $(FFLAGS) $(OFLAG_HIGH) $(INCS) -c $*$(SUFFIX)
$(OBJ_NOOPT):
	$(CPP)
	$(FC) $(FFLAGS) $(INCS) -c $*$(SUFFIX)

fft3dlib_f77.o: fft3dlib_f77.F
	$(CPP)
	$(F77) $(FFLAGS_F77) -c $*$(SUFFIX)

.F.o:
	$(CPP)
	$(FC) $(FFLAGS) $(OFLAG) $(INCS) -c $*$(SUFFIX)
.F$(SUFFIX):
	$(CPP)
$(SUFFIX).o:
	$(FC) $(FFLAGS) $(OFLAG) $(INCS) -c $*$(SUFFIX)

# special rules
#-----------------------------------------------------------------------
# these special rules are cummulative (that is once failed
#   in one compiler version, stays in the list forever)
# -tpp5|6|7 P, PII-PIII, PIV
# -xW use SIMD (does not pay of on PII, since fft3d uses double prec)
# all other options do no affect the code performance since -O1 is used
#-----------------------------------------------------------------------

fft3dlib.o : fft3dlib.F
	$(CPP)
	$(FC) -FR -lowercase -O1  -xW  -unroll0  -c $*$(SUFFIX)
#	$(FC) -FR -lowercase -O1  -c $*$(SUFFIX)
fft3dfurth.o : fft3dfurth.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)
fftmpi.o : fftmpi.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -xW -c $*$(SUFFIX)
fftmpi_map.o : fftmpi_map.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -xW -c $*$(SUFFIX)

fftw3d.o : fftw3d.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -xW -c $*$(SUFFIX)

radial.o : radial.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

symlib.o : symlib.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

symmetry.o : symmetry.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

dynbr.o : dynbr.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

broyden.o : broyden.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

us.o : us.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)

wave.o : wave.F
	$(CPP)
	$(FC) -FR -lowercase -O0 -c $*$(SUFFIX)

LDApU.o : LDApU.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)

