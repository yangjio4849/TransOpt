.SUFFIXES: .inc .f .f90 .F
#-----------------------------------------------------------------------
# Makefile for Intel Fortran compiler for Pentium/Athlon/Opteron 
# based systems
# we recommend this makefile for both Intel as well as AMD systems
# for AMD based systems appropriate BLAS (libgoto) and fftw libraries are
# however mandatory (whereas they are optional for Intel platforms)
# For Athlon we recommend 
#  ) to link against libgoto (and mkl as a backup for missing routines)
#  ) odd enough link in libfftw3xf_intel.a (fftw interface for mkl)
# feedback is greatly appreciated
#
# The makefile was tested only under Linux on Intel and AMD platforms
# the following compiler versions have been tested:
#  - ifc.7.1  works stable somewhat slow but reliably
#  - ifc.8.1  fails to compile the code properly
#  - ifc.9.1  recommended (both for 32 and 64 bit)
#  - ifc.10.1 partially recommended (both for 32 and 64 bit)
#             tested build 20080312 Package ID: l_fc_p_10.1.015
#             the gamma only mpi version can not be compiles
#             using ifc.10.1
#  - ifc.11.1 partially recommended (some problems with Gamma only and intel fftw)
#             Build 20090630 Package ID: l_cprof_p_11.1.046
#  - ifort.12.1 strongly recommended (we use this to compile vasp)
#             Version 12.1.5.339 Build 20120612
#
# it might be required to change some of library path ways, since
# LINUX installations vary a lot
#
# Hence check ***ALL*** options in this makefile very carefully
#-----------------------------------------------------------------------
#
# BLAS must be installed on the machine
# there are several options:
# 1) very slow but works:
#   retrieve the lapackage from ftp.netlib.org
#   and compile the blas routines (BLAS/SRC directory)
#   please use g77 or f77 for the compilation. When I tried to
#   use pgf77 or pgf90 for BLAS, VASP hang up when calling
#   ZHEEV  (however this was with lapack 1.1 now I use lapack 2.0)
# 2) more desirable: get an optimized BLAS 
#
# the two most reliable packages around are presently:
# 2a) Intels own optimised BLAS (PIII, P4, PD, PC2, Itanium)
#     http://developer.intel.com/software/products/mkl/
#   this is really excellent, if you use Intel CPU's
#
# 2b) probably fastest SSE2 (4 GFlops on P4, 2.53 GHz, 16 GFlops PD, 
#     around 30 GFlops on Quad core)
#   Kazushige Goto's BLAS
#   http://www.cs.utexas.edu/users/kgoto/signup_first.html
#   http://www.tacc.utexas.edu/resources/software/
#
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

# this release should be fpp clean
# we now recommend fpp as preprocessor
# if this fails go back to cpp
CPP_=fpp -f_com=no -free -w0 $*.F $*$(SUFFIX) 

#-----------------------------------------------------------------------
# possible options for CPP:
# NGXhalf             charge density   reduced in X direction
# wNGXhalf            gamma point only reduced in X direction
# avoidalloc          avoid ALLOCATE if possible
# PGF90               work around some for some PGF90 / IFC bugs
# CACHE_SIZE          1000 for PII,PIII, 5000 for Athlon, 8000-12000 P4, PD
# RPROMU_DGEMV        use DGEMV instead of DGEMM in RPRO (depends on used BLAS)
# RACCMU_DGEMV        use DGEMV instead of DGEMM in RACC (depends on used BLAS)
# tbdyn                 MD package of Tomas  Bucko
#-----------------------------------------------------------------------

CPP     = $(CPP_)  -DHOST=\"LinuxIFC\" \
          -DCACHE_SIZE=12000 -DPGF90 -Davoidalloc -DNGXhalf \
#          -DRPROMU_DGEMV  -DRACCMU_DGEMV

#-----------------------------------------------------------------------
# general fortran flags  (there must a trailing blank on this line)
# byterecl is strictly required for ifc, since otherwise
# the WAVECAR file becomes huge
#-----------------------------------------------------------------------

FFLAGS =  -FR -names lowercase -assume byterecl 

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

# ifc.9.1, ifc.10.1 recommended
OFLAG=-O2 -ip  

OFLAG_HIGH = $(OFLAG)
OBJ_HIGH = 
OBJ_NOOPT = 
DEBUG  = -FR -O0
INLINE = $(OFLAG)

#-----------------------------------------------------------------------
# the following lines specify the position of BLAS  and LAPACK
# we recommend to use mkl, that is simple and most likely 
# fastest in Intel based machines
#-----------------------------------------------------------------------

# mkl path for ifc 11 compiler
#MKL_PATH=$(MKLROOT)/lib/em64t

# mkl path for ifc 12 compiler
MKLROOT=/opt/intel/mkl
MKL_PATH=$(MKLROOT)/lib/intel64

MKL_FFTW_PATH=$(MKLROOT)/interfaces/fftw3xf

# BLAS
# setting -DRPROMU_DGEMV  -DRACCMU_DGEMV in the CPP lines usually speeds up program execution
#BLAS= -L/opt/intel/mkl/lib/intel64/ -lirc -lpthread -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -openmp
#BLAS= -Wl,--start-group -L$(MKL_PATH)/libmkl_intel_lp64 -L$(MKL_PATH)/libmkl_intel_thread -L$(MKL_PATH)/libmkl_core -L$(MKL_PATH)/libmkl_sequential -Wl,--end-group -lirc -lpthread -openmp #-lguide
BLAS= -L/opt/intel/mkl/lib/intel64/ -lirc -lpthread -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -openmp
# faster linking and available from at least version 11
#BLAS= -lguide  -mkl 


# LAPACK, use vasp.5.lib/lapack_double

LAPACK= ../vasp.5.lib/lapack_double.o -L/opt/intel/lib/intel64/libsvml

# LAPACK from mkl, usually faster and contains scaLAPACK as well

#LAPACK= -L$(MKL_PATH)/libmkl_intel_lp64 -L/opt/intel/lib/intel64/libsvml

# here a tricky version, link in libgoto and use mkl as a backup
# also needs a special line for LAPACK
# this is the best thing you can do on AMD based systems !!!!!!

#BLAS =  -Wl,--start-group /opt/libs/libgoto/libgoto.so $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group  -liomp5
#LAPACK= /opt/libs/libgoto/libgoto.so $(MKL_PATH)/libmkl_intel_lp64.a

#-----------------------------------------------------------------------

LIB  = -L../vasp.5.lib -ldmy \
     ../vasp.5.lib/linpack_double.o $(LAPACK) \
     $(BLAS)

# options for linking, nothing is required (usually)
LINK =

#-----------------------------------------------------------------------
# fft libraries:
# VASP.5.2 can use fftw.3.1.X (http://www.fftw.org)
# since this version is faster on P4 machines, we recommend to use it
#-----------------------------------------------------------------------

FFT3D   = fft3dfurth.o fft3dlib.o

# alternatively: fftw.3.1.X is slighly faster and should be used if available
#FFT3D   = fftw3d.o fft3dlib.o   /opt/libs/fftw-3.1.2/lib/libfftw3.a

# you may also try to use the fftw wrapper to mkl (but the path might vary a lot)
# it seems this is best for AMD based systems
#FFT3D   = fftw3d.o fft3dlib.o $(MKL_FFTW_PATH)/libfftw3xf_intel.a
#INCS = -I$(MKLROOT)/include/fftw 

#=======================================================================
# MPI section, uncomment the following lines until 
#    general  rules and compile lines
# presently we recommend OPENMPI, since it seems to offer better
# performance than lam or mpich
# 
# !!! Please do not send me any queries on how to install MPI, I will
# certainly not answer them !!!!
#=======================================================================
#-----------------------------------------------------------------------
# fortran linker for mpi
#-----------------------------------------------------------------------

FC=mpif90
#FCL=$(FC)

#-----------------------------------------------------------------------
# additional options for CPP in parallel version (see also above):
# NGZhalf             charge density   reduced in Z direction
# wNGZhalf            gamma point only reduced in Z direction
# scaLAPACK           use scaLAPACK (recommended if mkl is available)
# avoidalloc          avoid ALLOCATE if possible
# PGF90               work around some for some PGF90 / IFC bugs
# CACHE_SIZE          1000 for PII,PIII, 5000 for Athlon, 8000-12000 P4, PD
# RPROMU_DGEMV        use DGEMV instead of DGEMM in RPRO (depends on used BLAS)
# RACCMU_DGEMV        use DGEMV instead of DGEMM in RACC (depends on used BLAS)
# tbdyn                 MD package of Tomas  Bucko
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

CPP    = $(CPP_) -DMPI  -DHOST=\"LinuxIFC\" -DIFC \
     -DCACHE_SIZE=12000 -DPGF90 -Davoidalloc -DNGZhalf \
     -DMPI_BLOCK=8000   -DRPROMU_DGEMV  -DRACCMU_DGEMV

#-----------------------------------------------------------------------
# location of SCALAPACK
# if you do not use SCALAPACK simply leave this section commented out
#-----------------------------------------------------------------------

# usually simplest link in mkl scaLAPACK
#BLACS= -lmkl_blacs_openmpi_lp64
#SCA= $(MKL_PATH)/libmkl_scalapack_lp64.a $(BLACS)

#-----------------------------------------------------------------------
# libraries
#-----------------------------------------------------------------------

#LIB     = -L../vasp.5.lib -ldmy  \
#      ../vasp.5.lib/linpack_double.o \
#      $(SCA) $(LAPACK) $(BLAS)

#-----------------------------------------------------------------------
# parallel FFT
#-----------------------------------------------------------------------

# FFT: fftmpi.o with fft3dlib of Juergen Furthmueller
#FFT3D   = fftmpi.o fftmpi_map.o fft3dfurth.o fft3dlib.o 

# alternatively: fftw.3.1.X is slighly faster and should be used if available
#FFT3D   = fftmpiw.o fftmpi_map.o fftw3d.o fft3dlib.o  /opt/libs/fftw-3.1.2/lib/libfftw3.a

# you may also try to use the fftw wrapper to mkl (but the path might vary a lot)
# it seems this is best for AMD based systems
FFT3D   = fftmpiw.o fftmpi_map.o  fftw3d.o  fft3dlib.o   $(MKL_FFTW_PATH)/libfftw3xf_intel.a
#INCS = -I$(MKLROOT)/include/fftw 

#-----------------------------------------------------------------------
# general rules and compile lines
#-----------------------------------------------------------------------
BASIC=   symmetry.o symlib.o   lattlib.o  random.o   


SOURCE=  base.o     mpi.o      smart_allocate.o      xml.o  \
         constant.o jacobi.o   main_mpi.o  scala.o   \
         asa.o      lattice.o  poscar.o   ini.o  mgrid.o  xclib.o  vdw_nl.o  xclib_grad.o \
         radial.o   pseudo.o   gridq.o     ebs.o  \
         mkpoints.o wave.o     wave_mpi.o  wave_high.o  spinsym.o \
         $(BASIC)   nonl.o     nonlr.o    nonl_high.o dfast.o    choleski2.o \
         mix.o      hamil.o    xcgrad.o   xcspin.o    potex1.o   potex2.o  \
         constrmag.o cl_shift.o relativistic.o LDApU.o \
         paw_base.o metagga.o  egrad.o    pawsym.o   pawfock.o  pawlhf.o   rhfatm.o  hyperfine.o paw.o   \
         mkpoints_full.o       charge.o   Lebedev-Laikov.o  stockholder.o dipol.o    pot.o \
         dos.o      elf.o      tet.o      tetweight.o hamil_rot.o \
         chain.o    dyna.o     k-proj.o    sphpro.o    us.o  core_rel.o \
         aedens.o   wavpre.o   wavpre_noio.o broyden.o \
         dynbr.o    hamil_high.o  rmm-diis.o reader.o   writer.o   tutor.o xml_writer.o \
         brent.o    stufak.o   fileio.o   opergrid.o stepver.o  \
         chgloc.o   fast_aug.o fock_multipole.o  fock.o  mkpoints_change.o sym_grad.o \
         mymath.o   internals.o npt_dynamics.o   dynconstr.o dimer_heyden.o dvvtrajectory.o vdwforcefield.o \
         nmr.o      pead.o     subrot.o   subrot_scf.o \
         force.o    pwlhf.o    gw_model.o optreal.o  steep.o    davidson.o  david_inner.o \
         electron.o rot.o  electron_all.o shm.o    pardens.o  paircorrection.o \
         optics.o   constr_cell_relax.o   stm.o    finite_diff.o elpol.o    \
         hamil_lr.o rmm-diis_lr.o  subrot_cluster.o subrot_lr.o \
         lr_helper.o hamil_lrf.o   elinear_response.o ilinear_response.o \
         linear_optics.o \
         setlocalpp.o  wannier.o electron_OEP.o electron_lhf.o twoelectron4o.o \
         mlwf.o     ratpol.o screened_2e.o wave_cacher.o chi_base.o wpot.o \
         local_field.o ump2.o ump2kpar.o fcidump.o ump2no.o \
         bse_te.o bse.o acfdt.o chi.o sydmat.o dmft.o \
         rmm-diis_mlr.o  linear_response_NMR.o wannier_interpol.o linear_response.o getnabij.o

vasp: $(SOURCE) $(FFT3D) $(INC) main.o 
	rm -f vasp
	$(FCL) -o vasp main.o  $(SOURCE)   $(FFT3D) $(LIB) $(LINK)
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
	-rm -f *.g *.o *.L *.mod ; touch *.F

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
mkpoints.o: mkpoints.inc mkpoints.F
wave.o: wave.F
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
# these special rules have been tested for ifc.11 and ifc.12 only

fft3dlib.o : fft3dlib.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
fft3dfurth.o : fft3dfurth.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)
fftw3d.o : fftw3d.F
	$(CPP)
	$(FC) -FR -lowercase -O1 $(INCS) -c $*$(SUFFIX)
fftmpi.o : fftmpi.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)
fftmpiw.o : fftmpiw.F
	$(CPP)
	$(FC) -FR -lowercase -O1 $(INCS) -c $*$(SUFFIX)
wave_high.o : wave_high.F
	$(CPP)
	$(FC) -FR -lowercase -O1 -c $*$(SUFFIX)
# the following rules are probably no longer required (-O3 seems to work)
wave.o : wave.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
paw.o : paw.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
cl_shift.o : cl_shift.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
us.o : us.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
LDApU.o : LDApU.F
	$(CPP)
	$(FC) -FR -lowercase -O2 -c $*$(SUFFIX)
