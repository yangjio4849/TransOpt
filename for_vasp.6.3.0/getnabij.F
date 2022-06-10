#include "symbol.inc"
!************************* SUBROUTINE CALC_NABIJ ***********************
! these routines were written by Juergen Furthmueller:
! you can download the required postprocessing routines from
! http://pc06.physik.uni-jena.de/furth/pub/
!
! directory VASP/optics
! please ceck the README and Makefile
!
! You might learn more about the implementation in the following
! article:
! B. Adolph, J. Furthmueller, and F. Bechstedt, PRB 63, 125108 (2001).
!
!***********************************************************************
!
! main driver routine for the calculation of matrix elements of the
! nabla operator ( --> momentum operator --> velocity operator ... )
!
!***********************************************************************

      SUBROUTINE GET_NABIJ(NABIJ,W,WDES,P,KPOINTS,GRID_SOFT,LATT_CUR, &
                            IO,INFO,T_INFO,COMM,IU0,IU)
      USE prec
      USE base
      USE constant
      USE lattice
      USE poscar
      USE pseudo
      USE mkpoints
      USE mgrid
      USE wave
      
      USE mpimy
      USE reader_tags

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin)       W
      TYPE (wavedes)        WDES
      TYPE (type_info)      T_INFO
      TYPE (potcar)         P(T_INFO%NTYP)
      TYPE (kpoints_struct) KPOINTS
      TYPE (grid_3d)        GRID_SOFT
      TYPE (latt)           LATT_CUR
      TYPE (in_struct)      IO
      TYPE (info_struct)    INFO
      TYPE (communic)       COMM

      GDEF  NABIJ(WDES%NB_TOT,WDES%NB_TOT)
      CHARACTER (1) CHARAC
      LOGICAL LDUM

      REAL(kind=8),allocatable::opticyang(:,:,:)
      REAL(kind=8),allocatable::norm(:)
      CHARACTER*12 opticformat


#ifdef MPI
      IONODE  = WDES%COMM%IONODE
      NODE_ME = WDES%COMM%NODE_ME
      IF (WDES%COMM_KINTER%NCPU.NE.1) THEN
         CALL vtutor%alert("CALC_NABIJ: KPAR>1 not implemented. I will&
              & return now to the main routine without calculating any matrix&
              & elements (just returning zeros to save your day ...).")
         NABIJ=0
         RETURN
      END IF
#endif

#if defined(CRAY_MPP) && defined(MPI)
!-----------------------------------------------------------------------
! read in required parameters
!-----------------------------------------------------------------------
! CRAY C90 fails to read INCAR, if all nodes try to read it at the same time
      DO I=1,COMM%NCPU
      IF (COMM%NODE_ME == I) THEN
#endif
! if we cannot do all k-points specify the true number of k-points to be done
      NKOPT=KPOINTS%NKPTS
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'NKOPT', NKOPT, IERR, WRITEXMLINCAR)

! cannot be smaller than NKPTS - readjust (without warning) if invalid input
      NKOPT=MAX(NKOPT,KPOINTS%NKPTS)
! if we cannot do all k-points specify the current k-point counter offset
      NKOFF=0
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'NKOFFOPT', NKOFF, IERR, WRITEXMLINCAR)

! must be a positive number and smaller than the number of k-points to be done
! minus current number of k-points - readjust (without warning) if invalid input
      NKOFF=MAX(NKOFF,0)
      NKOFF=MIN(NKOFF,NKOPT-KPOINTS%NKPTS)
! number of valence bands written on file OPTIC
      NBVAL=(NINT(INFO%NELECT)+1)/2
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'NBVALOPT', NBVAL, IERR, WRITEXMLINCAR)

! must be larger than zero and smaller or equal to the total number of bands
! - readjust (without warning) if invalid input
      NBVAL=MIN(NBVAL,WDES%NB_TOT)
      NBVAL=MAX(NBVAL,1)
! number of conduction bands written of file OPTIC
      NBCON=WDES%NB_TOT-NBVAL
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'NBCONOPT', NBCON, IERR, WRITEXMLINCAR)

! must be larger than zero and smaller or equal to the total number of bands
! - readjust (without warning) if invalid input
      NBCON=MIN(NBCON,WDES%NB_TOT)
      NBCON=MAX(NBCON,1)

#if defined(CRAY_MPP) && defined(MPI)
      ENDIF
      CALL MPI_barrier( COMM%MPI_COMM, ierror )
      ENDDO
#endif
!-----------------------------------------------------------------------
! open files
!-----------------------------------------------------------------------

! we write with direct access, unformatted -> set the record length for OPTIC
      IRECLO=(MAX(NBVAL,2)+1)*IO%ICMPLX
! open OPTIC, write first record containing all data determining the layout
      io_begin
!      IF (IU>=0) THEN
!        IF (IO%LOPEN) THEN
!          OPEN(UNIT=IU,FILE='OPTIC',ACCESS='DIRECT', &
!               FORM='UNFORMATTED',RECL=IRECLO)
!        ELSE
!          OPEN(UNIT=IU,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=IRECLO)
!        ENDIF
!jF:  special trick! WDES%NRSPINORS=2 and ISPIN=1 results in "ISPIN" = -1 !
!     this allows later to detect whether WDES%NRSPINOR was 1 or 2 without the
!     need to write it explicitly onto file OPTIC (preserving record-structure
!     backward compatibility) -- needed because the original version was not
!     yet able to handle spinors and only wrote ISPIN as an output to file OPTIC
!     but no information about WDES%NRSPINORS; so I have changed argument #3
!     (in a way that for WDES%NRSPINORS=1 the output remains unchanged ...)
!        CALL OUTOPT_HEAD(WDES%NB_TOT,NKOPT,(INFO%ISPIN*(3-2*WDES%NRSPINORS)), &
!                         NBVAL,NBCON,IRECLO,IU)
!      ENDIF
      io_end
!-----------------------------------------------------------------------
! calculate < phi_i | nabla | phi_j >
!-----------------------------------------------------------------------

      open(unit=38,file="GROUPVEC",status="REPLACE")
      allocate(opticyang(3,INFO%ISPIN,WDES%NB_TOT))
      allocate(norm(INFO%ISPIN))

! loop over all (current) k-points, spin component and cartesian directions
      write(opticformat,"(A4,I1,A6)")"(I4,",INFO%ISPIN*4,"E18.8)"

      DO NK=1,KPOINTS%NKPTS

       opticyang=0

       DO ISP=1,INFO%ISPIN
        DO IDIR=1,3
! get the matrix element of the nabla-operator in eV/Angstroem units
         CALL NABIJ_SOFT_yang(NABIJ,IDIR,NK,ISP,W,WDES,GRID_SOFT,LATT_CUR)
         CALL NABIJ_AUG_ADD_yang(NABIJ,IDIR,NK,ISP,W,WDES,P,T_INFO)

         NB_TOT=WDES%NB_TOT

         CALLMPI( M_sum_g(WDES%COMM,NABIJ(1,1),NB_TOT*NB_TOT))

         SCALE_NABIJ_TO_VIJ=AUTOA    ! hbar/m_e * 1/length in whatever unit
         NABIJ=NABIJ*SCALE_NABIJ_TO_VIJ   ! rescaling to bohr radii/Hartree
         NABIJ=NABIJ*2*RYTOEV*(0._q,-1._q)
         DO I=1,NB_TOT
!           opticyang(IDIR,ISP,I)=real(NABIJ(I,I)*CONJG(NABIJ(I,I)))
           opticyang(IDIR,ISP,I)=real(NABIJ(I,I))
         ENDDO

        ENDDO    ! IDIR           
       ENDDO     ! ISP

       write(38,"(4F15.8)")(WDES%VKPT(I,NK),I=1,3),WDES%WTKPT(NK)
       DO I=1,WDES%NB_TOT

         norm=0.0
         DO IDIR=1,3
           DO ISP=1,INFO%ISPIN
             norm(ISP)=norm(ISP)+opticyang(IDIR,ISP,I)**2
           ENDDO
         ENDDO
         norm=sqrt(norm)

         write(38,opticformat)I,((opticyang(IDIR,ISP,I),ISP=1,INFO%ISPIN &
         ),IDIR=1,3),(norm(ISP),ISP=1,INFO%ISPIN)
       ENDDO
       write(38,*)""

      ENDDO      ! NK

      close(38)
      io_begin
      IF (IU>=0) CLOSE(IU)
      io_end

      RETURN

      END

!************************* SUBROUTINE NABIJ_SOFT ***********************
!
! calculates < Psi_i | \nabla | Psi_j > for the non-normconserving
! wave functions stored in CPTWFP, augmentation corrections elsewhere
!
!***********************************************************************

      SUBROUTINE NABIJ_SOFT_yang(NABIJ,IDIR,NK,ISP,W,WDES,GRID,LATT_CUR)
      USE prec
      USE constant
      USE lattice
      USE poscar
      USE mgrid
      USE wave
      USE dfast

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR

      GDEF       NABIJ(WDES%NB_TOT,WDES%NB_TOT)
      COMPLEX(q),ALLOCATABLE :: CBLOCK(:,:),GC(:)

      ALLOCATE(CBLOCK(NBLK,WDES%NB_TOT),GC(NBLK))

      NABIJ = zero

      NB_TOT=WDES%NB_TOT
      NPLDIM=WDES%NRPLWV
      NPL=WDES%NPLWKP(NK)

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
      block: DO IBLOCK=0,NPL-1,NBLK
       ILENPL=MIN(NBLK,NPL-IBLOCK)
       IADDPL=MIN(IBLOCK,NPL-1)
       ILENPL=MAX(ILENPL,0)

       DO I=1,ILENPL
        G1=WDES%IGX(I+IADDPL,NK)+WDES%VKPT(1,NK)
        G2=WDES%IGY(I+IADDPL,NK)+WDES%VKPT(2,NK)
        G3=WDES%IGZ(I+IADDPL,NK)+WDES%VKPT(3,NK)
        GC(I)=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))*TPI
       ENDDO
       DO NB=1,NB_TOT
        DO I =1,ILENPL
          CBLOCK(I,NB)=W%CPTWFP(I+IADDPL+ISPINOR*NPL,NB,NK,ISP)*GC(I)*(0._q,-1._q)
        ENDDO
       ENDDO

       CALL GGEMM(trans, 'N', NB_TOT, NB_TOT, m_ ILENPL, one, &
                  CBLOCK(1,1), m_ NBLK, W%CPTWFP(IADDPL+1+ISPINOR*NPL,1,NK,ISP), &
                  m_ NPLDIM, one, NABIJ(1,1), NB_TOT)

      ENDDO block
      ENDDO spinor

      DEALLOCATE(CBLOCK,GC)

      RETURN
      END


!********************* SUBROUTINE NABIJ_AUG_ADD ************************
!
! add augmentation corrections to matrix elements of \nabla
! WARNING: this routine is not optimized at all :-)
!
!***********************************************************************

      SUBROUTINE NABIJ_AUG_ADD_yang(NABIJ,IDIR,NK,ISP,W,WDES,P,T_INFO)
      USE prec
      USE poscar
      USE pseudo
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES

      GDEF       NABIJ(WDES%NB_TOT,WDES%NB_TOT),CSUM

!=======================================================================
! calculate augmentation correction to matrix elements
!=======================================================================

      LMBASE =0
      spinor: DO ISPINOR=0,WDES%NRSPINORS-1 
      NIS=1

      typ: DO NT=1,T_INFO%NTYP
      ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

      NIP=NI_LOCAL(NI, WDES%COMM_INB)  ! local storage index
      IF (NIP==0) CYCLE ion            ! projected wavefunction not on local node

      DO LM =1,P(NT)%LMMAX
      DO LMP=1,P(NT)%LMMAX

       DO NBP=1,WDES%NB_TOT
       DO NB =1,WDES%NB_TOT
         CSUM=GCONJG(W%CPROJ(LMBASE+LM,NB,NK,ISP))* &
              P(NT)%NABLA(IDIR,LM,LMP)* &
              W%CPROJ(LMBASE+LMP,NBP,NK,ISP)
         NABIJ(NB,NBP)=NABIJ(NB,NBP)+CSUM
       ENDDO
       ENDDO

      ENDDO
      ENDDO

      LMBASE=LMBASE+P(NT)%LMMAX

      ENDDO ion
      NIS = NIS+T_INFO%NITYP(NT)
      ENDDO typ
      ENDDO spinor

      RETURN
      END SUBROUTINE

