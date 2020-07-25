!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!input files of momentum matrix method:                          !!!!!!!!
!!!!!!!!  GROUPVEC EIGENVAL POSCAR SYMMETRY finale.input control.in     !!!!!!!!
!!!!!!!!input files of derivation method:                               !!!!!!!!
!!!!!!!!  GVEC EIGENVAL POSCAR SYMMETRY finale.input control.in OUTCAR  !!!!!!!!
!!!!!!!!Tips:                                                           !!!!!!!!
!!!!!!!!control.in is unnecessary, if there is no control.in file       !!!!!!!!
!!!!!!!!we will adopt the defalt value of parameters                    !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE params 
  implicit none
  integer,PARAMETER :: dp=kind(1.0d0)
  real(dp),PARAMETER :: tol8= 0.00000001_dp
  real(dp),PARAMETER :: pi=3.141592653589793238462643383279502884197_dp
  real(dp),PARAMETER :: zero=0._dp
  REAL(dp),PARAMETER   :: kb  = 8.6174103E-5_dp     ! eV/K
  REAL(dp),PARAMETER   :: hbar=1.05457148E-34_dp    ! Js
  REAL(dp),PARAMETER   :: e   = 1.602176462E-19_dp  ! ELECTRON CHARGE [C]
  REAL(dp),PARAMETER   :: me     = 9.1093826E-31_dp    ! ELECTRON MASS   [kg]
  REAL(dp),PARAMETER   :: A0_SI        = 0.5291772083E-10_dp ! Bohr radius [m]
  REAL(dp),PARAMETER   :: auprAA       = 0.5291772083_dp     ! Bohr radius [A]
  REAL(dp),PARAMETER   :: auprkg       = 1.66053886E-27_dp   ! atomic mass unit [kg] 
END MODULE params


MODULE Fermifunction
  use params
  implicit none

  CONTAINS 
REAL(dp) FUNCTION DFERMI(E,EF,T,spi)
  REAL(dp) :: e,ef,t
  REAL(dp) :: factor,efact
  integer spi
  factor=(e-ef)/(kb*T)
  if(factor>40) then
    DFERMI=0.0
  else if (factor<-40) then
    DFERMI=0.0
  else
    efact=EXP(factor)
    DFERMI=efact/(kb*T*(1.0+efact)**2)*2/spi   !!!!!!spi=nspin,when nspin=2,smaller0.5,when nspin=1,larger1.
  endif
END FUNCTION DFERMI

REAL(dp) FUNCTION FERMI(E,EF,T,spi)
  REAL(dp) :: e,ef,t
  REAL(dp) :: factor,efact
  integer spi
  factor=(e-ef)/(kb*T)
  if(factor>40) then
    FERMI=0.0
  else if (factor<-40) then
    FERMI=2.0/spi
  else
    efact=EXP(factor)
    FERMI=2.0/spi/(efact+1)
  endif
END FUNCTION FERMI  
END MODULE Fermifunction 

MODULE variable 
  use params 
  implicit none
  real(dp),save:: T,efermi
  real(dp),save,allocatable::irbandenergy(:,:,:),bandenergy(:,:,:),irkwpt(:),irkpt(:,:),kpt(:,:)
  real(dp),save:: scales,vec(3,3),volume
  integer,save,allocatable:: corresponding(:)   ! the correlation with full kmesh and irkmesh
  real(dp) kpttemp(3),kpttempijk(3),err
  integer i,j,k,kk,eof
  logical new
  integer iirk,ik,iband,ispin,ivec,isym
  integer,save:: nirk,nk,totale,nband,nspin,nsym,nss,nsoc
  logical ifsoc,ifderivation,ifprint
end MODULE variable

module ef
use variable
use Fermifunction
implicit none
interface efdetermin
  module procedure getef 
end interface

CONTAINS
SUBROUTINE getef(ferminorigin,fermaxorigin)
implicit none
real(kind=dp),intent(in),optional::ferminorigin,fermaxorigin
real(kind=dp) fermin,fermax
real(kind=dp) de
real(kind=dp) ne(nspin),netemp(nspin)
character*1 nomeaning

de=1.0
fermin=-100.0
fermax=100.0
if (present(ferminorigin).and.present(fermaxorigin)) then
fermin=ferminorigin
fermax=fermaxorigin
endif
do while(de.GE.tol8)
  efermi=(fermin+fermax)/2.0
  ne=0.0
  do iirk=1,nirk
    netemp=0.0
    do iband=1,nband
      do ispin=1,nspin
        netemp(ispin)=netemp(ispin)+FERMI(irbandenergy(ispin,iband,iirk),efermi,T,nss)
      enddo
    enddo
    do ispin=1,nspin
      ne(ispin)=ne(ispin)+netemp(ispin)*irkwpt(iirk)
    enddo
  enddo
  de=abs(sum(ne)-totale)
  if(sum(ne)>totale)fermax=efermi
  if(sum(ne)<totale)fermin=efermi
enddo
return
END SUBROUTINE getef 
END MODULE ef


MODULE variables
use params
implicit none
!! electronic structures related 
real(dp),save,allocatable::irvk(:,:,:,:),vk(:,:,:,:),irtau(:,:,:),alltau(:,:,:)
integer,save,allocatable:: nCBB(:)
real(dp),save,allocatable:: sis(:)        ! scissor for valence band
real(dp),save,allocatable:: sym(:,:,:)       
real(dp),save,allocatable::dos(:),dostemp(:),dosvk(:),dosvktemp(:),dostau(:),dostautemp(:),dosint(:),dosinttemp(:)
real(dp),save::defer       ! -df/de
real(dp),save::fermicur    ! current fermilevel
real(dp),save::carriertotale    ! totale for carrier calculations

!! transport related
logical,save,allocatable::dimens(:)
real(dp),save,allocatable::cond(:,:,:),condinv(:,:,:),cond_tau(:,:,:),cond_tau_inv(:,:,:),mu(:,:,:),mu_tau(:,:,:),kappa(:,:,:),kappa_tau(:,:,:),condtemp(:,:,:),condtemp_tau(:,:,:),mutemp(:,:,:),mutemp_tau(:,:,:),kappatemp(:,:,:),kappatemp_tau(:,:,:),vk2temp(:,:,:),vk2temp_tau(:,:,:),seebeck(:,:,:),seebeck_tau(:,:,:),ke(:,:,:),ke_tau(:,:,:)
real(dp),save:: estart,de,emu,Edef,elastic,sigma
integer,save:: iemu,nemu
logical,save:: CRTAonly

!! others
character*100 nomeaning,input
character*16 opticformat
character*16 bandformat
logical carrier, customizetau

!!for readTAU 2017.11.25
! for reading the TAU_fullmesh and fast calculating everything
!!2018.1.7 print temperature in TAU_fullmesh
!!2019.1 no more TAU_fullmesh, only TAU. All the full mesh info, e.g., alltau,
!can be from irtau and corresponding(nk)
logical alive
logical,save:: readTAU
real(dp),save:: sigma_origin,T_origin,Edef_origin,elastic_origin
!

!!2019.1 for converting band energies and group velocities in irriducible k
!space to the whole k space
logical newk

END MODULE variables


! Calculate the temperature curves for Seebeck, sigma, power factor, etc. 
! get the volume from POSCAR
! get the group velocity square from GROUPVEC
! 2017.11.25 get the tensor of transport parameters from full mesh calculations
! 2019.1.20 use the tensor 3:3 format, and fix the problem of off diagonal terms
program TransOpt
use ef 
use variables
implicit none

allocate(dimens(3))
!!!!!!!!!!!!read control.in
open(unit=37,file="control.in",status="old",iostat=eof)
if(eof/=0) then
  ifderivation=.FALSE.
  dimens=.TRUE.
  ifprint=.FALSE.
else 
  print*,"Attention, you build a control.in file. There are tree lines in this file"
  read(37,*)ifderivation
  read(37,*)dimens(:)
  read(37,*)ifprint
close(37)
endif

open(unit=38,file="EIGENVAL",status="old",iostat=eof)
if (eof/=0) stop "EIGENVAL not found!"
!if nspin==1, no spin; if nspin==2, with spin

if(ifprint) open(unit=39,file="EIGENVAL_full",status="replace")
read(38,*)nomeaning,nomeaning,nomeaning,nspin  !!!!!last is nspin
backspace(38)
read(38,"(A100)")nomeaning
if(ifprint) write(39,"(A100)")nomeaning
do i=1,4
  read(38,"(A100)") nomeaning
  if(ifprint) write(39,"(A100)") nomeaning
enddo
read(38,"(A100)")nomeaning
read(nomeaning,*,iostat=eof) totale,nirk,nband
if (eof/=0) stop "Too many kpoints. Number of electrons and kpoints should be seperated in EIGENVAL (line 6)! "
allocate(irkpt(3,nirk),irvk(3,nspin,nband,nirk),irbandenergy(nspin,nband,nirk),irkwpt(nirk))

do iirk=1,nirk
  read(38,*)irkpt(:,iirk),irkwpt(iirk)
  do iband=1,nband
    read(38,*)nomeaning,(irbandenergy(ispin,iband,iirk),ispin=1,nspin)
  enddo
enddo
close(38)

!! get the volume from real space vec from POSCAR
open(unit=38,file="POSCAR",status="old",iostat=eof)
if (eof/=0) stop "POSCAR not found!"
READ(38,*)nomeaning
READ(38,*)scales
READ(38,*)vec
close(38)
vec=vec*scales
volume=vec(1,1)*(vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3))+vec(2,1)*(vec(3,2)*vec(1,3)-vec(1,2)*vec(3,3))+vec(3,1)*(vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3))
volume=volume*1.0E-30


ifsoc=.false.
nsoc=1
if (ifderivation) then
  open(unit=38,file="OUTCAR",status="old",iostat=eof)
  if(eof/=0) stop "Need OUTCAR to determine whether there is spin-orbit coupling!"
  do while(.TRUE.)
    read(38,*)nomeaning
    if(nomeaning=="LSORBIT") then
      backspace(38)
      read(38,*) nomeaning , nomeaning, ifsoc
      exit 
    endif
  end do
  close(38)

    if (ifsoc) then
      nsoc=2
    else
      nsoc=1
    endif
endif
nss=nspin*nsoc


if (ifderivation) then
  open(unit=38,file="GVEC",status="old",iostat=eof)
  if (eof/=0) stop "GVEC not found!"
else
  print*,"!!!!!!Attention, you can only do the transport calculation without soc!!!!!"
  open(unit=38,file="GROUPVEC",status="old",iostat=eof)
  if (eof/=0) stop "GROUPVEC not found!"
endif

!! get the group velocities of the irriducible k mesh from GROUPVEC
do iirk=1,nirk
  read(38,*) nomeaning
  do iband=1,nband
    read(38,*)nomeaning,((irvk(ivec,ispin,iband,iirk),ispin=1,nspin),ivec=1,3)
  enddo
enddo
close(38)


!! set zero to some small coordinate values.
!forall (i=1:3,j=1:nirk,abs(irkpt(i,j))<1E-5)
!  irkpt(i,j)=0
!end forall


!! get the eigenvalue and vk on full k mesh 
write(*,*)"Get the eigenvalue and vk on full k mesh"


! read SYMMETRY
open(unit=38,file="SYMMETRY")
read(38,*)nsym
allocate(kpt(3,nsym*2*nirk),corresponding(nsym*2*nirk),sym(3,3,nsym*2),vk(3,nspin,nband,nsym*2*nirk),bandenergy(nspin,nband,nsym*2*nirk))
do isym=1,nsym
  do j=1,3
    read(38,*)(sym(j,k,isym),k=1,3)
  enddo
    sym(:,:,isym)=transpose(sym(:,:,isym))
    sym(:,:,nsym+isym)=sym(:,:,isym)*(-1.0)
enddo
nsym=nsym*2
close(38)

if(ifprint) open(unit=40,file="klist",status="replace")
nk=0
kpt=0.0
bandenergy=0.0
vk=0.0
do isym=1,nsym
  do iirk=1,nirk
    kpttemp(1)=irkpt(1,iirk)*sym(1,1,isym)+irkpt(2,iirk)*sym(1,2,isym)+irkpt(3,iirk)*sym(1,3,isym)
    kpttemp(2)=irkpt(1,iirk)*sym(2,1,isym)+irkpt(2,iirk)*sym(2,2,isym)+irkpt(3,iirk)*sym(2,3,isym)
    kpttemp(3)=irkpt(1,iirk)*sym(3,1,isym)+irkpt(2,iirk)*sym(3,2,isym)+irkpt(3,iirk)*sym(3,3,isym)

!! set zero to some small coordinate values.
    forall (i=1:3,abs(kpttemp(i))<1E-5)
      kpttemp(i)=0
    end forall

! Set the range of the symmetrized k point (-0.5,0.5] in the direct coordiate, Jiong Yang
    do i=1,3
      if (kpttemp(i)<=-0.5) then
        do while(.true.)
          kpttemp(i)=kpttemp(i)+1
          if (kpttemp(i)>-0.5) exit
        enddo
      endif
      if (kpttemp(i)>0.5) then
        do while(.true.)
          kpttemp(i)=kpttemp(i)-1
          if (kpttemp(i)<=0.5) exit
        enddo
      endif
    enddo

! Sometimes kpttemp(i) will be slightly larger than -0.5 by taking the above
! symmetry operation, so let them be 0.5 instead
    do i=1,3
      err=abs(-0.5-kpttemp(i))
      if (err<=1E-5) kpttemp(i)=0.5
    enddo

    newk=.TRUE.
    do ik=1,nk
      err=abs(kpttemp(1)-kpt(1,ik))+abs(kpttemp(2)-kpt(2,ik))+abs(kpttemp(3)-kpt(3,ik)) 
      if (err<=1E-4) then
        newk=.FALSE.
        exit
      endif
    enddo

    if (newk) then
      nk=nk+1
      corresponding(nk)=iirk
      kpt(:,nk)=kpttemp(:)
      if(ifprint) write(40,"(6F12.7,9F5.1)") kpt(:,nk),irkpt(:,iirk),(sym(:,i,isym),i=1,3)
      bandenergy(:,:,nk)=irbandenergy(:,:,iirk)
      do iband=1,nband
        do ispin=1,nspin
          call vkrotate(vk(:,ispin,iband,nk),irvk(:,ispin,iband,iirk),sym(:,:,isym),vec)
        enddo
      enddo
    endif
  enddo  ! iirk
enddo  ! isym

if(ifprint) close(40) ! klist


!! Write the eigenvalues and group velocities for full k mesh
if(ifprint) write(39,"(3I7)")totale,nk,nband

write(opticformat,"(A4,I1,A5,A6)")"(I4,",nspin*4,"E18.8",",F9.4)"
write(bandformat,"(A7,I1,A6)")"(I5,2X,",nspin,"F14.6)"
if(ifprint) then
  open(unit=40,file="GROUPVEC_full",status="replace")
  do ik=1,nk
    write(39,*)""
    write(39,"(4E15.7)")kpt(:,ik),1.0/nk
    write(40,"(4F15.8)")kpt(:,ik),1.0/nk
    do iband=1,nband
      write(39,bandformat)iband,bandenergy(:,iband,ik)
      write(40,opticformat)iband,((vk(ivec,ispin,iband,ik),ispin=1,nspin),ivec=1,3),(sqrt(sum(vk(:,ispin,iband,ik)**2)),ispin=1,nspin),sqrt(sum(vk(:,1,iband,ik)**2))/sqrt(sum(irvk(:,1,iband,corresponding(ik))**2))  
    enddo
    write(40,*)""
  enddo

  close(39)
  close(40)
endif
vk=vk*1.0E-10*hbar/e/me


!! for transport calculations
write(*,*)"Starting transport calculations"

allocate(alltau(nspin,nband,nk),irtau(nspin,nband,nirk),vk2temp(3,3,nspin),vk2temp_tau(3,3,nspin))
allocate(cond(3,3,nspin),condinv(3,3,nspin),cond_tau(3,3,nspin),cond_tau_inv(3,3,nspin),condtemp(3,3,nspin),condtemp_tau(3,3,nspin),mu(3,3,nspin),mu_tau(3,3,nspin),mutemp(3,3,nspin),mutemp_tau(3,3,nspin),kappa(3,3,nspin),kappa_tau(3,3,nspin),kappatemp(3,3,nspin),kappatemp_tau(3,3,nspin),seebeck(3,3,nspin),seebeck_tau(3,3,nspin),ke(3,3,nspin),ke_tau(3,3,nspin),dos(nspin),dostemp(nspin),dosint(nspin),dosinttemp(nspin),dosvk(nspin),dostau(nspin),dosvktemp(nspin),dostautemp(nspin),nCBB(nspin),sis(nspin))

!! get the input info from finale.input.
open(unit=38,file="finale.input",status="old",iostat=eof)
if (eof/=0) stop "finale.input not found!"
read(38,"(A80)") input
read(input,*,iostat=eof)efermi,estart,de,nemu,T,nCBB(:),sis(:),carriertotale
carrier=.FALSE.
if(eof==0) carrier=.TRUE.
read(38,*,iostat=eof) Edef,elastic,sigma   !  in the unit eV, GPa
CRTAonly=.false.
if (eof/=0) CRTAonly=.true. !then
!  Edef=10
!  elastic=100
!  sigma=0.5
!endif
close(38)

!! aplly the scissors difference for valence band
do ispin=1,nspin
  bandenergy(ispin,1:nCBB(ispin)-1,:)=bandenergy(ispin,1:nCBB(ispin)-1,:)+sis(ispin)
  irbandenergy(ispin,1:nCBB(ispin)-1,:)=irbandenergy(ispin,1:nCBB(ispin)-1,:)+sis(ispin)
enddo

!!for readTAU 2017.11.25
! for reading the TAU_fullmesh and fast calculating everything
!!2018.1.7 print temperature in TAU_fullmesh
!!2019.1 no more TAU_fullmesh, only TAU. All the full mesh info, e.g., alltau,
!can be from irtau and corresponding(nk)
if (.not.CRTAonly) then
  irtau=0.0
  alltau=0.0
  readTAU=.false.
  inquire(file="TAU",exist=alive)
  if (alive) then
    open(unit=38,file="TAU",status="old")
    read(38,*,iostat=eof)nomeaning,sigma_origin
    close(38)
    if (eof==0.and.nomeaning=='Gauss_smearing='.and.abs(sigma-sigma_origin)<1E-4) readTAU=.true.
  endif
  if(readTAU) then
    write(*,"(A)")"Caution: Existing TAU file has been read. If this is not your intention, please delete TAU and rerun."
    open(unit=38,file="TAU",status="old")
    read(38,*)nomeaning,nomeaning,nomeaning,T_origin,nomeaning,Edef_origin,nomeaning,elastic_origin
    do iirk=1,nirk
      read(38,*)nomeaning
      do iband=1,nband
        read(38,*)nomeaning,irtau(:,iband,iirk)
      enddo
    enddo
    close(38)
    irtau=irtau*1E-14/elastic_origin*Edef_origin*Edef_origin*T_origin*elastic/Edef/Edef/T
  else
    call getTAU()
  endif

!! copy the irtau info to alltau
  do ik=1,nk
    alltau(:,:,ik)=irtau(:,:,corresponding(ik))
  enddo
endif

!
open (unit=38,file="CRTA-trace-e.txt",status="replace")
open (unit=39,file="CRTA-tensor-e.txt",status="replace")
if (.not.CRTAonly) then
  open (unit=138,file="RTA-trace-e.txt",status="replace")
  open (unit=139,file="RTA-tensor-e.txt",status="replace")
endif
if (nspin==1) then
  if (carrier) then
    write(38,"(A10,10A18,A33)")"emu(eV)","n(10^20/cm3)","DOS","intDOS","DOSvk2","|vk|(m/s)","sigma","ke","L(1E-8V2K-2)","S(uV/K)","PF","EFF(10^-20 W^5/3 m s^-1/3 K^-2)"
    write(39,"(A10,7A)")"emu(eV)","  n(10^20/cm3)","  sigma(1,1)-(3,3)","  ke(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1,1)-(3,3)","  EFF(10^-20 W^5/3 m s^-1/3 K^-2)(1,1)-(3,3)"
    if (.not.CRTAonly) then
      write(138,"(A10,10A18)")"emu(eV)","n(10^20/cm3)","DOS","intDOS","tau(1E-14s)","|vk|(m/s)","sigma(S/m)","ke(Wm-1K-1)","L(1E-8V2K-2)","S(uV/K)","PF(1E-4Wm-1K-2)"
      write(139,"(A10,6A)")"emu(eV)","  n(10^20/cm3)","  sigma(S/m)(1,1)-(3,3)","  ke(Wm-1K-1)(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1E-4Wm-1K-2)(1,1)-(3,3)"
    endif
  else
    write(38,"(A10,9A18,A33)")"emu(eV)","DOS","intDOS","DOSvk2","|vk|(m/s)","sigma","ke","L(1E-8V2K-2)","S(uV/K)","PF","EFF(10^-20 W^5/3 m s^-1/3 K^-2)"
    write(39,"(A10,6A)")"emu(eV)","  sigma(1,1)-(3,3)","  ke(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1,1)-(3,3)","  EFF(10^-20 W^5/3 m s^-1/3 K^-2)(1,1)-(3,3)"
    if (.not.CRTAonly) then
      write(138,"(A10,9A18)")"emu(eV)","DOS","intDOS","tau(1E-14s)","|vk|(m/s)","sigma(S/m)","ke(Wm-1K-1)","L(1E-8V2K-2)","S(uV/K)","PF(1E-4Wm-1K-2)"
      write(139,"(A10,5A)")"emu(eV)","  sigma(S/m)(1,1)-(3,3)","  ke(Wm-1K-1)(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1E-4Wm-1K-2)(1,1)-(3,3)"
    endif
  endif
else
  if (carrier) then
    write(38,"(A10,14A18,A33)")"emu(eV)","n(10^20/cm3)","DOS_up","DOS_dn","intDOS_up","intDOS_dn","DOSvk2_up","DOSvk2_dn","|vk|_up(m/s)","|vk|_dn(m/s)","sigma","ke","L(1E-8V2K-2)","S(uV/K)","PF","EFF(10^-20 W^5/3 m s^-1/3 K^-2)"
    write(39,"(A10,7A)")"emu(eV)","  n(10^20/cm3)","  sigma(1,1)-(3,3)","  ke(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1,1)-(3,3)","  EFF(10^-20 W^5/3 m s^-1/3 K^-2)(1,1)-(3,3)"
    if (.not.CRTAonly) then
      write(138,"(A10,14A18)")"emu(eV)","n(10^20/cm3)","DOS_up","DOS_dn","intDOS_up","intDOS_dn","tau_up(1E-14s)","tau_dn(1E-14s)","|vk|_up(m/s)","|vk|_dn(m/s)","sigma(S/m)","ke(Wm-1K-1)","L(1E-8V2K-2)","S(uV/K)","PF(1E-4Wm-1K-2)"
      write(139,"(A10,6A)")"emu(eV)","  n(10^20/cm3)","  sigma(S/m)(1,1)-(3,3)","  ke(Wm-1K-1)(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1E-4Wm-1K-2)(1,1)-(3,3)"
    endif
  else
    write(38,"(A10,13A18,A33)")"emu(eV)","DOS_up","DOS_dn","intDOS_up","intDOS_dn","DOSvk2_up","DOSvk2_dn","|vk|_up(m/s)","|vk|_dn(m/s)","sigma","ke","L(1E-8V2K-2)","S(uV/K)","PF","EFF(10^-20 W^5/3 m s^-1/3 K^-2)"
    write(39,"(A10,6A)")"emu(eV)","  sigma(1,1)-(3,3)","  ke(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1,1)-(3,3)","  EFF(10^-20 W^5/3 m s^-1/3 K^-2)(1,1)-(3,3)"
    if (.not.CRTAonly) then
      write(138,"(A10,13A18)")"emu(eV)","DOS_up","DOS_dn","intDOS_up","intDOS_dn","tau_up(1E-14s)","tau_dn(1E-14s)","|vk|_up(m/s)","|vk|_dn(m/s)","sigma(S/m)","ke(Wm-1K-1)","L(1E-8V2K-2)","S(uV/K)","PF(1E-4Wm-1K-2)"
      write(139,"(A10,5A)")"emu(eV)","  sigma(S/m)(1,1)-(3,3)","  ke(Wm-1K-1)(1,1)-(3,3)","  L(1E-8V2K-2)(1,1)-(3,3)","  S(uV/K)(1,1)-(3,3)","  PF(1E-4Wm-1K-2)(1,1)-(3,3)"
    endif
  endif
endif
close(38)
close(39)
if (.not.CRTAonly) then
  close(138)
  close(139)
endif

!!main loop
do iemu=1,nemu
  emu=estart+(iemu-1)*de
  fermicur=efermi+emu
  call gettrans()
  call getdos()
  call output()
enddo
end program TransOpt 


SUBROUTINE getTAU()
use variables
use variable
!$ use omp_lib
implicit none
real(dp) gauss,C
integer ikprime,ibandprime
real(dp) deltaenergy

irtau=0

!!call OMP_SET_NUM_THREADS(24)
!$OMP parallel do
do iirk=1,nirk
  do iband=1,nband
    do ispin=1,nspin
      do ikprime=1,nirk
        do ibandprime=1,nband
! assume deltaenergy> 5 to set gauss 0
          deltaenergy=abs((irbandenergy(ispin,ibandprime,ikprime)-irbandenergy(ispin,iband,iirk))/sigma)
          if(deltaenergy > 5.0) then
            gauss=0.0
          else
            gauss=exp(deltaenergy*deltaenergy/(-2.0))/sqrt(2*pi)/(sigma*e)*irkwpt(ikprime)
          endif
          irtau(ispin,iband,iirk)=irtau(ispin,iband,iirk)+gauss
        enddo
      enddo
    enddo
  enddo
enddo
!$OMP end parallel do

C=2*pi*kb*Edef*Edef*e*e*e/(elastic*1E9)/hbar

irtau=1/(irtau*C*T/volume)

open(unit=38,file="TAU",status="replace")

!!for readTAU 2017.11.25
! for reading the TAU_fullmesh and fast calculating everything
!!2018.1.7 print temperature in TAU_fullmesh
!!2019.1 no more TAU_fullmesh, only TAU. All the full mesh info, e.g., alltau,
!can be from irtau and corresponding(nk)
write(38,"(A,F10.4,A,F10.4,A,F10.4,A,F10.4)") "Gauss_smearing= ",sigma,"   T(K)= ",T, "   Edef(eV)= ",Edef,"   Elastic(GPa)= ",elastic
write(38,*)""
!

do iirk=1,nirk
  write(38,"(4F15.8)")irkpt(:,iirk),irkwpt(iirk)
  do iband=1,nband
    if (nspin==1) then
      write(38,"(I5,F12.6,A)")iband,irtau(:,iband,iirk)*1E14, "  *10^-14 s"
    else
      write(38,"(I5,2F12.6,A)")iband,irtau(:,iband,iirk)*1E14, "  *10^-14 s"
    endif
  enddo
  write(38,*)" "
enddo
close(38)

END SUBROUTINE getTAU


SUBROUTINE getdos()
use variables
use variable
use Fermifunction
implicit none
dos=0.0
dosvk=0.0
dostau=0.0
dosint=0.0
do ik=1,nk
  dostemp=0.0
  dosvktemp=0.0
  dostautemp=0.0
  dosinttemp=0.0
  do iband=1,nband
    do ispin=1,nspin
      defer=DFERMI(bandenergy(ispin,iband,ik),fermicur,T,nss)
      dostemp(ispin)=dostemp(ispin)+defer
      dosvktemp(ispin)=dosvktemp(ispin)+defer*sum(vk(:,ispin,iband,ik)**2)
      if (.not.CRTAonly) dostautemp(ispin)=dostautemp(ispin)+defer*alltau(ispin,iband,ik)
      dosinttemp(ispin)=dosinttemp(ispin)+FERMI(bandenergy(ispin,iband,ik),fermicur,T,nss)
    enddo
  enddo
  do ispin=1,nspin
    dos(ispin)=dos(ispin)+dostemp(ispin)/(1.0*nk)
    dosvk(ispin)=dosvk(ispin)+dosvktemp(ispin)/(1.0*nk)
    if (.not.CRTAonly) dostau(ispin)=dostau(ispin)+dostautemp(ispin)/(1.0*nk)
    dosint(ispin)=dosint(ispin)+dosinttemp(ispin)/(1.0*nk)
  enddo
enddo
END SUBROUTINE getdos 

! 2019.1.20  fix the tensor problem
SUBROUTINE gettrans()
use variables
use variable
use Fermifunction
implicit none
cond=0.0
mu=0.0
kappa=0.0
if (.not.CRTAonly) cond_tau=0.0
if (.not.CRTAonly) mu_tau=0.0
if (.not.CRTAonly) kappa_tau=0.0
do ik=1,nk
  condtemp=0.0
  mutemp=0.0
  kappatemp=0.0
  if (.not.CRTAonly) condtemp_tau=0.0
  if (.not.CRTAonly) mutemp_tau=0.0
  if (.not.CRTAonly) kappatemp_tau=0.0
  do iband=1,nband
    do ispin=1,nspin
      defer=DFERMI(bandenergy(ispin,iband,ik),fermicur,T,nss)/e
      do i=1,3
        do j=1,3
          vk2temp(i,j,ispin)=vk(i,ispin,iband,ik)*vk(j,ispin,iband,ik)*defer*e*e
          condtemp(i,j,ispin)=condtemp(i,j,ispin)+vk2temp(i,j,ispin)
          mutemp(i,j,ispin)=mutemp(i,j,ispin)+vk2temp(i,j,ispin)*(bandenergy(ispin,iband,ik)-fermicur)*e
          kappatemp(i,j,ispin)=kappatemp(i,j,ispin)+vk2temp(i,j,ispin)*(bandenergy(ispin,iband,ik)-fermicur)**2*e*e
          if (.not.CRTAonly) vk2temp_tau(i,j,ispin)=vk2temp(i,j,ispin)*alltau(ispin,iband,ik)
          if (.not.CRTAonly) condtemp_tau(i,j,ispin)=condtemp_tau(i,j,ispin)+vk2temp_tau(i,j,ispin)
          if (.not.CRTAonly) mutemp_tau(i,j,ispin)=mutemp_tau(i,j,ispin)+vk2temp_tau(i,j,ispin)*(bandenergy(ispin,iband,ik)-fermicur)*e
          if (.not.CRTAonly) kappatemp_tau(i,j,ispin)=kappatemp_tau(i,j,ispin)+vk2temp_tau(i,j,ispin)*(bandenergy(ispin,iband,ik)-fermicur)**2*e*e
        enddo
      enddo
    enddo
  enddo
  do ispin=1,nspin
    do i=1,3
      do j=1,3
        cond(i,j,ispin)=cond(i,j,ispin)+condtemp(i,j,ispin)/(1.0*nk)
        mu(i,j,ispin)=mu(i,j,ispin)+mutemp(i,j,ispin)/(1.0*nk)
        kappa(i,j,ispin)=kappa(i,j,ispin)+kappatemp(i,j,ispin)/(1.0*nk)
        if (.not.CRTAonly) cond_tau(i,j,ispin)=cond_tau(i,j,ispin)+condtemp_tau(i,j,ispin)/(1.0*nk)
        if (.not.CRTAonly) mu_tau(i,j,ispin)=mu_tau(i,j,ispin)+mutemp_tau(i,j,ispin)/(1.0*nk)
        if (.not.CRTAonly) kappa_tau(i,j,ispin)=kappa_tau(i,j,ispin)+kappatemp_tau(i,j,ispin)/(1.0*nk)
      enddo
    enddo
  enddo
enddo

cond=cond/volume
mu=mu/(-1*e*T*volume)
kappa=kappa/(e*e*T*volume)
if (.not.CRTAonly) cond_tau=cond_tau/volume
if (.not.CRTAonly) mu_tau=mu_tau/(-1*e*T*volume)
if (.not.CRTAonly) kappa_tau=kappa_tau/(e*e*T*volume)
do i=1,3
  if (.NOT.dimens(i)) then 
    do j=1,3
      do ispin=1,nspin
      cond(i,j,ispin)=1E-30
      cond(j,i,ispin)=1E-30
      if (.not.CRTAonly) cond_tau(i,j,ispin)=1E-30
      if (.not.CRTAonly) cond_tau(j,i,ispin)=1E-30
      enddo
    enddo
  endif
enddo

do ispin=1,nspin
  call inverse(cond(:,:,ispin),condinv(:,:,ispin))
  condinv(:,:,ispin)=transpose(condinv(:,:,ispin))
  if (.not.CRTAonly) then
    call inverse(cond_tau(:,:,ispin),cond_tau_inv(:,:,ispin))
    cond_tau_inv(:,:,ispin)=transpose(cond_tau_inv(:,:,ispin))
  endif
enddo

!write(*,*)cond,condinv,mu

seebeck=0.0
if (.not.CRTAonly) seebeck_tau=0.0
ke=0.0
if (.not.CRTAonly) ke_tau=0.0
do ispin=1,nspin
  do i=1,3
    do j=1,3
!      do k=1,3
!        seebeck(i,j,ispin)=seebeck(i,j,ispin)-condinv(k,i,ispin)*mu(k,j,ispin)
!        if (.not.CRTAonly) seebeck_tau(i,j,ispin)=seebeck_tau(i,j,ispin)-cond_tau_inv(k,i,ispin)*mu_tau(k,j,ispin)
!      enddo
      seebeck(i,j,ispin)=condinv(i,j,ispin)*mu(i,j,ispin)
      if (.not.CRTAonly) seebeck_tau(i,j,ispin)=cond_tau_inv(i,j,ispin)*mu_tau(i,j,ispin)
      ke(i,j,ispin)=kappa(i,j,ispin)-T*mu(i,j,ispin)**2*condinv(i,j,ispin)
      if (.not.CRTAonly) ke_tau(i,j,ispin)=kappa_tau(i,j,ispin)-T*mu_tau(i,j,ispin)**2*cond_tau_inv(i,j,ispin)
    enddo
  enddo
enddo

END SUBROUTINE gettrans


SUBROUTINE output()
use variables
use variable

implicit none
real(dp) outputseebeck_tensor(3,3),outputseebeck_tautensor(3,3),outputcond_tensor(3,3),outputcond_tautensor(3,3),outputke_tensor(3,3),outputke_tautensor(3,3)
real(dp) outputseebeck_trace,outputseebeck_tautrace,outputcond_trace,outputcond_tautrace,outputke_trace,outputke_tautrace

do i=1,3
  do j=1,3
    outputcond_tensor(i,j)=sum(cond(i,j,:))
    outputke_tensor(i,j)=sum(ke(i,j,:))
    outputseebeck_tensor(i,j)=(seebeck(i,j,1)*cond(i,j,1)+seebeck(i,j,nspin)*cond(i,j,nspin))/(cond(i,j,1)+cond(i,j,nspin))
    if (.not.CRTAonly) then
      outputcond_tautensor(i,j)=sum(cond_tau(i,j,:))
      outputke_tautensor(i,j)=sum(ke_tau(i,j,:))
      outputseebeck_tautensor(i,j)=(seebeck_tau(i,j,1)*cond_tau(i,j,1)+seebeck_tau(i,j,nspin)*cond_tau(i,j,nspin))/(cond_tau(i,j,1)+cond_tau(i,j,nspin))
    endif
  enddo
enddo

outputcond_trace=(outputcond_tensor(1,1)+outputcond_tensor(2,2)+outputcond_tensor(3,3))/3.0
outputke_trace=(outputke_tensor(1,1)+outputke_tensor(2,2)+outputke_tensor(3,3))/3.0
outputseebeck_trace=(outputseebeck_tensor(1,1)+outputseebeck_tensor(2,2)+outputseebeck_tensor(3,3))/3.0
if (.not.CRTAonly) then
outputcond_tautrace=(outputcond_tautensor(1,1)+outputcond_tautensor(2,2)+outputcond_tautensor(3,3))/3.0
outputke_tautrace=(outputke_tautensor(1,1)+outputke_tautensor(2,2)+outputke_tautensor(3,3))/3.0
outputseebeck_tautrace=(outputseebeck_tautensor(1,1)+outputseebeck_tautensor(2,2)+outputseebeck_tautensor(3,3))/3.0
endif

open (unit=38,file="CRTA-trace-e.txt",status="old",position="append")
open (unit=39,file="CRTA-tensor-e.txt",status="old",position="append")
if (.not.CRTAonly) then
  open (unit=138,file="RTA-trace-e.txt",status="old",position="append")
  open (unit=139,file="RTA-tensor-e.txt",status="old",position="append")
endif

if (nspin==1) then
  if (carrier) then
    write(38,"(F10.5,11G18.8)")emu,(sum(dosint)-carriertotale)/volume/1.0E26,dos(:),dosint(:),dosvk(:),sqrt(dosvk(:)/dos(:)),outputcond_trace,outputke_trace,outputke_trace/outputcond_trace/T*1E8,outputseebeck_trace*1E6,outputseebeck_trace**2*outputcond_trace/1E10,outputseebeck_trace**2*outputcond_trace/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    write(39,"(F10.5,G18.8,54G)") emu,(sum(dosint)-carriertotale)/volume/1.0E26,outputcond_tensor,outputke_tensor,outputke_tensor/outputcond_tensor/T*1E8,outputseebeck_tensor*1E6,outputseebeck_tensor**2*outputcond_tensor/1E10,outputseebeck_tensor**2*outputcond_tensor/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    if (.not.CRTAonly) then
      write(138,"(F10.5,10G18.8)")emu,(sum(dosint)-carriertotale)/volume/1.0E26,dos(:),dosint(:),dostau(:)/dos(:)*1E14,sqrt(dosvk(:)/dos(:)),outputcond_tautrace,outputke_tautrace,outputke_tautrace/outputcond_tautrace/T*1E8,outputseebeck_tautrace*1E6,outputseebeck_tautrace**2*outputcond_tautrace*1E4
      write(139,"(F10.5,G18.8,45G)") emu,(sum(dosint)-carriertotale)/volume/1.0E26,outputcond_tautensor,outputke_tautensor,outputke_tautensor/outputcond_tautensor/T*1E8,outputseebeck_tautensor*1E6,outputseebeck_tautensor**2*outputcond_tautensor*1E4
    endif
  else
    write(38,"(F10.5,10G18.8)")emu,dos(:),dosint(:),dosvk(:),sqrt(dosvk(:)/dos(:)),outputcond_trace,outputke_trace,outputke_trace/outputcond_trace/T*1E8,outputseebeck_trace*1E6,outputseebeck_trace**2*outputcond_trace/1E10,outputseebeck_trace**2*outputcond_trace/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    write(39,"(F10.5,54G)") emu,outputcond_tensor,outputke_tensor,outputke_tensor/outputcond_tensor/T*1E8,outputseebeck_tensor*1E6,outputseebeck_tensor**2*outputcond_tensor/1E10,outputseebeck_tensor**2*outputcond_tensor/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    if (.not.CRTAonly) then
      write(138,"(F10.5,9G18.8)")emu,dos(:),dosint(:),dostau(:)/dos(:)*1E14,sqrt(dosvk(:)/dos(:)),outputcond_tautrace,outputke_tautrace,outputke_tautrace/outputcond_tautrace/T*1E8,outputseebeck_tautrace*1E6,outputseebeck_tautrace**2*outputcond_tautrace*1E4
      write(139,"(F10.5,45G)") emu,outputcond_tautensor,outputke_tautensor,outputke_tautensor/outputcond_tautensor/T*1E8,outputseebeck_tautensor*1E6,outputseebeck_tautensor**2*outputcond_tautensor*1E4
    endif
  endif

else

  if (carrier) then
    write(38,"(F10.5,15G18.8)")emu,(sum(dosint)-carriertotale)/volume/1.0E26,dos(:),dosint(:),dosvk(:),sqrt(dosvk(:)/dos(:)),outputcond_trace,outputke_trace,outputke_trace/outputcond_trace/T*1E8,outputseebeck_trace*1E6,outputseebeck_trace**2*outputcond_trace/1E10,outputseebeck_trace**2*outputcond_trace/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    write(39,"(F10.5,G18.8,54G)") emu,(sum(dosint)-carriertotale)/volume/1.0E26,outputcond_tensor,outputke_tensor,outputke_tensor/outputcond_tensor/T*1E8,outputseebeck_tensor*1E6,outputseebeck_tensor**2*outputcond_tensor/1E10,outputseebeck_tensor**2*outputcond_tensor/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    if (.not.CRTAonly) then
      write(138,"(F10.5,14G18.8)")emu,(sum(dosint)-carriertotale)/volume/1.0E26,dos(:),dosint(:),dostau(:)/dos(:)*1E14,sqrt(dosvk(:)/dos(:)),outputcond_tautrace,outputke_tautrace,outputke_tautrace/outputcond_tautrace/T*1E8,outputseebeck_tautrace*1E6,outputseebeck_tautrace**2*outputcond_tautrace*1E4
    write(139,"(F10.5,G18.8,45G)") emu,(sum(dosint)-carriertotale)/volume/1.0E26,outputcond_tautensor,outputke_tautensor,outputke_tautensor/outputcond_tautensor/T*1E8,outputseebeck_tautensor*1E6,outputseebeck_tautensor**2*outputcond_tautensor*1E4
    endif
  else
    write(38,"(F10.5,14G18.8)")emu,dos(:),dosint(:),dosvk(:),sqrt(dosvk(:)/dos(:)),outputcond_trace,outputke_trace,outputke_trace/outputcond_trace/T*1E8,outputseebeck_trace*1E6,outputseebeck_trace**2*outputcond_trace/1E10,outputseebeck_trace**2*outputcond_trace/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    write(39,"(F10.5,54G)") emu,outputcond_tensor,outputke_tensor,outputke_tensor/outputcond_tensor/T*1E8,outputseebeck_tensor*1E6,outputseebeck_tensor**2*outputcond_tensor/1E10,outputseebeck_tensor**2*outputcond_tensor/(sum(dos)/e/volume)**(2.0/3.0)*1E20
    if (.not.CRTAonly) then
      write(138,"(F10.5,13G18.8)")emu,dos(:),dosint(:),dostau(:)/dos(:)*1E14,sqrt(dosvk(:)/dos(:)),outputcond_tautrace,outputke_tautrace,outputke_tautrace/outputcond_tautrace/T*1E8,outputseebeck_tautrace*1E6,outputseebeck_tautrace**2*outputcond_tautrace*1E4
      write(139,"(F10.5,45G)") emu,outputcond_tautensor,outputke_tautensor,outputke_tautensor/outputcond_tautensor/T*1E8,outputseebeck_tautensor*1E6,outputseebeck_tautensor**2*outputcond_tautensor*1E4
    endif
  endif
endif
close (38)
close (39)
if (.not.CRTAonly) then
  close (138)
  close (139)
endif
END SUBROUTINE output


subroutine vkrotate(finalvk,originalvk,sym,vec)
real(kind=8)::finalvk(3),originalvk(3),sym(3,3),vec(3,3)
real(kind=8)::vecb(3,3),symvecb(3,3),vecsymvecb(3,3)
call inverse(transpose(vec),vecb)
call matmulthreethree(transpose(sym),vecb,symvecb)
call matmulthreethree(vec,symvecb,vecsymvecb)
call matmulonethree(originalvk,vecsymvecb,finalvk)
end subroutine 


subroutine inverse(a,b)
real(kind=8)::a(3,3),b(3,3)
real(kind=8)::debt
debt=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(2,1)*a(3,2)*a(1,3)-a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
b(2,1)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
b(3,1)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
b(1,2)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
b(3,2)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
b(1,3)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
b(2,3)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
b=b/debt
end subroutine

subroutine matmulonethree(a,b,c)
real(kind=8)::a(3),b(3,3),c(3)
c(1)=a(1)*b(1,1)+a(2)*b(2,1)+a(3)*b(3,1)
c(2)=a(1)*b(1,2)+a(2)*b(2,2)+a(3)*b(3,2)
c(3)=a(1)*b(1,3)+a(2)*b(2,3)+a(3)*b(3,3)
end subroutine

subroutine matmulthreethree(a,b,c)
real(kind=8)::a(3,3),b(3,3),c(3,3)
c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
c(1,3)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)+a(2,3)*b(3,1)
c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
c(2,3)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
c(3,1)=a(3,1)*b(1,1)+a(3,2)*b(2,1)+a(3,3)*b(3,1)
c(3,2)=a(3,1)*b(1,2)+a(3,2)*b(2,2)+a(3,3)*b(3,2)
c(3,3)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
end subroutine
