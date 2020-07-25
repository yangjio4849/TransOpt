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

MODULE variable 
  use params 
  implicit none
  real(dp),save:: T,efermi
  integer(kind=8),save::tota
  real(dp),save:: scales,vec(3,3),volume
  integer(kind=8) i,j,k
  integer l,eof
  integer ik,iband,ispin,ivec,totale
  integer,save:: nbands,nspins
end MODULE variable

recursive subroutine quick_sort(k,left,right,corresponding)     
  use params 
  implicit none                                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=8),intent(in)::left,right                         !!!!!written by cilly_lan(xinli),2018.03.16                 !!!!!!!!!!!!! 
  real(kind=dp),intent(inout)::k(:,:) !!!k(3,nk)                  !!!!!three dimension quick sort with the z dimension first.!!!!!!!!!!!!!!
  real(kind=dp)::ref(3)     !!!save the reference data            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=8)::intref,intrefvec
  real(kind=dp)::ktmp(3)
  integer(kind=8),intent(inout)::corresponding(:)  !!!save ik for the energy(ispin,iband,corresponding(ik))
  integer(kind=8)::L       !!!l indicates the left element of the area which is been searching,l=1,nk
  integer(kind=8)::r       !!!r indicates the right element of the area which is been searching,r=1,nk
  integer(kind=8)::temp    !!!save a temp corresponding
  if(right<=left)then  !!! only have one or zero element 
    return
  elseif(right-left==1)then  !!!only have two elements
    if(abs(k(1,left)-k(1,right))<1E-5.and.abs(k(2,left)-k(2,right))<1E-5.and.k(3,left)>k(3,right).or.abs(k(1,left)-k(1,right))<1E-5.and.k(2,left)-k(2,right)>1E-5.or.k(1,left)-k(1,right)>1E-5)then
      ktmp=k(:,left);k(:,left)=k(:,right);k(:,right)=ktmp
      temp=corresponding(left);corresponding(left)=corresponding(right);corresponding(right)=temp
    endif
    return
  endif
  L=left;r=right+1
  ktmp=k(:,L);k(:,L)=k(:,(L+r)/2);k(:,(L+r)/2)=ktmp 
  temp=corresponding(L);corresponding(L)=corresponding((L+r)/2);corresponding((L+r)/2)=temp  
  ref=k(:,left)      !!!center element is the reference data
  intref=corresponding(left)
  do
    if(L>=r)then
      exit
    endif
    do L=L+1,r-1   !!!!find an element which is larger than the reference data
      if(abs(k(1,L)-ref(1))<1E-5.and.abs(k(2,L)-ref(2))<1E-5.and.k(3,L)>ref(3).or.abs(k(1,L)-ref(1))<1E-5.and.k(2,L)-ref(2)>1E-5.or.k(1,L)-ref(1)>1E-5)then
        exit
      endif
    enddo
    
    do r=r-1,L,-1  !!!!find an element which is smaller than the reference data
      if(abs(k(1,r)-ref(1))<1E-5.and.abs(k(2,r)-ref(2))<1E-5.and.k(3,r)<ref(3).or.abs(k(1,r)-ref(1))<1E-5.and.k(2,r)-ref(2)<-1E-5.or.k(1,r)-ref(1)<-1E-5)then
        exit
      endif
    enddo
    
    if(L<r)then
    ktmp=k(:,L);k(:,L)=k(:,r);k(:,r)=ktmp
    temp=corresponding(L);corresponding(L)=corresponding(r);corresponding(r)=temp
    endif
  enddo
  
  if(left<r)then  !!!!divided to two areas
    k(:,left)=k(:,r)
    corresponding(left)=corresponding(r)
    k(:,r)=ref
    corresponding(r)=intref
  endif
  
  call quick_sort(k,left,r-1,corresponding)   !!!!sort the left area
  call quick_sort(k,r+1,right,corresponding)  !!!!sort the right area
end subroutine quick_sort
subroutine derivation(Energy,x,y,z,lenth,rtvecb,corresponding,derva,dervb,dervc,dervx,dervy,dervz)
  implicit none
  integer,PARAMETER :: dp=kind(1.0d0)
  real(kind=dp),intent(in)::Energy(:,:,:)
  integer(kind=8),intent(in)::corresponding(:)
  real(dp),intent(in)::lenth(3)
  real(dp),intent(in)::rtvecb(3,3)
  real(dp),intent(out)::derva(:,:,:),dervb(:,:,:),dervc(:,:,:),dervx(:,:,:),dervy(:,:,:),dervz(:,:,:)
  integer(kind=8),intent(in)::x,y,z
  integer(kind=8)::ix,iy,iz,i
  integer::lz,ly,lx,lxx,lyy,lzz 

  do ix=1,x
    do iy=1,y
      do iz=1,z
      i=(ix-1)*y*z+(iy-1)*z+iz  
      lx=0;lxx=1
      ly=0;lyy=1
      lz=0;lzz=1
        if (iz-1==0) lz=1
        if (iy-1==0) ly=1
        if (ix-1==0) lx=1
        if (iz==z) lzz=0
        if (iy==y) lyy=0
        if (ix==x) lxx=0
!!!!!!!vc=(ix-1)*y*z+(iy-1)*z+iz+1  (ix-1)*y*z+(iy-1)*z+iz-1
!!!!!!!vb=(ix-1)*y*z+(iy)*z+iz  (ix-1)*y*z+(iy-2)*z+iz
!!!!!!!va=(ix)*y*z+(iy-1)*z+iz  (ix-2)*y*z+(iy-1)*z+iz
        dervc(:,:,corresponding(i))=(Energy(:,:,corresponding((ix-1)*y*z+(iy-1)*z+iz*lzz+1))-Energy(:,:,corresponding((ix-1)*y*z+(iy-1)*z+iz-1+z*lz)))/2/lenth(3)
        dervb(:,:,corresponding(i))=(Energy(:,:,corresponding((ix-1)*y*z+(iy*lyy)*z+iz))-Energy(:,:,corresponding((ix-1)*y*z+((iy-1+y*ly)-1)*z+iz)))/2/lenth(2)
        derva(:,:,corresponding(i))=(Energy(:,:,corresponding((ix*lxx)*y*z+(iy-1)*z+iz))-Energy(:,:,corresponding(((ix-1+x*lx)-1)*y*z+(iy-1)*z+iz)))/2/lenth(1)
      dervx(:,:,corresponding(i))=derva(:,:,corresponding(i))*rtvecb(1,1)+dervb(:,:,corresponding(i))*rtvecb(2,1)+dervc(:,:,corresponding(i))*rtvecb(3,1)
      dervy(:,:,corresponding(i))=derva(:,:,corresponding(i))*rtvecb(1,2)+dervb(:,:,corresponding(i))*rtvecb(2,2)+dervc(:,:,corresponding(i))*rtvecb(3,2)
      dervz(:,:,corresponding(i))=derva(:,:,corresponding(i))*rtvecb(1,3)+dervb(:,:,corresponding(i))*rtvecb(2,3)+dervc(:,:,corresponding(i))*rtvecb(3,3)
      enddo
    enddo
  enddo
end subroutine derivation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!   writen by cilly_lan(xinli)  2018.Mar.07           !!!!!!!!!!!!!
!!!!!!   get vk efmass                                     !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main 
  use variable
  implicit none
  interface
  recursive subroutine quick_sort(k,left,right,corresponding)
  integer,PARAMETER :: dp=kind(1.0d0)
    real(kind=dp),intent(inout)::k(:,:)
    integer(kind=8),intent(in)::left,right
    integer(kind=8),intent(inout)::corresponding(:)
  end subroutine quick_sort
  subroutine derivation(Energy,x,y,z,lenth,rtvecb,corresponding,derva,dervb,dervc,dervx,dervy,dervz)
  integer,PARAMETER :: dp=kind(1.0d0)
  real(kind=dp),intent(in)::Energy(:,:,:)
  integer(kind=8),intent(in)::corresponding(:),x,y,z
  real(dp),intent(out)::derva(:,:,:),dervb(:,:,:),dervc(:,:,:),dervx(:,:,:),dervy(:,:,:),dervz(:,:,:)
  real(dp),intent(in)::lenth(3)
  real(dp),intent(in)::rtvecb(3,3)
  end subroutine derivation
  end interface
!!!!some counter parameter
integer(kind=8) ix,iy,iz
!!!! parameters in POSCAR part
real(dp) vecscale,vecb(3,3),rtvecb(3,3),volume_per_atom,vecba,vecbb,vecbc,vcomp(3),avev(3),dLa,dLb,dLc,dL
real(dp)::lenth(3),derlenth(3)
!!!! parameters in SYMMETRY part
integer nsym,isym
integer,allocatable::sym(:,:,:)
!!!!OUTCAR part
character*100 efline
!!!!EIGENVAL part
character*100 nomeaning
integer(kind=8) nkpts  !!nbands already defined in module
integer(kind=8) iirk
real(kind=dp),allocatable::irkpt(:,:),eig(:,:,:),eigtemp(:,:),eigeig(:,:,:),kpt(:,:)  !!!!!irkpt and eig is correspond to calculate kpoints
real(kind=dp),allocatable::fulleigcal(:,:,:),fullkptcal(:,:)   !!!!!correspond to the full mesh of calculate kpoints
  !!!!main loop part
integer(kind=8) nkpt
real(kind=dp) err
real(kind=dp) kpttemp(3)
logical new,newx,newy,newz
integer(kind=8),allocatable::corresponding(:)
integer(kind=8) x,y,z
  !!!!derivation part
real(dp),allocatable::va(:,:,:),vb(:,:,:),vc(:,:,:),maa(:,:,:),mab(:,:,:),mac(:,:,:),mba(:,:,:),mbb(:,:,:),mbc(:,:,:),mca(:,:,:),mcb(:,:,:),mcc(:,:,:),mtemp(:,:,:)
real(dp),allocatable::vx(:,:,:),vy(:,:,:),vz(:,:,:),mxx(:,:,:),mxy(:,:,:),mxz(:,:,:),myx(:,:,:),myy(:,:,:),myz(:,:,:),mzx(:,:,:),mzy(:,:,:),mzz(:,:,:)
integer(kind=8),allocatable::veccorres(:),veccorrestemp(:),vecveccorres(:),vecsymm(:),record(:),num(:)

!!!!OUTCAR part
open(unit=41,file="OUTCAR",status="old")
do while (.TRUE.)
  read(41,"(A100)")nomeaning
  if (nomeaning(1:8)==" E-fermi") then
    efline=nomeaning
    exit
  endif
enddo
close(41)

! read in the scale, vec
!
open(unit=39,file="POSCAR",status="old",iostat=eof)
if(eof/=0) stop "no POSCAR"
read(39,*) nomeaning 
read(39,*)vecscale
do ivec=1,3
  read(39,*)vec(:,ivec)
enddo
close(39)

!
! generate vecb from vec & scale
!
vec=vec*vecscale
volume=vec(1,1)*(vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3))+vec(2,1)*(vec(3,2)*vec(1,3)-vec(1,2)*vec(3,3))+vec(3,1)*(vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3))

vecb(1,1)=vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3)
vecb(2,1)=vec(1,3)*vec(3,2)-vec(1,2)*vec(3,3)
vecb(3,1)=vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3)

vecb(1,2)=vec(2,3)*vec(3,1)-vec(3,3)*vec(2,1)
vecb(2,2)=vec(3,3)*vec(1,1)-vec(1,3)*vec(3,1)
vecb(3,2)=vec(1,3)*vec(2,1)-vec(2,3)*vec(1,1)

vecb(1,3)=vec(2,1)*vec(3,2)-vec(2,2)*vec(3,1)
vecb(2,3)=vec(1,2)*vec(3,1)-vec(1,1)*vec(3,2)
vecb(3,3)=vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)

rtvecb=vecb/volume   !!!!!!!!!(i,j),i=colum,j=line
do i=1,3
  dL=sqrt(dot_product(rtvecb(:,i),rtvecb(:,i)))
  rtvecb(:,i)=rtvecb(:,i)/dL
enddo
!rtvecb=transpose(rtvecb)  !!!!!!need a transposed matrix, so now i=line,j=colum
call inverse(rtvecb)

vecb=2*pi*vecb/volume


volume=volume*1E-30

! read SYMMETRY
open(unit=38,file="SYMMETRY")
read(38,*)nsym
allocate(sym(3,3,nsym*2))
do isym=1,nsym
  do j=1,3
    read(38,*)(sym(j,k,isym),k=1,3)
  enddo
    sym(:,:,isym)=transpose(sym(:,:,isym))
    sym(:,:,isym+nsym)=sym(:,:,isym)*(-1.0)
enddo
nsym=nsym*2
close(38)

! read EIGENVAL
open(unit=38,file="EIGENVAL",status="old")
read(38,*)nomeaning,nomeaning,nomeaning,nspins
backspace(38)
read(38,"(A100)")nomeaning
do i=1,4
  read(38,"(A100)") nomeaning
enddo
read(38,"(A100)")nomeaning
read(nomeaning,*,iostat=eof)totale,nkpts,nbands
if (eof/=0) stop "Too many kpoints. Number of electrons and kpoints should be seperated in EIGENVAL (line 6)! "
if (nkpts*nsym<8) stop "kpoints are not enough. Please increasing the Kmesh"
allocate(irkpt(4,nkpts),eig(nspins,nbands,nkpts),eigtemp(nspins,nbands),eigeig(nspins,nbands,nsym*nkpts*2),kpt(3,nsym*nkpts*2),veccorrestemp(nsym*nkpts*2),vecsymm(nsym*nkpts*2),vecveccorres(nsym*nkpts*2)) 
do iirk=1,nkpts
  read(38,*)(irkpt(j,iirk),j=1,4)
  kpt(:,iirk)=irkpt(1:3,iirk)   
  do iband=1,nbands
    read(38,*)nomeaning,eig(:,iband,iirk)
    eigeig(:,iband,iirk)=eig(:,iband,iirk)
  enddo
enddo
close(38)

! main loop  
nkpt=0
kpt=0.0
do isym=1,nsym
  do iirk=1,nkpts
    vecveccorres(iirk)=0
    kpttemp(1)=irkpt(1,iirk)*sym(1,1,isym)+irkpt(2,iirk)*sym(1,2,isym)+irkpt(3,iirk)*sym(1,3,isym)
    kpttemp(2)=irkpt(1,iirk)*sym(2,1,isym)+irkpt(2,iirk)*sym(2,2,isym)+irkpt(3,iirk)*sym(2,3,isym)
    kpttemp(3)=irkpt(1,iirk)*sym(3,1,isym)+irkpt(2,iirk)*sym(3,2,isym)+irkpt(3,iirk)*sym(3,3,isym)

!! set zero to some small coordinate values.
    forall (i=1:3,j=1:nkpts,abs(irkpt(i,j))<1E-5)
      irkpt(i,j)=0;kpt(i,j)=0
    end forall

! Set the range of the symmetrized k point (-0.5,0.5] in the direct coordiate,
! Jiong Yang
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

    new=.TRUE.
    do ik=1,nkpt
      err=abs(kpttemp(1)-kpt(1,ik))+abs(kpttemp(2)-kpt(2,ik))+abs(kpttemp(3)-kpt(3,ik))   
      if (err<=1E-4) then
        new=.FALSE.
        exit
      endif
    enddo
    if (new) then
      nkpt=nkpt+1
      kpt(:,nkpt)=kpttemp(:)
      if (nkpt==1) then
        x=1
        y=1
        z=1
      else
        do ik=1,nkpt-1
          newx=.TRUE.
          if (abs(kpt(1,nkpt)-kpt(1,ik))<=1E-5) then
            newx=.FALSE.
            exit
          endif
        enddo
        if (newx) x=x+1
        do ik=1,nkpt-1
          newy=.TRUE.
          if (abs(kpt(2,nkpt)-kpt(2,ik))<=1E-5) then
            newy=.FALSE.
            exit
          endif
        enddo
        if (newy) y=y+1
        do ik=1,nkpt-1
          newz=.TRUE.
          if (abs(kpt(3,nkpt)-kpt(3,ik))<=1E-5) then
            newz=.FALSE.
            exit
          endif
        enddo
        if (newz) z=z+1
      endif
      eigeig(:,:,nkpt)=eig(:,:,iirk)
      vecsymm(nkpt)=isym
      if (vecveccorres(iirk)==0) then
        veccorrestemp(nkpt)=nkpt
        vecveccorres(iirk)=nkpt
      else
        veccorrestemp(nkpt)=vecveccorres(iirk)
      endif
    endif
  enddo  ! iirk
enddo  ! isym
!deallocate (sym);deallocate(eigtemp)
deallocate(vecveccorres)
allocate (fulleigcal(nspins,nbands,nkpt),fullkptcal(3,nkpt),corresponding(nkpt),veccorres(nkpt))
do i=1,nkpt
  fulleigcal(:,:,i)=eigeig(:,:,i);fullkptcal(:,i)=kpt(:,i);corresponding(i)=i;veccorres(i)=veccorrestemp(i)
end do
deallocate (eigeig);deallocate(kpt);deallocate(veccorrestemp)
call quick_sort(fullkptcal,1,nkpt,corresponding)  !!!!sort the full kmesh of calculate kpoints
allocate(record(nkpt),num(nkpt),vecveccorres(nkpt))
do j=1,nkpt
  record(veccorres(j))=j
  do i=1,nkpt
    if (veccorres(j)==veccorres(i)) num(j)=num(j)+1
  enddo
enddo
do j=1,nkpt
  vecveccorres(corresponding(j))=j
enddo
tota=x*y*z
allocate(va(nspins,nbands,tota),vb(nspins,nbands,tota),vc(nspins,nbands,tota))
allocate(maa(nspins,nbands,tota),mab(nspins,nbands,tota),mac(nspins,nbands,tota))
allocate(mba(nspins,nbands,tota),mbb(nspins,nbands,tota),mbc(nspins,nbands,tota))
allocate(mca(nspins,nbands,tota),mcb(nspins,nbands,tota),mcc(nspins,nbands,tota))
allocate(vx(nspins,nbands,tota),vy(nspins,nbands,tota),vz(nspins,nbands,tota))
allocate(mxx(nspins,nbands,tota),mxy(nspins,nbands,tota),mxz(nspins,nbands,tota))
allocate(myx(nspins,nbands,tota),myy(nspins,nbands,tota),myz(nspins,nbands,tota))
allocate(mzx(nspins,nbands,tota),mzy(nspins,nbands,tota),mzz(nspins,nbands,tota))
allocate(mtemp(nspins,nbands,tota))

vecba=sqrt(DOT_PRODUCT(vecb(:,1),vecb(:,1)))
vecbb=sqrt(DOT_PRODUCT(vecb(:,2),vecb(:,2)))
vecbc=sqrt(DOT_PRODUCT(vecb(:,3),vecb(:,3)))

derlenth(1)=vecba/x
derlenth(2)=vecbb/y
derlenth(3)=vecbc/z

lenth(1)=1.0000000/x
lenth(2)=1.0000000/y
lenth(3)=1.0000000/z

call derivation(fulleigcal,x,y,z,derlenth,rtvecb,corresponding,va,vb,vc,vx,vy,vz)
call derivation(vx,x,y,z,derlenth,rtvecb,corresponding,maa,mab,mac,mxx,mxy,mxz)
call derivation(vy,x,y,z,derlenth,rtvecb,corresponding,mba,mbb,mbc,myx,myy,myz)
call derivation(vz,x,y,z,derlenth,rtvecb,corresponding,mca,mcb,mcc,mzx,mzy,mzz)
call derivation(va,x,y,z,derlenth,rtvecb,corresponding,maa,mab,mac,mtemp,mtemp,mtemp)
call derivation(vb,x,y,z,derlenth,rtvecb,corresponding,mba,mbb,mbc,mtemp,mtemp,mtemp)
call derivation(vc,x,y,z,derlenth,rtvecb,corresponding,mca,mcb,mcc,mtemp,mtemp,mtemp)
vx=vx*e*e*me/hbar/hbar
vy=vy*e*e*me/hbar/hbar
vz=vz*e*e*me/hbar/hbar
deallocate(mtemp)

open(unit=42,file="klist_int",status="replace")
open(unit=47,file="GVEC",status="replace")
open(unit=48,file="EFMASS",status="replace") !!!!!!scale is (a/2pi)^2
do iirk=1,nkpts
  do i=1,nkpt
    if (abs(irkpt(1,iirk)-fullkptcal(1,i))+abs(irkpt(2,iirk)-fullkptcal(2,i))+abs(irkpt(3,iirk)-fullkptcal(3,i))<1E-4) then   
      write(42,"(3F20.14)") fullkptcal(1:3,i)
      write(47,"(3F20.14)") fullkptcal(1:3,i)
      write(48,"(3F20.14)") fullkptcal(1:3,i)
      if (nspins==1) then
        do iband=1,nbands
          write(47,"(I5,5X,4E14.6)")iband,vx(:,iband,corresponding(i)),vy(:,iband,corresponding(i)),vz(:,iband,corresponding(i)),sqrt(vx(:,iband,corresponding(i))**2+vy(:,iband,corresponding(i))**2+vz(:,iband,corresponding(i))**2)
!!          write(47,"(I5,5X,4E14.6)")iband,vx(:,iband,corresponding(i)),vy(:,iband,corresponding(i)),vz(:,iband,corresponding(i)),sqrt(vx(:,iband,corresponding(i))**2+vy(:,iband,corresponding(i))**2+vz(:,iband,corresponding(i))**2)/sqrt(vx(:,iband,veccorres(corresponding(i)))**2+vy(:,iband,veccorres(corresponding(i)))**2+vz(:,iband,veccorres(corresponding(i)))**2)
          write(48,"(I5,5X,9E14.6)")iband,1/mxx(:,iband,corresponding(i)),1/mxy(:,iband,corresponding(i)),1/mxz(:,iband,corresponding(i)),1/myx(:,iband,corresponding(i)),1/myy(:,iband,corresponding(i)),1/myz(:,iband,corresponding(i)),1/mzx(:,iband,corresponding(i)),1/mzy(:,iband,corresponding(i)),1/mzz(:,iband,corresponding(i))
        enddo
        write(47,*)""
        write(48,*)""
      else
        do iband=1,nbands
          write(47,"(I5,5X,8E14.6)")iband,vx(:,iband,corresponding(i)),vy(:,iband,corresponding(i)),vz(:,iband,corresponding(i)),sqrt(vx(:,iband,corresponding(i))**2+vy(:,iband,corresponding(i))**2+vz(:,iband,corresponding(i))**2)
!!          write(47,"(I5,5X,8E14.6)")iband,vx(:,iband,corresponding(i)),vy(:,iband,corresponding(i)),vz(:,iband,corresponding(i)),sqrt(vx(:,iband,corresponding(i))**2+vy(:,iband,corresponding(i))**2+vz(:,iband,corresponding(i))**2)/sqrt(vx(:,iband,veccorres(corresponding(i)))**2+vy(:,iband,veccorres(corresponding(i)))**2+vz(:,iband,veccorres(corresponding(i)))**2)
          write(48,"(I5,5X,18E14.6)")iband,1/mxx(:,iband,corresponding(i)),1/mxy(:,iband,corresponding(i)),1/mxz(:,iband,corresponding(i)),1/myx(:,iband,corresponding(i)),1/myy(:,iband,corresponding(i)),1/myz(:,iband,corresponding(i)),1/mzx(:,iband,corresponding(i)),1/mzy(:,iband,corresponding(i)),1/mzz(:,iband,corresponding(i))
        enddo
        write(47,*)""
        write(48,*)""
      endif
    exit
    endif
  enddo
enddo
close(42)
close(47)
close(48)
 
end

SUBROUTINE error()
implicit none
stop
end subroutine error


subroutine inverse(a)
integer,PARAMETER :: dp=kind(1.0d0)
real(dp)::a(3,3),b(3,3)
real(dp)::debt
debt=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(2,1)*a(3,2)*a(1,3)-a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
b=b/debt
a=b
end subroutine
