module mod_met_iniuttau
implicit none
contains
      subroutine met_iniuttau( &
                 l,mu,mut,v,equat, &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 sn,vol, &
                 ncbd,ncin,mnc, &
                 mnr,xnr,ynr,znr, &
                 nxn,nyn,nzn, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                 ztemp)
!-----------------------------------------------------------------------
!
!     initialisation de utau
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use schemanum
      use boundary
      use proprieteflu
use mod_zfluto
use mod_zgrad2
use mod_rbtr
use mod_zgrad
use mod_rbte
use mod_zvismo
use mod_rbtc
implicit none
integer :: l
double precision :: v
double precision :: txxf5x
double precision :: txyf5y
double precision :: txzf5z
double precision :: tyyf6x
double precision :: tyzf6y
double precision :: tzzf6z
double precision :: qcxts5
double precision :: qcyts6
double precision :: qcz000
double precision :: sn
double precision :: vol
integer :: ncbd
integer :: ncin
integer :: mnc
integer :: mnr
double precision :: xnr
double precision :: ynr
double precision :: znr
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
double precision :: ztemp
integer :: lgsnlt
integer :: mf
integer :: mfc
integer :: mfr
integer :: npsn
!
!-----------------------------------------------------------------------
!
      character(len=7 ) :: equat
      real nxn,nyn,nzn
      real mu,mut
      dimension txxf5x(ip12),txyf5y(ip12),txzf5z(ip12), &
                tyyf6x(ip12),tyzf6y(ip12),tzzf6z(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42), &
                ncbd(ip41),ncin(ip41)
      dimension sn(ip31*ndir),vol(ip11),ztemp(ip11)
      dimension mu(ip12),mut(ip12),mnc(ip43)
      dimension v(ip11,ip60)
      dimension mnr(ip44),xnr(ip44),ynr(ip44),znr(ip44)
      dimension qcxts5(ip12),qcyts6(ip12),qcz000(ip12)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00), &
                dvyz(ip00),dvzx(ip00),dvzy(ip00),dvzz(ip00)
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: dtdx,dtdy,dtdz
      ALLOCATE(dtdx(ip00),dtdy(ip00),dtdz(ip00))
!
       npsn=ndir*npfb(l)+1
       lgsnlt=nnn(l)
!
       call zvismo(l,mu,v,ztemp)
!
          if(ischema.eq.1) then
             call zgrad( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 v,ztemp, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dtdx,dtdy,dtdz)
            else
             call zgrad2( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 v,ztemp, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dtdx,dtdy,dtdz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
            endif
!
         call zfluto( &
                 l,mu,mut,txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dtdx,dtdy,dtdz)
!
      do mf=1,mtnx
         lbd(mf)=nfbn(mf)
      enddo
      nbd=mtnx
      call rbte( &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 ncbd,ncin)
!
      do mf=1,mtax
         lbd(mf)=nfba(mf)
      enddo
      nbd=mtax
            call rbte( &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 ncbd,ncin)
!
      do mfc=1,mtcx
         lbd(mfc)=nfbc(mfc)
      enddo
      nbd=mtcx
            call rbtc( &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 ncbd,ncin,mnc)
!
      do mfr=1,mtrx
         lbd(mfr)=nfbr(mfr)
      enddo
      nbd=mtrx
            call rbtr( &
                 ncbd,mnr,xnr,ynr,znr, &
                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                 qcxts5,qcyts6,qcz000, &
                 ncin)

DEALLOCATE(dtdx,dtdy,dtdz)

      return
      end
end module
