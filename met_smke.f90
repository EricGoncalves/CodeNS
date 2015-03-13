module mod_met_smke
implicit none
contains
      subroutine met_smke( &
                 l, &
                 sn, &
                 vol,v,mu,mut,tprod,cfke, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 t,dtdx,dtdy,dtdz,bark,bare, &
                 qcxts5,qcyts6, &
                 tn1,tn2,tn3, &
                 cson, & 
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_H   DATE_C : Eric Goncalves
!
!     ACT
!_A   Calcul du terme source des equations pour k-epsilon.
!_A   Modele de Jones-Laundeur ou Launder-Sharma.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_L    tprod      : arg real(ip00     )  ; production de k
!_L    dvxx       : arg real(ip00     )  ; grad(V)  vx,x
!_L    dvxy       : arg real(ip00     )  ; grad(V)  vx,y
!_L    dvxz       : arg real(ip00     )  ; grad(V)  vx,z
!_L    dvyx       : arg real(ip00     )  ; grad(V)  vy,x
!_L    dvyy       : arg real(ip00     )  ; grad(V)  vy,y
!_L    dvyz       : arg real(ip00     )  ; grad(V)  vy,z
!_L    dvzx       : arg real(ip00     )  ; grad(V)  vz,x
!_L    dvzy       : arg real(ip00     )  ; grad(V)  vz,y
!_L    dvzz       : arg real(ip00     )  ; grad(V)  vz,z
!_L    bark       : arg real(ip00     )  ; bas Reynolds pour k
!_L    bare       : arg real(ip00     )  ; bas Reynolds pour e
!
!     OUT
!_/    qcxts5     : arg real(ip12     )  ; terme source equation pour k
!_/    qcyts6     : arg real(ip12     )  ; terme source seconde equation
!-O
!     LOC
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  t,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  t,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  t,z
!_L    tn18       : arg real(ip00     )  ; variable de travail
!_L    tn19       : arg real(ip00     )  ; variable de travail
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use chainecarac
use mod_met_smkec
use mod_met_smkesas
use mod_met_bark
use mod_met_bare
use mod_met_smker
use mod_met_smkes
use mod_met_exgr
implicit none
integer :: l
double precision :: sn
double precision :: vol
double precision :: v
double precision :: tprod
double precision :: cfke
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: t
double precision :: dtdx
double precision :: dtdy
double precision :: dtdz
double precision :: bark
double precision :: bare
double precision :: qcxts5
double precision :: qcyts6
double precision :: tn1
double precision :: tn2
double precision :: tn3
double precision :: cson
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
!
!-----------------------------------------------------------------------
!
      double precision mu,mut
      dimension v(ip11,ip60)
      dimension mut(ip12),mu(ip12),qcxts5(ip12),qcyts6(ip12)
      dimension cfke(ip13),vol(ip11),cson(ip11)
      dimension sn(ip31*ndir)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00), &
                dtdx(ip00),dtdy(ip00),dtdz(ip00), &
                t(ip00),tprod(ip00),bark(ip00),bare(ip00), &
                tn1(ip00),tn2(ip00),tn3(ip00)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
!        ----------------------------------------------------------
!com     teq_exgr(grad(v)) --> grad(v) sur les points fictifs
         call met_exgr(l,dvxx,dvxy,dvxz)
         call met_exgr(l,dvyx,dvyy,dvyz)
         call met_exgr(l,dvzx,dvzy,dvzz)
!
!     ---------------------------------------------------------------
!com  met_bark --> calcul de ||grad(racine(k)||**2 (terme bas-reynolds)
!
      call met_bark( &
                 l, &
                 equat, &
                 sn, &
                 vol,v,mu, &
                 t,dtdx,dtdy,dtdz,bark, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!     --------------------------------------------------------------
!com  met_bare --> calcul de  (d2v/dn2)**2 (terme bas-reynolds)
!com               somme des 27 termes (a_ijk)**2 ; a_ijk = grad( grad(v) )
!
      call met_bare( &
                 l, &
                 equat, &
                 sn, &
                 vol,v,mu,mut, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 dtdx,dtdy,dtdz,bare, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!     ----------------------------------------------------------
!com  met_sour --> calcul des termes sources
!
      if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
!       modele de base et avec correction SST
        call met_smkes( &
                 l, &
                 mu,v,cfke, &
                 tprod,bark,bare,qcxts5,qcyts6)
!
      elseif(equatt(4:4).eq.'R') then
!      terme source pour JL realisable
        call met_smker( &
                 l, &
                 v,cfke, &
                 tprod,qcxts5,qcyts6, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
      elseif(equatt(4:4).eq.'C') then      ! Jean Decaix 06/2010
!      terme source pour k-e compressible
        call met_smkec( &
                 l, &
                 mu,v,cfke, &
                 cson, &
                 tprod,bark,bare,qcxts5,qcyts6)
      elseif(equatt(4:4).eq.'L') then      ! Jean Decaix 06/2010 
!      terme source pour k-e SAS
        call met_smkesas( &
                 l, &
                 sn, &
                 vol,mu,v,cfke, &
                 tprod,bark,bare,qcxts5,qcyts6, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 t,dtdx,tn1,tn2,tn3, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
      endif
!
      return
      end subroutine
end module
