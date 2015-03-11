module mod_lpke
implicit none
contains
      subroutine lpke( &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,l, &
                 mnpar,fgam,ncyc, &
                 tprod,tp, &
                 ztemp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 1999-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Lois de paroi : modele k-epsilon de Jones Launder
!_A    - production de k imposee a la paroi
!_A    - epsilon imposee a la paroi
!_A    2 cas : parois adiabatiques et parois isothermes
!
!     INP
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
use mod_lpke1
use mod_lpke2
implicit none
double precision :: v
double precision :: dist
integer :: ncin
integer :: ncbd
integer :: l
integer :: mnpar
double precision :: fgam
integer :: ncyc
double precision :: tprod
double precision :: tp
double precision :: ztemp
integer :: ldom
integer :: mf
integer :: mfb
integer :: no
!
!-----------------------------------------------------------------------
!
      real mu,mut,nxn,nyn,nzn
!
      dimension mu(ip12),mut(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncin(ip41),ncbd(ip41)
      dimension v(ip11,ip60),dist(ip12)
      dimension mnpar(ip12),fgam(ip42)
      dimension tprod(ip00),tp(ip40)
      dimension ztemp(ip11)
!
      nbd=0
      do no=1,mtbx
!     boucle sur toutes les frontieres
        mfb=nba(no)
        ldom=ndlb(mfb)
        if((cl(mfb)(1:2).eq.'lp').and.(l.eq.ldom)) then
!       la frontiere est une paroi et appartient au domaine en cours de traitement
          nbd=nbd+1
          lbd(nbd)=mfb
        endif
      enddo
!
      do mf=1,nbd
!     boucle sur les frontieres a traiter (parois)
        mfb=lbd(mf)
!
        if(cl(mfb)(1:3).eq.'lp2') then
!       parois adiabatiques
           call lpke1( &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfb,l, &
                 mnpar,fgam,ncyc, &
                 tprod,ztemp)
!
        else if(cl(mfb)(1:3).eq.'lp3') then
!       parois isothermes
          call lpke2( &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfb,l, &
                 mnpar,fgam,ncyc, &
                 tprod,tp,ztemp)
        endif
!
      enddo
!
      return
      end subroutine

end module
