module mod_met_mtcorf1
implicit none
contains
      subroutine met_mtcorf1( &
                 l,ncin, &
                 dist,mnpar,frac)
!
!***********************************************************************
!
!     ACT
!_A   Correction de la fonction de raccordement de Menter
!_A   l=0  : toutes les parois de tous les domaines
!_A   l>0  : parois du domaine "l"
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    dist       : arg real(ip12      ) ; distance a la paroi
!_I    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_I                                        stockees du point de rattach normale
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    nbdko      : com int              ; nombre de parois a traiter
!_I                                        defini dans "atcaldist"
!_I    lbdko      : com int (mtt       ) ; numero interne des parois
!_I                                        (necessaire pour k-omega)
!
!     I/O
!_/    frac       : arg real(ip12     )  ; fonction de raccord des modeles
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use modeleturb
use mod_atindnor
implicit none
integer :: l
integer :: ncin
double precision :: dist
integer :: mnpar
double precision :: frac
integer :: isdm3
integer :: isens3
integer :: lpar
integer :: m
integer :: m10
integer :: m1max
integer :: m1min
integer :: m20
integer :: m2max
integer :: m2min
integer :: m3
integer :: m30
integer :: m3lim
integer :: m3max
integer :: m3min
integer :: m3mxx
integer :: mfbi
integer :: mpar
integer :: mpb0
integer :: mpn0
integer :: mpn1
integer :: n
integer :: nfac
integer :: nmin
integer :: npar
!
!-----------------------------------------------------------------------
!
      integer dm1,dm2,dm3
      dimension dist(ip12),mnpar(ip12),frac(ip12)
      dimension ncin(ip41)
!
      do mpar=1,nbdko
!       boucle sur les parois
!
        mfbi=lbdko(mpar)
        lpar=ndlb(mfbi)
        if(l.eq.0 .or. l.eq.lpar) then
          mpn0=mpn(mfbi)
          mpn1=mpn0+mmb(mfbi)
          mpb0=mpb(mfbi)
!
!                          mpn0 : pointeur dernier point frontiere precedente
!                                 dans tableaux frontiers a normales stockees
!                          mpn1 : pointeur dernier point de la frontiere
!
!                          nmin : pointeur cellule adjacente a la paroi
!                         m3mxx : dernier point de la normale rattache a la
!                                 paroi
!
!com      calculs des indices des directions parallele et normale a la paroi
          call atindnor( &
                 mfbi, &
                 m10,m20,m30, &
                 m1min,m1max,m2min,m2max,m3min,m3max, &
                 dm1,dm2,dm3,isens3)
!
          isdm3=isens3*dm3
!
          do m=1,mmb(mfbi)
!           boucle sur les points de la frontiere
!
            nfac =mpb0+m
            nmin =ncin(nfac)
!
!
!           F1 force a 1 (Wilcox) pres de la paroi
!
            m3mxx=m3min+isens3*(imxclko+1)
            do m3=m3min,m3mxx,isens3
!             boucle sur les "imxclko+1" points pres de la paroi
!
              n=nmin+(m3-m3min)*dm3*isens3
              frac(n)=1.
            end do
!
!           recherche du dernier maximum a 1 de frac proche de la
!           frontiere externe de la couche limite
            m3lim=0
            do m3=m3min,m3max,isens3
!             boucle depuis la paroi suivant la normale
!
              n=nmin+(m3-m3min)*isdm3
              npar=mnpar(n)
              if(npar.gt.mpn0 .and. npar.le.mpn1) then
!
!               la cellule est rattachee a la frontiere en cours
!               de traitement
!
                m3mxx=m3
                if(frac(n).gt.0.99) then
                  m3lim=m3
                end if
              end if
!             fin de boucle sur la normale a la paroi
            end do
!
!           "frac=1" entre la paroi et "m3lim"
!           "frac=1" decroissante au dela de "m3lim" jusqu'a "m3mxx"
!
            if(isens3.gt.0) then
              if(m3mxx.gt.m3min) then
                do m3=m3min,m3mxx,isens3
                  n=nmin+(m3-m3min)*isdm3
                  if(m3.le.m3lim) then
                    frac(n)=1.
                  else
                    frac(n)=min(frac(n),frac(n-isdm3))
                  end if
!               fin boucle sur la normale
                end do
              end if
            else
              if(m3mxx.lt.m3min) then
                do m3=m3min,m3mxx,isens3
                  n=nmin+(m3-m3min)*isdm3
                  if(m3.ge.m3lim) then
                    frac(n)=1.
                  else
                    frac(n)=min(frac(n),frac(n-isdm3))
                  end if
!               fin boucle sur la normale
                end do
              end if
            end if
!
!           fin de boucle sur les points de la frontiere
          end do
        end if
!
!       fin de boucle sur parois
      end do
!
      return
      end subroutine
end module
