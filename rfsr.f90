module mod_rfsr
implicit none
contains
      subroutine rfsr( &
                 t, &
                 ncbd,mnr,xnr,ynr,znr)
!
!***********************************************************************
!
!     ACT
!_A    Calcul du scalaire (t) dans des mailles fictives adjacentes a des
!_A    frontieres recouvertes (valeurs interpolees dans le domaine recouvrant).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpr        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front recouvertes
!_I    ndrr       : com int (mtb       ) ; numero de domaine recouvrant
!_I    srotr      : com real(mtb       ) ; rotation amenant une front periodique
!_I                                        dans sa position recouverte, sinus
!_I    crotr      : com real(mtb       ) ; rotation amenant une front periodique
!_I                                        dans sa position recouverte, cosinus
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use maillage
implicit none
double precision :: t
integer :: ncbd
integer :: mnr
double precision :: xnr
double precision :: ynr
double precision :: znr
integer :: l
integer :: m
integer :: mf
integer :: mfb
integer :: mfbm
integer :: ml
integer :: mr
integer :: mt
integer :: nd
integer :: nid
integer :: njd
integer :: nr
!
!-----------------------------------------------------------------------
!
      dimension t(ip11)
      dimension ncbd(ip41)
      dimension mnr(ip44),xnr(ip44),ynr(ip44),znr(ip44)
!
      do mf=1,nbd
!
      mfbm=lbd(mf)
      mt=mmb(mfbm)
!
      mfb =mod(mfbm,mtb)
      if (mfb.eq.0) mfb=mtb
      l  = ndrr(mfb)
!
      nid=id2(l)-id1(l)+1
      njd=jd2(l)-jd1(l)+1
!
!DEC$ IVDEP
      do m=1,mt
      mr=mpr(mfbm)+m
      nr=mnr(mr)
      ml=mpb(mfbm)+m
      nd=ncbd(ml)
!
!     definition du scalaire aux points fictifs
!
      t(nd)= &
             (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              )+ &
             (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd)+ &
             (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        )+ &
             (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd)+ &
                xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            )+ &
                xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd)+ &
                xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        )+ &
                xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd)
!
      enddo
      enddo
!
      return
      end subroutine
end module
