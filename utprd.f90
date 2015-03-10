      subroutine utprd( &
                 bceqt, &
                 mfl,pres,icyc, &
                 mmb,mpb)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme utilisateur de preparation des donnees pour le
!_A    sous-programme clprd.
!_A    Il doit remplir le tableau pres des pressions.
!_A
!_I    mfl        : arg int              ; numero de la frontiere
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!
!     OUT
!_O    pres       : arg real(ip40      ) ; pression statique
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
implicit none
double precision :: bceqt
integer :: mfl
double precision :: pres
integer :: icyc
integer :: mmb
integer :: mpb
integer :: m
integer :: ml
integer :: mt
!
!-----------------------------------------------------------------------
!
      dimension bceqt(ip41,neqt)
      dimension pres(ip40)
      dimension mmb(mtt),mpb(mtt)
!
      mt=mmb(mfl)
      do m=1,mt
       ml=mpb(mfl)+m
       pres(m)=bceqt(ml,1)
      enddo
!
      end
