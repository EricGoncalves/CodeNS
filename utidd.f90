module mod_utidd
implicit none
contains
      subroutine utidd( &
                 bceqt, &
                 mfl,rpi,rti,d0x,d0y,d0z,ncbd, &
                 mmb,mpb)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme de preparation des donnees pour clidd.
!_A    Il doit remplir les tableaux rpi,rti,d0x,d0y,d0z pour
!_A     direction de la vitesse absolue imposee, parallele a (d0x,d0y,d0z)
!_A     pression d'arret imposee , pi=pa1 * rpi
!_A     enthalpie d'arret imposee ,  hi=ha1 * rti .
!_A
!_A     Utilisation des rapports de pressions et de temperatures fournis
!_A     dans utdon (urpi, urti), et direction imposee parallele a l'axe
!_A     des x (1,0,0).
!
!     INP
!_I    mfl        : arg int              ; numero de la frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    bceqt      : arg real(ip41,neqt ) ; pres d'arret / pres d'arret etat de
!_I                                        ref utilisateur
!
!     OUT
!_O    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_O                                        ref utilisateur a imposer
!_O    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_O                                        ref utilisateur a imposer
!_O    d0x        : arg real(ip40      ) ; composante en x d'une
!_O                                        direction de l'ecoulement
!_O    d0y        : arg real(ip40      ) ; composante en y d'une
!_O                                        direction de l'ecoulement
!_O    d0z        : arg real(ip40      ) ; composante en z d'une
!_O                                        direction de l'ecoulement
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
implicit none
double precision :: bceqt
integer :: mfl
double precision :: rpi
double precision :: rti
double precision :: d0x
double precision :: d0y
double precision :: d0z
integer :: ncbd
integer :: mmb
integer :: mpb
integer :: m
integer :: ml
integer :: mt
!
!-----------------------------------------------------------------------
!
      dimension bceqt(ip41,neqt)
      dimension ncbd(ip41)
      dimension rpi(ip40),rti(ip40),d0x(ip40),d0y(ip40),d0z(ip40)
      dimension mmb(mtt),mpb(mtt)
!
      mt=mmb(mfl)
      do m=1,mt
       ml=mpb(mfl)+m
       rpi(m)=bceqt(ml,1)
       rti(m)=bceqt(ml,2)
       d0x(m)=bceqt(ml,3)
       d0y(m)=bceqt(ml,4)
       d0z(m)=bceqt(ml,5)
      enddo
!
      return
      end subroutine
end module
