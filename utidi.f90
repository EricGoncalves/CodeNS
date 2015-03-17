module mod_utidi
  implicit none
contains
  subroutine utidi( &
       bceqt, &
       mfl,rpi,rti,ncbd, &
       mmb,mpb)
!
!***********************************************************************
!
!     ACT
!_A    preparation des donnees pour le sous-programme clidi.
!_A    Il doit remplir les tableaux rpi,rti pour
!_A     pression d'arret imposee , pi=pa1 * rpi
!_A     enthalpie d'arret imposee ,  hi=ha1 * rti
!_A
!_A     Utilisation des rapports de pressions et de temperatures fournis
!_A     dans utdon (urpi, urti).
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
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
    integer          ::    m, mfl,  ml, mmb, mpb
    integer          ::   mt,ncbd
    double precision :: bceqt,  rpi,  rti
!
!-----------------------------------------------------------------------
!
    dimension bceqt(ip41,neqt)
    dimension ncbd(ip41)
    dimension rpi(ip40),rti(ip40)
    dimension mmb(mtt),mpb(mtt)
!
    mt=mmb(mfl)
    do m=1,mt
       ml=mpb(mfl)+m
       rpi(m)=bceqt(ml,1)
       rti(m)=bceqt(ml,2)
    enddo
!
    return
  end subroutine utidi
end module mod_utidi
