module mod_clchoc
  implicit none
contains
  subroutine clchoc( &
       mfb,mmb,mpb, &
       ncbd,v)
!
!***********************************************************************
!
!_DA  DATE_C : fevrier 2002 -- AUTEUR : Eric Goncalves / SINUMEF
!
!     ACT
!_A    Calcul des variables a la frontiere apres une onde de choc
!_A    Relation de saut de Rankine-Hugoniot.
!
!     INP
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    mfb        : arg int              ; numero de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use constantes
    use definition
    use proprieteflu
    use schemanum
    implicit none
    integer          ::          m,        mb,       mfb,        ml,  mmb(mtt)
    integer          ::   mpb(mtt),        mt,ncbd(ip41),        nl
    double precision ::          a22,       angle,     cmach22,        gam6,           p
    double precision ::           p2,        rmn2,          ro,         ro2,        thet
    double precision ::           tt,v(ip11,ip60),         v22,          xx
!
!-----------------------------------------------------------------------
!
!
    gam6=gam+1.
    ml=mpb(mfb)+1
!     angle du choc oblique
!      angle=30.8*degrad
!      angle=34.*degrad
    angle=35.*degrad
!
!     rmn2     -> mach normal amont au carre
!     indice 2 -> aval choc
!
    rmn2  =(rm0*sin(angle))**2.
    ro    =roa1/(1.+gam2*rm0**2)**gam4
    tt    =(gam6*rmn2)/(gam1*rmn2+2.)
    ro2   =ro*tt
    p     =pa1 /(1.+gam2*rm0**2)**(gam/gam1)
    p2    =p*(2.*gam*rmn2-gam1)/gam6
    a22   =gam*p2/ro2
    xx    =atan(tan(angle)/tt)
    thet  =-angle+xx
    cmach22=((1.+gam2*rmn2)/(gam*rmn2-gam2))/(sin(xx)**2)
    v22   =a22*cmach22
!
    mt=mmb(mfb)
    do m=1,mt
       mb=mpb(mfb)+m
       nl=ncbd(mb)
       v(nl,1)=ro2
       v(nl,2)=ro2*sqrt(v22)*cos(thet)
       v(nl,3)=ro2*sqrt(v22)*sin(thet)
       v(nl,4)=0.
       v(nl,5)=p2/gam1+0.5*ro2*v22
    enddo
!
    return
  end subroutine clchoc
end module mod_clchoc
