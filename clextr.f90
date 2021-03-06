module mod_clextr
  implicit none
contains
  subroutine clextr( &
       mfb, &
       ncbd,v, &
       mmb,mpb,ncin, &
       pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    une extrapolation des variables de l'interieur du domaine.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
    integer          ::          m,        mb,       mfb,  mmb(mtt),  mpb(mtt)
    integer          ::         mt,ncbd(ip41),ncin(ip41),        ni,        nl
    double precision ::     cson(ip11),pression(ip11),    temp(ip11),  v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
    mt=mmb(mfb)
!
    do m=1,mt
       mb=mpb(mfb)+m
       nl=ncbd(mb)
       ni=ncin(mb)
       v(nl,1)=v(ni,1)
       v(nl,2)=v(ni,2)
       v(nl,3)=v(ni,3)
       v(nl,4)=v(ni,4)
       v(nl,5)=v(ni,5)
       pression(nl)=pression(ni)
       temp(nl)=temp(ni)
       cson(nl)=cson(ni)
    enddo
!
    return
  end subroutine clextr
end module mod_clextr
