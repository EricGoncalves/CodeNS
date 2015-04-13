module mod_clpari
  implicit none
contains
  subroutine clpari( &
       mfb,tp, &
       ncbd,v, &
       mmb,mpb,ncin, &
       pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition d'adherence.
!_A    pour une paroi isotherme.
!_A    Derivee de p : nulle suivant les lignes de maillage
!_A                   qui "sortent" de la paroi.
!
!_I    mfb        : arg int              ; numero de frontiere
!_I    tpar       : arg real             ; temp paroi a imposer
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    cv         : com real             ; chal spec a vol cst adim
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    implicit none
    integer          ::          m,        mb,       mfb,  mmb(mtt),  mpb(mtt)
    integer          ::         mt,ncbd(ip41),ncin(ip41),        ni,        nl
    double precision ::     cson(ip11),          epar,pression(ip11),         rhoe1,    temp(ip11)
    double precision ::       tp(ip40),  v(ip11,ip60)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
    mt=mmb(mfb)
!
!$OMP SIMD
    do m=1,mt
       mb=mpb(mfb)+m
       nl=ncbd(mb)
       ni=ncin(mb)
!
       epar=cv*tp(m)
       rhoe1=v(ni,5)-.5*(v(ni,2)**2+v(ni,3)**2+v(ni,4)**2)/v(ni,1)
!
       v(nl,5)=rhoe1
       v(nl,1)=v(nl,5)/epar
       v(nl,2)=0.
       v(nl,3)=0.
       v(nl,4)=0.
!
       pression(nl)=pression(ni)      !hypothese dPdn=0
       temp(nl)=tp(m)
       cson(nl)=sqrt(temp(nl))
    enddo
!
!$OMP END MASTER
    return
  end subroutine clpari
end module mod_clpari
