module mod_rbte
  implicit none
contains
  subroutine rbte( &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       ncbd,ncin)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des composantes du tenseur visqueux et des flux de chaleur
!_A    sur des facettes frontieres par extrapolation.
!
!
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    implicit none
  integer          ::          m,        mb,        mf,       mfb,        mt
  integer          :: ncbd(ip41),ncin(ip41),        nd,       ndm
  double precision ::  qcx(ip12), qcy(ip12), qcz(ip12),toxx(ip12),toxy(ip12)
  double precision :: toxz(ip12),toyy(ip12),toyz(ip12),tozz(ip12)
!
!-----------------------------------------------------------------------
!

!
    do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
!!$OMP SIMD
       do m=1,mt
          mb=mpb(mfb)+m
          nd=ncbd(mb)
          ndm=ncin(mb)
          toxx(nd)=toxx(ndm)
          toxy(nd)=toxy(ndm)
          toxz(nd)=toxz(ndm)
          toyy(nd)=toyy(ndm)
          toyz(nd)=toyz(ndm)
          tozz(nd)=tozz(ndm)
          qcx(nd)=qcx(ndm)
          qcy(nd)=qcy(ndm)
          qcz(nd)=qcz(ndm)
       enddo
    enddo
!
    return
  end subroutine rbte
end module mod_rbte
