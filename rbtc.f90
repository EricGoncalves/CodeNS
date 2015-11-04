module mod_rbtc
  implicit none
contains
  subroutine rbtc( &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des composantes du tenseur visqueux et des flux de chaleur
!_A    pour des facettes frontieres coincidentes (continuite par moyennes).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux front coinc
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!     OUT
!
!     I/O
!_/    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_/                                        contraintes visqueuses
!_/    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_/                                        contraintes visqueuses
!_/    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_/                                        contraintes visqueuses
!_/    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_/                                        contraintes visqueuses
!_/    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_/                                        contraintes visqueuses
!_/    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_/                                        contraintes visqueuses
!_/    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_/    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_/    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use definition
    use mod_mpi
    implicit none
    integer          ::          m,        mb,        mf,       mfb
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),ncin(ip41)
    integer          ::         nd,       ndm
    double precision ::         cr, qcx(ip12), qcy(ip12)
    double precision ::  qcz(ip12),        sr,toxx(ip12),toxy(ip12)
    double precision :: toxz(ip12),toyy(ip12),toyz(ip12),tozz(ip12)
    double precision,allocatable :: buff(:,:,:,:)
    integer :: req(nbd,2),other,me,bcg_to_mf(num_bcg)
!
!-----------------------------------------------------------------------
!
!
    req=MPI_REQUEST_NULL
    mt=0
    do mf=1,nbd
       mfb=lbd(mf)
       mt=max(mt,mmb(mfb))
       me=bcl_to_bcg(mfb)
       bcg_to_mf(me)=mf
    enddo
    allocate(buff(9,mt,nbd,2))
!
!
    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
!
!     we have to exchange the globally numbered me boundary with the owner of the globally numbered other boundary
!
       sr=-sin(real(mper(mfb))*protat)
       cr= cos(real(mper(mfb))*protat)
!
       do m=1,mt
          mb=mpb(mfb)+m
          nc=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          buff(1,m,mf,1)=toxx(nc)
          buff(2,m,mf,1)=cr*toxy(nc)-sr*toxz(nc)
          buff(3,m,mf,1)=cr*toxz(nc)+sr*toxy(nc)
          buff(4,m,mf,1)=cr*cr*toyy(nc)-2.*cr*sr*toyz(nc)+sr*sr*tozz(nc)
          buff(5,m,mf,1)=-cr*sr*(tozz(nc)-toyy(nc))+(2.*cr*cr-1.)*toyz(nc)
          buff(6,m,mf,1)=sr*sr*toyy(nc)+2.*cr*sr*toyz(nc)+cr*cr*tozz(nc)
          buff(7,m,mf,1)=qcx(nc)
          buff(8,m,mf,1)=cr*qcy(nc)-sr*qcz(nc)
          buff(9,m,mf,1)=cr*qcz(nc)+sr*qcy(nc)
!
       enddo
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
         call MPI_itrans2(buff(:,1:mt,mf,1),bcg_to_proc(me),bcg_to_proc(other),req(mf,1),me) ! send
         call MPI_itrans2(buff(:,1:mt,mf,2),bcg_to_proc(other),bcg_to_proc(me),req(mf,2),other) ! recv
       endif
    enddo

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
         call WAIT_MPI(req(mf,2))  ! waiting for the message to be received
       else
         buff(:,1:mt,mf,2)=buff(:,1:mt,bcg_to_mf(other),1)
       endif
!
       do m=1,mt
          mb =mpb(mfb)+m
          nd =ncbd(mb)
          ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          toxx(nd) = 0.5*( toxx(ndm)+buff(1,m,mf,2))
          toxy(nd) = 0.5*( toxy(ndm)+buff(2,m,mf,2))
          toxz(nd) = 0.5*( toxz(ndm)+buff(3,m,mf,2))
          toyy(nd) = 0.5*( toyy(ndm)+buff(4,m,mf,2))
          toyz(nd) = 0.5*( toyz(ndm)+buff(5,m,mf,2))
          tozz(nd) = 0.5*( tozz(ndm)+buff(6,m,mf,2))
          qcx(nd) = 0.5*(  qcx(ndm)+ buff(7,m,mf,2))
          qcy(nd) = 0.5*(  qcy(ndm)+ buff(8,m,mf,2))
          qcz(nd) = 0.5*(  qcz(ndm)+ buff(9,m,mf,2))
!
       enddo
    enddo

    do mf=1,nbd
       mfb=lbd(mf)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) &
           call WAIT_MPI(req(mf,1))  ! waiting for all the messages to be sent
    enddo
    deallocate(buff)
!
    return
  end subroutine rbtc
end module mod_rbtc
