module mod_rbvc
  implicit none
contains
  subroutine rbvc(t,ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) sur des facettes frontieres coincidentes
!_A    (continuite par moyennes).
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
!     I/O
!_/    t          : arg real(ip11,ip60 ) ; variables de calcul
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use definition
    use mod_mpi
    implicit none
    integer          ::          m,        mb,        mc,        mf,       mfb
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),ncin(ip41)
    integer          ::         nd,       ndm
    double precision :: t(ip11,ip60),        tper
    double precision,allocatable :: buff(:,:,:,:)
    integer :: req(nbd,2),other,me
!
!-----------------------------------------------------------------------
!
!
    mt=0
    do mf=1,nbd
       mfb=lbd(mf)
       mt=max(mt,mmb(mfb))
    enddo
    allocate(buff(5,mt,nbd,2))

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
!
!     we have to exchange the globally numbered me boundary with the owner of the globally numbered other boundary
!
       tper=protat*real(mper(mfb))
!
       do m=1,mt
          mb=mpb(mfb)+m
          nc=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          buff(1,m,mf,1)=t(nc,1)
          buff(2,m,mf,1)=t(nc,2)
          buff(3,m,mf,1)=t(nc,3)*cos(tper)+t(nc,4)*sin(tper)
          buff(4,m,mf,1)=t(nc,4)*cos(tper)-t(nc,3)*sin(tper)
          buff(5,m,mf,1)=t(nc,5)
!
       enddo
       call MPI_itrans2(buff(:,1:mt,mf,1),bcg_to_proc(me),bcg_to_proc(other),req(mf,1)) ! send
       call MPI_itrans2(buff(:,1:mt,mf,2),bcg_to_proc(other),bcg_to_proc(me),req(mf,2)) ! recv
    enddo

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
!
       call WAIT_MPI(req(mf,2))  ! waiting for the message to be received
!       buff(:,1:mt,mf,2)=buff(:,1:mt,mf,1)!
!
       do m=1,mt
          mb=mpb(mfb)+m
          nd=ncbd(mb)
          ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          t(nd,1) = 0.5*( t(ndm,1)+buff(1,m,mf,2))
          t(nd,2) = 0.5*( t(ndm,2)+buff(2,m,mf,2))
          t(nd,3) = 0.5*( t(ndm,3)+buff(3,m,mf,2))
          t(nd,4) = 0.5*( t(ndm,4)+buff(4,m,mf,2))
          t(nd,5) = 0.5*( t(ndm,5)+buff(5,m,mf,2))
!
       enddo
    enddo

    do mf=1,nbd
       call WAIT_MPI(req(mf,1))  ! waiting for all the messages to be sent
    enddo
    deallocate(buff)
!
    return
  end subroutine rbvc
end module mod_rbvc
