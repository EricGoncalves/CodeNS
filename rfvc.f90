module mod_rfvc
  implicit none
contains
  subroutine rfvc( &
       t,ncbd,mnc, &
       ps,temp,cson,ncin)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) dans des mailles fictives adjacentes a des
!_A    frontieres coincidentes (valeurs dans le domaine coincident).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
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
!_/    rfvc       t          : arg real(ip11,ip60 ) ; variables de calcul
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
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),        nd,ncin(ip41)
    double precision ::   cson(ip11),    ps(ip11),t(ip11,ip60),  temp(ip11),        tper
    double precision,allocatable :: buff(:,:,:,:)
    integer :: req(nbd,2),other,me
!
!-----------------------------------------------------------------------
!
!
    req=MPI_REQUEST_NULL
    mt=0
    do mf=1,nbd
       mfb=lbd(mf)
       mt=max(mt,mmb(mfb))
    enddo
    allocate(buff(8,mt,nbd,2))

    do mf=1,nbd
       mfb=lbd(mf)
       mt =mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
!
!     we have to exchange the globally numbered me boundary with the owner of the globally numbered other boundary
!
       tper=protat*real(mper(mfb))
       do m=1,mt
          mb=mpb(mfb)+m
          nc=ncin(mb)

!       definition des variables aux points fictifs

          buff(1,m,mf,1)= t(nc,1) ! we fill a buffer, so we can send bigger messages simultaneously
          buff(2,m,mf,1)= t(nc,2)
          buff(3,m,mf,1)= t(nc,3)*cos(tper)+t(nc,4)*sin(tper)
          buff(4,m,mf,1)= t(nc,4)*cos(tper)-t(nc,3)*sin(tper)
          buff(5,m,mf,1)= t(nc,5)
          buff(6,m,mf,1)= ps(nc)
          buff(7,m,mf,1)=temp(nc)
          buff(8,m,mf,1)=cson(nc)
       enddo
       call MPI_itrans2(buff(:,1:mt,mf,1),bcg_to_proc(me),bcg_to_proc(other),req(mf,1)) ! send
       call MPI_itrans2(buff(:,1:mt,mf,2),bcg_to_proc(other),bcg_to_proc(me),req(mf,2)) ! recv
        endif
    enddo

    do mf=1,nbd
       mfb=lbd(mf)
       mt =mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
       if (bcg_to_proc(me)/=bcg_to_proc(other)) then
!
       call WAIT_MPI(req(mf,2))  ! waiting for the message to be received
!       buff(:,1:mt,mf,2)=buff(:,1:mt,mf,1)

       do m=1,mt
          mb=mpb(mfb)+m
          nd=ncbd(mb)

!       definition des variables aux points fictifs

          t(nd,1)= buff(1,m,mf,2)
          t(nd,2)= buff(2,m,mf,2)
          t(nd,3)= buff(3,m,mf,2)
          t(nd,4)= buff(4,m,mf,2)
          t(nd,5)= buff(5,m,mf,2)
          ps(nd) = buff(6,m,mf,2)
          temp(nd)=buff(7,m,mf,2)
          cson(nd)=buff(8,m,mf,2)

       enddo
        else
       tper=protat*real(mper(mfb))
       do m=1,mt
          mc=mpc(mfb)+m
          nc=mnc(mc)
          mb=mpb(mfb)+m
          nd=ncbd(mb)
!       definition des variables aux points fictifs
          t(nd,1)= t(nc,1)
          t(nd,2)= t(nc,2)
          t(nd,3)= t(nc,3)*cos(tper)+t(nc,4)*sin(tper)
          t(nd,4)= t(nc,4)*cos(tper)-t(nc,3)*sin(tper)
          t(nd,5)= t(nc,5)
          ps(nd) = ps(nc)
          temp(nd)=temp(nc)
          cson(nd)=cson(nc)
       enddo
        endif
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
  end subroutine rfvc
end module mod_rfvc
