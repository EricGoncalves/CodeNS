module mod_rfsc
  implicit none
contains
  subroutine rfsc(t,ncbd,mnc)
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
!_/    t          : arg real(ip11      ) ; variables de calcul
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use mod_mpi
    implicit none
    integer          ::          m,        mb,        mc,        mf,       mfb
    integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),        nd
    double precision :: t(ip11)
    double precision,allocatable :: buff(:,:,:)
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
    allocate(buff(mt,nbd,2))


    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       other=ndcc(mfb)
       me=bcl_to_bcg(mfb)
!
!     we have to exchange the globally numbered me boundary with the owner of the globally numbered other boundary
!
       do m=1,mt
          mc=mpc(mfb)+m
          nc=mnc(mc)
!
!     definition d'un scalaire aux points fictifs
!
           buff(m,mf,1)=t(nc) ! we fill a buffer, so we can send bigger messages simultaneously
!
       enddo
       call MPI_itrans2(buff(1:mt,mf,1),bc_to_proc(me),bc_to_proc(other),req(mf,1)) ! send
       call MPI_itrans2(buff(1:mt,mf,2),bc_to_proc(other),bc_to_proc(me),req(mf,2)) ! recv
       
    enddo

    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
!
       call WAIT_MPI(req(mf,2))  ! waiting for the message to be received
!       buff(1:mt,mf,2)=buff(1:mt,mf,1)

       do m=1,mt
          mb=mpb(mfb)+m
          nd=ncbd(mb)
!
!     definition d'un scalaire aux points fictifs
!
          t(nd)=buff(m,mf,2)
          write(500+rank,*) bcl_to_bcg(mfb),t(nd),buff(m,mf,2)
!
       enddo
    enddo

    do mf=1,nbd
       call WAIT_MPI(req(mf,1))  ! waiting for all the messages to be sent
    enddo
    deallocate(buff)
!
    return
  end subroutine rfsc
end module mod_rfsc
