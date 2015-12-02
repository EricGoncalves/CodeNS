module mod_tcmd_dfpmcfg
  implicit none
contains
  subroutine tcmd_dfpmcfg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmcfg.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use kcle
    use chainecarac
    use maillage
    use boundary
    use mod_valenti
    use mod_mpi
    implicit none
    integer          ::      icmt,imot(nmx),       nm,     nmot,       no,tmp,n1,n2
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    if(knba.eq.2) knba=3
!
    if(nmot.eq.2)then
       comment=cb
       call synterr(mot,imot,2,comment)
    endif
!
    if(nmot.gt.2) then
       nm=2
       n1=0
       do no=1,maxval(bcg_to_bci)
          nm=nm+1
          if(nmot.lt.nm) then
             comment=ci
             call synterr(mot,imot,nmot,comment)
          else
             call valenti(mot,imot,nm,tmp,knba)
             do n2=1,mtb
             if(bcg_to_bci(bcl_to_bcg(n2))==tmp) then
               n1=n1+1
               nba(n1)=n2
             endif
             enddo
          endif
       enddo
       do n2=1,mtb
       if(bcg_to_bci(bcl_to_bcg(n2))==0) then
         n1=n1+1
         nba(n1)=n2
       endif
       enddo
    else
       comment=cs
       call synterr(mot,imot,nm,comment)
    endif
!
    return
  end subroutine tcmd_dfpmcfg
end module mod_tcmd_dfpmcfg
