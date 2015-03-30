module mod_tcmd_inbdb
  implicit none
contains
  subroutine tcmd_inbdb( &
       mot,imot,nmot, &
       lmfb,lmfbd,clmf,kibdb, &
       ibdcst,ibdcfl,ibddim,nvbc,vbc)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action inbdb.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use maillage
    use kcle
    use mod_valenti
    use mod_valreel
    use mod_vallent
    implicit none
    integer          ::    ibdcfl,   ibdcst,   ibddim,     icmt,       im
    integer          :: imot(nmx),    kibdb,     kval,lmfb(mtb),    lmfbd
    integer          ::        nm,     nmot,      nmr,     nvbc
    double precision :: vbc(ista*lsta)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
    character(len=4 ) :: clmf
!
!$OMP SIMD
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    kval=0
    ibdcst=0
    ibdcfl=0
    ibddim=0
    nvbc=0
!
    nm=3
    nm=nm+1
    if(nmot.lt.nm) then
       comment=cm
       call synterr(mot,imot,nmot,comment)
    else
       call vallent(mot,imot,nm,lmfb,lmfbd,mtbx,kmtbx)
    endif
!
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ch
       call synterr(mot,imot,nmot,comment)
    else
       if(imot(nm).gt.4)then
          comment=cc
          call synterr(mot,imot,nm,comment)
       else
          clmf(1:4)='    '
!$OMP SIMD
          do im=1,imot(nm)
             clmf(im:im)=mot(nm)(im:im)
          enddo
       endif
    endif
!
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,kibdb,kval)
    endif
!
    nm=nm+1
    if(nmot.lt.nm) return
!
    if((imot(nm).eq.4).and.(mot(nm).eq.'file')) then
       nm=nm+1
       call valenti(mot,imot,nm,ibdcfl,kval)
    else
       nm=nm-1
    endif
!
    nm=nm+1
    if(nmot.lt.nm) return
!
    if((imot(nm).eq.5).and.(mot(nm).eq.'state')) then
       ibddim=1
       nm=nm+1
       call valenti(mot,imot,nm,ibdcst,kval)
!
    else if((imot(nm).eq.7).and.(mot(nm).eq.'val/usi')) then
       ibddim=1
       nm=nm+1
       nvbc=0
       do nmr=nm,nmot
          nvbc=nvbc+1
          call valreel(mot,imot,nmr,vbc(nvbc),kval)
       enddo
       nm=nmot
!
    else if((imot(nm).eq.7).and.(mot(nm).eq.'val/ref')) then
       ibddim=0
       nm=nm+1
       nvbc=0
       do nmr=nm,nmot
          nvbc=nvbc+1
          call valreel(mot,imot,nmr,vbc(nvbc),kval)
       enddo
       nm=nmot
!
    else
       comment=cb
       call synterr(mot,imot,nm,comment)
    endif
!
    return
  end subroutine tcmd_inbdb
end module mod_tcmd_inbdb
