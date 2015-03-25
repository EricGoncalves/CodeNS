module mod_tcmd_inbdc
  implicit none
contains
  subroutine tcmd_inbdc( &
       mot,imot,nmot, &
       krr,mfbea,mfbeb,kibdc,epsmsh, &
       iba,jba,kba,tvi,tvj,tvk)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action inbdc.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use mod_valenti
    use mod_valreel
    implicit none
  integer          ::       iba,     icmt,       im,imot(nmx),      jba
  integer          ::       kba,    kibdc,      krr,     kval,    mfbea
  integer          ::     mfbeb,       nm,     nmot
  double precision :: epsmsh
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
    character(len=2 ) :: tvi,tvj,tvk
!
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
    kval=0
!
    nm=3
    nm=nm+1
    if(nmot.lt.nm) then
       comment=ci
       call synterr(mot,imot,nmot,comment)
    else
       call valenti(mot,imot,nm,mfbea,kval)
    endif
!
    nm=nm+1
    if((imot(nm).eq.3).and.(mot(nm).eq.'frc')) then
       nm=nm+1
       if(nmot.lt.nm) then
          comment=ci
          call synterr(mot,imot,nmot,comment)
       else
          call valenti(mot,imot,nm,mfbeb,kval)
       endif
    else
       comment=cb
       call synterr(mot,imot,nm,comment)
    endif
!
    nm=nm+1
    if((imot(nm).eq.5).and.(mot(nm).eq.'kibdc')) then
       nm=nm+1
       if(nmot.lt.nm) then
          comment=ci
          call synterr(mot,imot,nmot,comment)
       else
          call valenti(mot,imot,nm,kibdc,kval)
       endif
    else
       comment=cb
       call synterr(mot,imot,nm,comment)
    endif
!
    if (kibdc.ne.0) then
!
!
       nm=nm+1
       if((imot(nm).eq.3).and.(mot(nm).eq.'krr')) then
          nm=nm+1
          if(nmot.lt.nm) then
             comment=ci
             call synterr(mot,imot,nmot,comment)
          else
             call valenti(mot,imot,nm,krr,kval)
          endif
       else
          comment=cb
          call synterr(mot,imot,nm,comment)
       endif
!
       if (krr.eq.0) then
!
          nm=nm+1
          if((imot(nm).eq.3).and.(mot(nm).eq.'ptc')) then
!
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,iba,kval)
             endif
!
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,jba,kval)
             endif
!
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,kba,kval)
             endif
!
          else
             comment=cb
             call synterr(mot,imot,nm,comment)
          endif
!
          nm=nm+1
          if((imot(nm).eq.3).and.(mot(nm).eq.'dir')) then
!
             nm=nm+1
             if(imot(nm).gt.2)then
                comment=cc
                call synterr(mot,imot,nm,comment)
             else
                tvi(1:2)='  '
                do im=1,imot(nm)
                   tvi(im:im)=mot(nm)(im:im)
                enddo
             endif
!
             nm=nm+1
             if(imot(nm).gt.2)then
                comment=cc
                call synterr(mot,imot,nm,comment)
             else
                tvj(1:2)='  '
                do im=1,imot(nm)
                   tvj(im:im)=mot(nm)(im:im)
                enddo
             endif
!
             nm=nm+1
             if(imot(nm).gt.2)then
                comment=cc
                call synterr(mot,imot,nm,comment)
             else
                tvk(1:2)='  '
                do im=1,imot(nm)
                   tvk(im:im)=mot(nm)(im:im)
                enddo
             endif
!
          else
             comment=cb
             call synterr(mot,imot,nm,comment)
          endif
!
       else if (krr.eq.1) then
!
          nm=nm+1
          if((imot(nm).eq.6).and.(mot(nm).eq.'epsmsh')) then
!
             nm=nm+1
             if(nmot.lt.nm) then
                comment=cr
                call synterr(mot,imot,nmot,comment)
             else
                call valreel(mot,imot,nm,epsmsh,kval)
             endif
!
          else
             comment=cb
             call synterr(mot,imot,nm,comment)
          endif
!
       endif
!
    endif
!
    return
  end subroutine tcmd_inbdc
end module mod_tcmd_inbdc
