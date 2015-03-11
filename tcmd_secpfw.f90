module mod_tcmd_secpfw
implicit none
contains
      subroutine tcmd_secpfw( &
                 mot,imot,nmot, &
                 lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action secpfw.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use maillage
      use kcle
      use schemanum
use mod_synterr
use mod_valenti
use mod_vallent
implicit none
integer :: imot
integer :: nmot
integer :: lgr
integer :: lgrd
integer :: icmt
integer :: ient
integer :: im
integer :: kient
integer :: ng
integer :: ngr
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
      dimension lgr(nobj)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      if(kkvn.eq.2)     kkvn=3
      if(kncyresi.eq.2) kncyresi=3
      if(kncysave.eq.2) kncysave=3
      if(kncyexpl.eq.2) kncyexpl=3
      if(kdiscsv.eq.2)  kdiscsv=3
      do ngr=1,lg
      if(kncycle(ngr).eq.2)  kncycle(ngr)=3
      enddo
!
      lgrd=lg
      do ng=1,lgrd
      lgr(ng)=ng
      enddo
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
       nm=2
       do while(nm.lt.nmot)
        nm=nm+1
        if((imot(nm).eq.8).and.(mot(nm).eq.'verifmet')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          else
          call valenti(mot,imot,nm,kvn,kkvn)
          endif
        else if((imot(nm).eq.10).and.(mot(nm).eq.'utpostfreq')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          else
          call valenti(mot,imot,nm,ncyexpl,kncyexpl)
          endif
        else if((imot(nm).eq.8).and.(mot(nm).eq.'svfwfreq')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          else
          call valenti(mot,imot,nm,ncysave,kncysave)
          endif
        else if((imot(nm).eq.8).and.(mot(nm).eq.'svfwdisc')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ch
            call synterr(mot,imot,nmot,comment)
          else
            if(imot(nm).gt.4)then
              comment=cc
              call synterr(mot,imot,nm,comment)
            else
              discsv(1:4)='    '
              do im=1,imot(nm)
              discsv(im:im)=mot(nm)(im:im)
              enddo
              kdiscsv=2
            endif
          endif
        else if((imot(nm).eq.10).and.(mot(nm).eq.'residufreq')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          else
          call valenti(mot,imot,nm,ncyresi,kncyresi)
          endif
        else if((imot(nm).eq.11).and.(mot(nm).eq.'realisation')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ch
            call synterr(mot,imot,nmot,comment)
          else if((imot(nm).eq.3).and.(mot(nm).eq.'lgr')) then
            nm=nm+1
            if(nmot.lt.nm) then
              comment=cm
              call synterr(mot,imot,nmot,comment)
            else
              call vallent(mot,imot,nm,lgr,lgrd,lgx,klgx)
              nm=nm+1
              if((imot(nm).eq.6).and.(mot(nm).eq.'ncycle')) then
                nm=nm+1
                if(nmot.lt.nm) then
                  comment=ci
                  call synterr(mot,imot,nmot,comment)
                else
                  call valenti(mot,imot,nm,ient,kient)
                  do ng=1,lgrd
                  ncycle(lgr(ng))=ient
                  kncycle(lgr(ng))=kient
                  enddo
                endif
              else
                comment=cs
                call synterr(mot,imot,nm,comment)
              endif
            endif
          else
            comment=cs
            call synterr(mot,imot,nm,comment)
          endif
        else if(imot(nm).eq.0) then
          comment=cs
          call synterr(mot,imot,nm,comment)
        else
          comment=cb
          call synterr(mot,imot,nm,comment)
        end if
       enddo
      endif
!
      return
      end subroutine
end module
