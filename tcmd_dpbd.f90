module mod_tcmd_dpbd
implicit none
contains
      subroutine tcmd_dpbd( &
                 mot,imot,nmot, &
                 lmfb,lmfbd, &
                 lgr,lgrd, &
                 typdat)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dpbd.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use maillage 
      use kcle
use mod_synterr
use mod_vallent
implicit none
integer :: imot
integer :: nmot
integer :: lmfb
integer :: lmfbd
integer :: lgr
integer :: lgrd
integer :: icmt
integer :: im
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      character(len=32) ::  typdat
      dimension imot(nmx)
      dimension lmfb(mtb)
      dimension lgr(nobj)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      lgrd=1
      lgr(1)=1
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      nm=2
      nm=nm+1
      if(nmot.lt.nm) then
        comment=ch
        call synterr(mot,imot,nmot,comment)
      else
        typdat(1:32)='                                '
        do im=1,imot(nm)
        typdat(im:im)=mot(nm)(im:im)
        enddo
        nm=nm+1
        if((imot(nm).eq.4).and.(mot(nm).eq.'lmfb')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=cm
            call synterr(mot,imot,nmot,comment)
          else
            call vallent(mot,imot,nm,lmfb,lmfbd,mtbx,kmtbx)
            nm=nm+1
            if((imot(nm).eq.3).and.(mot(nm).eq.'lgr')) then
              nm=nm+1
              if(nmot.lt.nm) then
                comment=cm
                call synterr(mot,imot,nmot,comment)
              else
                call vallent(mot,imot,nm,lgr,lgrd,lgx,klgx)
              endif
            else
              comment=cb
              call synterr(mot,imot,nm,comment)
            endif
          endif
        else
          comment=cb
          call synterr(mot,imot,nm,comment)
        endif
      endif
!
      return
      end
end module
