module mod_tcmd_dfgm
implicit none
contains
      subroutine tcmd_dfgm(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfgm.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use kcle
      use definition
use mod_synterr
use mod_valreel
implicit none
integer :: imot
integer :: nmot
integer :: icmt
integer :: im
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      if(kconfig.eq.2) kconfig=3
      if(kperio.eq.2) kperio=3
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
        if((imot(nm).eq.6).and.(mot(nm).eq.'config')) then
          nm=nm+1
          if(imot(nm).gt.4)then
            comment=cc
            call synterr(mot,imot,nm,comment)
          else
            config(1:4)='    '
            do im=1,imot(nm)
            config(im:im)=mot(nm)(im:im)
            enddo
            kconfig=2
          endif
        else if((imot(nm).eq.5).and.(mot(nm).eq.'perio')) then
          nm=nm+1
            call valreel(mot,imot,nm,perio,kperio)
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
