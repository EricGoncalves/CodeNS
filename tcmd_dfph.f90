module mod_tcmd_dfph
implicit none
contains
      subroutine tcmd_dfph(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfph.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use kcle
      use proprieteflu
use mod_synterr
use mod_valreel
implicit none
integer :: imot
integer :: nmot
integer :: icmt
integer :: kval
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
      kval=0
!
      if(kgam.eq.2) kgam=3
      if(kpr.eq.2) kpr=3
      if(kprt.eq.2) kprt=3
      if(kreynz.eq.2) kreynz=3
!
      if(nmot.eq.2) then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
       nm=2
        do while(nm.lt.nmot)
        nm=nm+1
        if((imot(nm).eq.6).and.(mot(nm).eq.'gammal')) then
          nm=nm+1
            call valreel(mot,imot,nm,gam,kgam)
        else if((imot(nm).eq.5).and.(mot(nm).eq.'rgazl')) then
          nm=nm+1
            call valreel(mot,imot,nm,rd,krd)
        else if((imot(nm).eq.5).and.(mot(nm).eq.'pinfl')) then
          nm=nm+1
            call valreel(mot,imot,nm,pinfl,kpinfl)
        else if((imot(nm).eq.2).and.(mot(nm).eq.'ql')) then
          nm=nm+1
            call valreel(mot,imot,nm,ql,kql)
        else if((imot(nm).eq.8).and.(mot(nm).eq.'prandtll')) then
          nm=nm+1
            call valreel(mot,imot,nm,pr,kpr)
        else if((imot(nm).eq.8).and.(mot(nm).eq.'prandtlt')) then
          nm=nm+1
            call valreel(mot,imot,nm,prt,kprt)
        else if((imot(nm).eq.6).and.(mot(nm).eq.'reyref')) then
          nm=nm+1
            call valreel(mot,imot,nm,reynz,kreynz)
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
      end
end module
