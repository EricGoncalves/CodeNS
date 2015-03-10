module mod_tcmd_dfpmdtd
implicit none
contains
      subroutine tcmd_dfpmdtd( &
                 mot,imot,nmot, &
                 ldom,ldomd, &
                 lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmdtd.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use kcle
      use maillage
      use schemanum
use mod_synterr
use mod_valreel
use mod_vallent
implicit none
integer :: imot
integer :: nmot
integer :: ldom
integer :: ldomd
integer :: lgr
integer :: lgrd
double precision :: et
integer :: icmt
integer :: img
integer :: ket
integer :: l
integer :: lm
integer :: ng
integer :: nl
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
!
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension lgr(nobj)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
      do l=1,lt
      if(keta(l).eq.2) keta(l)=3
      enddo
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      nm=2
      nm=nm+1
      if((imot(nm).eq.4).and.(mot(nm).eq.'ldom')) then
        nm=nm+1
        if(nmot.lt.nm) then
          comment=cm
          call synterr(mot,imot,nmot,comment)
        else
          call vallent(mot,imot,nm,ldom,ldomd,lzx,klzx)
          nm=nm+1
          if((imot(nm).eq.3).and.(mot(nm).eq.'lgr')) then
            nm=nm+1
            if(nmot.lt.nm) then
              comment=cm
              call synterr(mot,imot,nmot,comment)
            else
              call vallent(mot,imot,nm,lgr,lgrd,lgx,klgx)
              nm=nm+1
              if((imot(nm).eq.3).and.(mot(nm).eq.'eta')) then
                nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
                  call valreel(mot,imot,nm,et,ket)
                  do nl=1,ldomd
                  l=ldom(nl)
                  do ng=1,lgrd
                  img=lgr(ng)
                  lm=l+(img-1)*lz
                  eta(lm)=et
                  keta(lm)=ket
                  enddo
                  enddo
                endif
                else if(imot(nm).eq.0) then
                  comment=cs
                  call synterr(mot,imot,nm,comment)
              else
                comment=cb
                call synterr(mot,imot,nm,comment)
              end if
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
!
      return
      end
end module
