module mod_tcmd_dfpmimd
  implicit none
contains
  subroutine tcmd_dfpmimd( &
       mot,imot,nmot, &
       ldom,ldomd, &
       lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a l'action dfpmimd.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use schemanum
    use kcle
    use chainecarac
    use maillage
    use mod_valenti
    use mod_vallent
    implicit none
  integer          ::       icmt,      ient,       img, imot(nmx),     kient
  integer          ::          l,ldom(nobj),     ldomd, lgr(nobj),      lgrd
  integer          ::         lm,        ng,        nl,        nm,      nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    do l=1,lt
       if(kkmf(l).eq.2) kkmf(l)=3
    enddo
!
    if(nmot.eq.2)then
       comment=cb
       call synterr(mot,imot,2,comment)
    endif
!
    nm=2
!
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
!
                do while(nm.lt.nmot)
                   nm=nm+1
                   if((imot(nm).eq.3).and.(mot(nm).eq.'kmf')) then
                      nm=nm+1
                      if(nmot.lt.nm) then
                         comment=ci
                         call synterr(mot,imot,nmot,comment)
                      else
                         call valenti(mot,imot,nm,ient,kient)
                         do nl=1,ldomd
                            l=ldom(nl)
                            do ng=1,lgrd
                               img=lgr(ng)
                               lm=l+(img-1)*lz
                               kmf(lm)=ient
                               kkmf(lm)=kient
                            enddo
                         enddo
                      endif
                   elseif((imot(nm).eq.4).and.(mot(nm).eq.'lmax')) then
                      nm=nm+1
                      if(nmot.lt.nm) then
                         comment=ci
                         call synterr(mot,imot,nmot,comment)
                      else
                         call valenti(mot,imot,nm,ient,kient)
                         do nl=1,ldomd
                            l=ldom(nl)
                            do ng=1,lgrd
                               img=lgr(ng)
                               lm=l+(img-1)*lz
                               lmax(lm)=ient
                               klmax(lm)=kient
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
                enddo
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
  end subroutine tcmd_dfpmimd
end module mod_tcmd_dfpmimd
