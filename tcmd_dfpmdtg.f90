module mod_tcmd_dfpmdtg
  implicit none
contains
  subroutine tcmd_dfpmdtg(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmdtg.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use kcle
    use schemanum
    use mod_valenti
    use mod_valreel
    implicit none
    integer          ::      icmt,imot(nmx),       nm,     nmot
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
    if(kkdtl.eq.2) kkdtl=3
    if(kicychr0.eq.2) kicychr0=3
    if(kncychro.eq.2) kncychro=3
    if(kdt1min.eq.2) kdt1min=3
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
          if((imot(nm).eq.4).and.(mot(nm).eq.'kdtl')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,kdtl,kkdtl)
             endif
          else if((imot(nm).eq.7).and.(mot(nm).eq.'icychr0')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,icychr0,kicychr0)
             endif
          else if((imot(nm).eq.7).and.(mot(nm).eq.'ncychro')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,ncychro,kncychro)
             endif
          else if((imot(nm).eq.6).and.(mot(nm).eq.'dt1min')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=cr
                call synterr(mot,imot,nmot,comment)
             else
                call valreel(mot,imot,nm,dt1min,kdt1min)
             endif
          else if(imot(nm).eq.0) then
             comment=cs
             call synterr(mot,imot,nm,comment)
          else
             comment=cb
             call synterr(mot,imot,nm,comment)
          endif
       enddo
    endif
!
    return
  end subroutine tcmd_dfpmdtg
end module mod_tcmd_dfpmdtg
