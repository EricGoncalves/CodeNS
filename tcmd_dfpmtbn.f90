module mod_tcmd_dfpmtbn
  implicit none
contains
  subroutine tcmd_dfpmtbn(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfpmtbn.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use modeleturb
    use kcle
    use mod_valenti
    implicit none
    integer          ::      icmt,imot(nmx),       nm,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
    character(len=32) ::  mot(nmx)
!
!
!$OMP SIMD
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
    if(kicytur0.eq.2) kicytur0=3
    if(kncyturb.eq.2) kncyturb=3
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
          if((imot(nm).eq.7).and.(mot(nm).eq.'icytur0')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,icytur0,kicytur0)
             endif
          else if((imot(nm).eq.7).and.(mot(nm).eq.'ncyturb')) then
             nm=nm+1
             if(nmot.lt.nm) then
                comment=ci
                call synterr(mot,imot,nmot,comment)
             else
                call valenti(mot,imot,nm,ncyturb,kncyturb)
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
  end subroutine tcmd_dfpmtbn
end module mod_tcmd_dfpmtbn
