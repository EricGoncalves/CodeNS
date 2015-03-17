module mod_tcmd_dfnm
  implicit none
contains
  subroutine tcmd_dfnm(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a
!_A    l'action dfnm.
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use chainecarac
    use kcle
    use schemanum
    use maillage
    use proprieteflu
    use mod_valenti
    use mod_valreel
    implicit none
    integer          :: icmt,imot,  nm,nmot
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
    if(knfi.eq.2) knfi=3
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
          if((imot(nm).eq.3).and.(mot(nm).eq.'nfi')) then
             nm=nm+1
             call valenti(mot,imot,nm,nfi,knfi)
          else if((imot(nm).eq.3).and.(mot(nm).eq.'fmg')) then
             nm=nm+1
             call valenti(mot,imot,nm,kfmg,kkfmg)
             nm=nm+1
             call valenti(mot,imot,nm,lgx,klgx)
             nm=nm+1
             call valenti(mot,imot,nm,kcg,kkcg)
          else if((imot(nm).eq.6).and.(mot(nm).eq.'schema')) then
             nm=nm+1
             call valenti(mot,imot,nm,ischema,kischema)
          else if((imot(nm).eq.5).and.(mot(nm).eq.'muscl')) then
             nm=nm+1
             call valenti(mot,imot,nm,muscl,kmuscl)
             nm=nm+1
             call valreel(mot,imot,nm,xk,kxk)
             nm=nm+1
             call valenti(mot,imot,nm,ilim,kilim)
          else if((imot(nm).eq.4).and.(mot(nm).eq.'prcd')) then
             nm=nm+1
             call valenti(mot,imot,nm,kprec,kkprec)
             nm=nm+1
             call valreel(mot,imot,nm,cte,kcte)
             nm=nm+1
             call valenti(mot,imot,nm,kvisq,kkvisq)
          else if((imot(nm).eq.4).and.(mot(nm).eq.'dual')) then
             nm=nm+1
             call valenti(mot,imot,nm,kdualns,kkdualns)
             nm=nm+1
             call valreel(mot,imot,nm,tol,ktol)
             nm=nm+1
             call valenti(mot,imot,nm,niter,kniter)
             nm=nm+1
             call valreel(mot,imot,nm,tolke,ktolke)
             nm=nm+1
             call valenti(mot,imot,nm,nitur,knitur)
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
  end subroutine tcmd_dfnm
end module mod_tcmd_dfnm
