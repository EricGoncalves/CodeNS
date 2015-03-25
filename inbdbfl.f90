module mod_inbdbfl
  implicit none
contains
  subroutine inbdbfl(ibdcfl,mfl,bceqt,mpb)
!
!***********************************************************************
!
!     ACT
!_A    Utilisation d'un etat thermodynamique pour appliquer les
!_A    conditions aux limites
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use mod_inbdbad
    implicit none
  integer          ::   ibdcfl,     img,       m,     mfl,    mflm
  integer          ::       ml,mpb(mtt),      mt,  nbidon,      nd
  integer          ::       ne,     nif,     njf,     nkf,      nr
  integer          ::       nt,      nv,     nvp,     nvt
  double precision ::             adim,bceqt(ip41,neqt),          rbidon
!
!-----------------------------------------------------------------------
!
!
    character(len=20) ::  chvar,frmt
    character(len=80) ::  titre
!
! lecture du fichier de conditions aux limites
!
    rewind(kfa)
!
    read(kfa,1000) nvt
!
    do nvp=1,nvt
!
       frmt='                    '
       read(kfa,1020) chvar,frmt
!
       if(chvar(1:2).eq.'CH'.or.chvar(1:2).eq.'ch') then
          read(kfa,1010) nt
          do m=1,nt
             read(kfa,1030) titre
          enddo
!
       elseif (chvar(1:3).eq.'INT'.or.chvar(1:3).eq.'int') then
          read(kfa,1010) ne
          if (frmt.eq.'                    ') frmt(1:6)='(10i6)'
          read(kfa,frmt) (nbidon,m=1,ne)
!
       elseif (chvar(1:4).eq.'REEL'.or.chvar(1:4).eq.'reel') then
          read(kfa,1010) nr
          if (frmt.eq.'                    ') frmt(1:8)='(6e12.5)'
          read(kfa,frmt) (rbidon,m=1,nr)
!
       else
          read(kfa,1010) nd,nif,njf,nkf,img
          write(999,1010) nd,nif,njf,nkf,img
          if (img.eq.0) img=1
          mt=nif*njf*nkf
          if (frmt.eq.'                    ') frmt(1:13)='(6e12.5)'
          if (nd.eq.ibdcfl) then
             mflm = mfl + (img-1)*mtb
             write(999,1010) mt,mfl,mflm,mtb,mpb(mflm)
!
             call inbdbad( &
                  mfl,nv,adim, &
                  chvar)
!
             read(kfa,5555) (bceqt(mpb(mflm)+m,nv),m=1,mt)
             do m=1,mt
                ml=mpb(mflm)+m
                bceqt(ml,nv)=bceqt(ml,nv)/adim
             enddo
             write(999,frmt) (bceqt(mpb(mflm)+m,nv),m=1,mt),adim
          else
             read(kfa,frmt) (rbidon,m=1,mt)
          endif
!
       endif
!
    enddo
!
1000 format (i5)
1010 format (5i6)
1020 format (2a20)
1030 format (a)
5555 format (6(e12.5))
!
    return
  end subroutine inbdbfl
end module mod_inbdbfl
