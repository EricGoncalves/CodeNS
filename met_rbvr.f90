module mod_met_rbvr
  implicit none
contains
  subroutine met_rbvr( &
       ncbd,mnr,xnr,ynr,znr, &
       t, &
       ncin)
!
!***********************************************************************
!
!     ACT
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    implicit none
    integer          ::          l,         m,        mf,       mfb,        ml
    integer          ::  mnr(ip44),        mr,        mt,ncbd(ip41),ncin(ip41)
    integer          ::         nd,       nid,       njd,        nr
    double precision ::           cr,          sr,t(ip11,ip60),   xnr(ip44),   ynr(ip44)
    double precision ::    znr(ip44)
!
!-----------------------------------------------------------------------
!
!
    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
       sr=-srotr(mfb)
       cr=crotr(mfb)
       l=ndrr(mfb)
       nid=id2(l)-id1(l)+1
       njd=jd2(l)-jd1(l)+1
!
!!$OMP SIMD
!$OMP SIMD
       do m=1,mt
          mr=mpr(mfb)+m
          nr=mnr(mr)
          ml=mpb(mfb)+m
          nd=ncbd(ml)
!
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          t(nd,6)= &
               ((1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,6)+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,6)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,6)+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,6)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,6)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,6)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,6)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,6) &
               +t(ncin(ml),6))*.5
          t(nd,7)= &
               ((1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*t(nr              ,7)+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *t(nr      +nid*njd,7)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*t(nr  +nid        ,7)+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *t(nr  +nid+nid*njd,7)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*t(nr+1            ,7)+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *t(nr+1    +nid*njd,7)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*t(nr+1+nid        ,7)+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *t(nr+1+nid+nid*njd,7) &
               +t(ncin(ml),7))*.5
!
       enddo
!
    enddo
!
    return
  end subroutine met_rbvr
end module mod_met_rbvr
