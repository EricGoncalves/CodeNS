module mod_met_rbsr
  implicit none
contains
  subroutine met_rbsr( &
       ncbd,mnr,xnr,ynr,znr, &
       s1x,s1y,s1z,s2x,s2y,s2z, &
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
    integer          ::    l,   m,  mb,  mf, mfb
    integer          ::  mnr,  mr,  mt,ncbd,ncin
    integer          ::   nd, nid, njd,  nr
    double precision ::   cr, s1x,s1xr, s1y,s1yr
    double precision ::  s1z,s1zr, s2x,s2xr, s2y
    double precision :: s2yr, s2z,s2zr,  sr, xnr
    double precision ::  ynr, znr
!
!-----------------------------------------------------------------------
!
    dimension s1x(ip12),s1y (ip12),s1z(ip12), &
         s2x(ip12),s2y (ip12),s2z(ip12)
    dimension ncbd(ip41)
    dimension mnr(ip44),xnr(ip44),ynr(ip44),znr(ip44)
    dimension ncin(ip41)
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
       do m=1,mt
          mr=mpr(mfb)+m
          nr=mnr(mr)
          mb=mpb(mfb)+m
          nd=ncbd(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
          s1x(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s1x(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s1x(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s1x(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s1x(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s1x(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s1x(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s1x(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s1x(nr+1+nid+nid*njd)
          s2x(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s2x(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s2x(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s2x(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s2x(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s2x(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s2x(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s2x(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s2x(nr+1+nid+nid*njd)
          s2z(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s2z(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s2z(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s2z(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s2z(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s2z(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s2z(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s2z(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s2z(nr+1+nid+nid*njd)
          s1y(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s1y(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s1y(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s1y(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s1y(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s1y(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s1y(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s1y(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s1y(nr+1+nid+nid*njd)
!
          s1z(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s1z(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s1z(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s1z(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s1z(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s1z(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s1z(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s1z(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s1z(nr+1+nid+nid*njd)
!
          s2y(nd)= &
               (1-xnr(mr))*(1-ynr(mr))*(1-znr(mr))*s2y(nr              )+ &
               (1-xnr(mr))*(1-ynr(mr))*   znr(mr) *s2y(nr      +nid*njd)+ &
               (1-xnr(mr))*   ynr(mr) *(1-znr(mr))*s2y(nr  +nid        )+ &
               (1-xnr(mr))*   ynr(mr) *   znr(mr) *s2y(nr  +nid+nid*njd)+ &
               xnr(mr) *(1-ynr(mr))*(1-znr(mr))*s2y(nr+1            )+ &
               xnr(mr) *(1-ynr(mr))*   znr(mr) *s2y(nr+1    +nid*njd)+ &
               xnr(mr) *   ynr(mr) *(1-znr(mr))*s2y(nr+1+nid        )+ &
               xnr(mr) *   ynr(mr) *   znr(mr) *s2y(nr+1+nid+nid*njd)
!
!
!  attention : le sens de rotation est le sens direct.
!  ----------
!
!                   1      0      0 
!          p^-1  =  0    cos   -sin 
!                   0    sin    cos 
!
!       x                x 
!       y      = p^-1 *  y 
!       z new            z old
!
!
!        a       = p^-1 * a * p
!
          s1xr=s1x(nd)
          s1yr=cr*s1y(nd)-sr*s1z(nd)
          s1zr=cr*s1z(nd)+sr*s1y(nd)
          s2xr=cr*cr*s2x(nd)-2.*cr*sr*s2y(nd)+sr*sr*s2z(nd)
          s2yr=-cr*sr*(s2z(nd)-s2x(nd))+(2.*cr*cr-1.)*s2y(nd)
          s2zr=sr*sr*s2x(nd)+2.*cr*sr*s2y(nd)+cr*cr*s2z(nd)
!
          s1x(nd)=(s1xr+s1x(ncin(mb)))*.5
          s1y(nd)=(s1yr+s1y(ncin(mb)))*.5
          s1z(nd)=(s1zr+s1z(ncin(mb)))*.5
          s2x(nd)=(s2xr+s2x(ncin(mb)))*.5
          s2y(nd)=(s2yr+s2y(ncin(mb)))*.5
          s2z(nd)=(s2zr+s2z(ncin(mb)))*.5
!
       enddo
!
    enddo
!
    return
  end subroutine met_rbsr
end module mod_met_rbsr
