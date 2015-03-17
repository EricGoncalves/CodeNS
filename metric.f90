module mod_metric
  implicit none
contains
  subroutine metric( &
       x,y,z,r,exs1,exs2, &
       sn,vol, &
       ncbd,mnc, &
       mnr,xnr,ynr,znr, &
       tn1,tn2,tn3)
!
!***********************************************************************
!
!     ACT
!_A    Appel d'un calcul de metrique.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    exs1       : arg real             ; premier coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    exs2       : arg real             ; deuxieme coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    lzx        : com int              ; nbr total de domaines
!
!_O    r          : arg real(ip11      ) ; distance a l'axe
!_O    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use schemanum
    use chainecarac 
    use mod_metrics
    use mod_rfsc
    use mod_rfsr
    implicit none
    integer          ::  img,   l,  lm, mfc, mfr
    integer          ::  mnc, mnr,ncbd
    double precision :: exs1,exs2,   r,  sn, tn1
    double precision ::  tn2, tn3, vol,   x, xnr
    double precision ::    y, ynr,   z, znr
!
!-----------------------------------------------------------------------
!
    dimension x(ip21),y(ip21),z(ip21),r(ip11),vol(ip11)
    dimension sn(ip31*ndir)
    dimension ncbd(ip41),mnc(ip43)
    dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
    dimension tn1(ip00),tn2(ip00),tn3(ip00)
!
    do l=1,lzx
!
       do img=1,lgx
          lm=l+(img-1)*lz
          if ((img.eq.1).or.(kcg.eq.1)) then
             call metrics( &
                  lm,x,y,z,r,exs1,exs2, &
                  sn,vol, &
                  tn1,tn2,tn3)
          else
             stop 'img ne 1 et kcg ne 1 !!! '
          endif
       enddo
!
    enddo
!
    if (mtcx.gt.0) then
       do img=1,lgx
          do mfc=1,mtcx
             lbd(mfc)=nfbc(mfc)+(img-1)*mtb
          enddo
          nbd=mtcx
          call rfsc(vol,ncbd,mnc)
!
       enddo
    endif
!
    if (mtrx.gt.0) then
       do img=1,lgx
          do mfr=1,mtrx
             lbd(mfr)=nfbr(mfr)+(img-1)*mtb
          enddo
          nbd=mtrx
          call rfsr(vol,ncbd,mnr,xnr,ynr,znr)
!
       enddo
    endif
!
    return
  end subroutine metric
end module mod_metric
