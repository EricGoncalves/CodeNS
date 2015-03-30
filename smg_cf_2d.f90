module mod_smg_cf_2d
  implicit none
contains
  subroutine smg_cf_2d( &
       imgc,imgf, &
       vol, &
       vv,vc)
!
!-----------------------------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use schemanum
    implicit none
    integer          ::        i,     i1c,     i2c,      ic,      id
    integer          ::       if,    imgc,    imgf,  incc_i,  incc_j
    integer          ::   incc_k,  incf_i, incf_ij,incf_ijk, incf_ik
    integer          ::   incf_j, incf_jk,  incf_k,       j,     j1c
    integer          ::      j2c,      jc,      jd,      jf,       k
    integer          ::      k1c,     k2c,      kc,      kd,      kf
    integer          ::        l,     lmc,     lmf,    n0cf,    n0nc
    integer          ::      nc1,     nc2,     nc3,     nc4,     nf1
    integer          ::      nf2,     nf3,     nf4,     nic,     nif
    integer          ::     nijc,    nijf,     njc,     njf
    double precision ::           cc1,          cc2,          cc3,          cc4,          sss
    double precision :: vc(ip11,ip60),    vol(ip11),vv(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!      i,j+1            i+1,j+1
!        X----------------X
!        |        |       |
!        |    4   |   3   |
!        |        |       |
!        |--------|-------|
!        |        |       |
!        |    1   |   2   |
!        |        |       |
!        X----------------X
!       i,j             i+1,j
!
!
!




!
    do l = 1,lzx
!
       lmc  = l+(imgc-1)*lz
       lmf  = l+(imgf-1)*lz
!
       n0nc = npn(lmc)
       nic  = id2(lmc)-id1(lmc)+1
       njc  = jd2(lmc)-jd1(lmc)+1
       nijc = nic*njc
!
       n0cf = npc(lmf)
       nif  = id2(lmf)-id1(lmf)+1
       njf  = jd2(lmf)-jd1(lmf)+1
       nijf = nif*njf
!
       i1c =ii1(lmc)
       i2c =ii2(lmc)
       j1c =jj1(lmc)
       j2c =jj2(lmc)
       k1c =kk1(lmc)
       k2c =kk2(lmc)
!
       incc_i  = inc_c(1,0,0)
       incc_j  = inc_c(0,1,0)
       incc_k  = inc_c(0,0,1)
!
       incf_i   = inc_f(1,0,0)
       incf_j   = inc_f(0,1,0)
       incf_k   = inc_f(0,0,1)
       incf_ij  = inc_f(1,1,0)
       incf_ik  = inc_f(1,0,1)
       incf_jk  = inc_f(0,1,1)
       incf_ijk = inc_f(1,1,1)
!
       if (ktrans.eq.2) then
! ----- Interpolation bilineaire
          sss = 16.
          cc1 = 9./sss
          cc2 = 3./sss
          cc3 = 1./sss
          cc4 = 3./sss
       else if (ktrans.eq.1) then
! ----- Interpolation distance Centre aux 4 Noeuds
          sss  = (4*sqrt(5.)+6)/(3*sqrt(5.))
          cc1 = 1./sss
          cc2 = 1./(sqrt(5.)*sss)
          cc3 = 1./(3.*sss)
          cc4 = 1./(sqrt(5.)*sss)
       endif
!
       do kc = k1c,k2c-1
          kf =  kc
          do jc = j1c,j2c-1
             jf = 2*jc -1
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 1                              --- C
! ------------------------------------------------------------------------ C
!
!$OMP SIMD
             do ic = i1c,i2c-1
!
                nc1 = indn_c(ic,jc,kc)
                nc2 = nc1 + incc_i
                nc4 = nc1 + incc_j
                nc3 = nc4 + incc_i
!
                if  = 2*ic - 1
!
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf4 = nf1 + incf_j
                nf3 = nf4 + incf_i
! --------------------------------------------------
!          if (ktrans.eq.3) then
! --     Interplotion Bilineaire Volumique
!             cix = .5*vol(nf1)/(vol(nf1)+vol(nf2))
!             ciy = .5*vol(nf1)/(vol(nf1)+vol(nf4))
!
!             cixm1 = 1-cix
!             ciym1 = 1-ciy
!
!             cc1 = cixm1*ciym1
!             cc2 = cix  *ciym1
!             cc3 = cix  *ciy
!             cc4 = cixm1*ciy
!          endif
! --------------------------------------------------
!
                vc(nf1,1)=(cc1*vv(nc1,1)+cc2*vv(nc2,1)+cc3*vv(nc3,1)+cc4*vv(nc4,1))
                vc(nf1,2)=(cc1*vv(nc1,2)+cc2*vv(nc2,2)+cc3*vv(nc3,2)+cc4*vv(nc4,2))
                vc(nf1,3)=(cc1*vv(nc1,3)+cc2*vv(nc2,3)+cc3*vv(nc3,3)+cc4*vv(nc4,3))
                vc(nf1,4)=(cc1*vv(nc1,4)+cc2*vv(nc2,4)+cc3*vv(nc3,4)+cc4*vv(nc4,4))
                vc(nf1,5)=(cc1*vv(nc1,5)+cc2*vv(nc2,5)+cc3*vv(nc3,5)+cc4*vv(nc4,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 2                              --- C
! ------------------------------------------------------------------------ C
!
!$OMP SIMD
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
!
                if = 2*ic - 1
!
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf4 = nf1 + incf_j
                nf3 = nf4 + incf_i
!
!          if (ktrans.eq.3) then
! --     Interplotion Bilineaire Volumique
!             cix = ( vol(nf1) + .5*vol(nf2))/(vol(nf2)+vol(nf1))
!             ciy = (            .5*vol(nf2))/(vol(nf2)+vol(nf3))
!
!             cixm1 = 1-cix
!             ciym1 = 1-ciy
!
!             cc4 = cixm1*ciym1
!             cc1 = cix  *ciym1
!             cc2 = cix  *ciy
!             cc3 = cixm1*ciy
!          endif
! --------------------------------------------------
!
                vc(nf2,1)=(cc4*vv(nc1,1)+cc1*vv(nc2,1)+cc2*vv(nc3,1)+cc3*vv(nc4,1))
                vc(nf2,2)=(cc4*vv(nc1,2)+cc1*vv(nc2,2)+cc2*vv(nc3,2)+cc3*vv(nc4,2))
                vc(nf2,3)=(cc4*vv(nc1,3)+cc1*vv(nc2,3)+cc2*vv(nc3,3)+cc3*vv(nc4,3))
                vc(nf2,4)=(cc4*vv(nc1,4)+cc1*vv(nc2,4)+cc2*vv(nc3,4)+cc3*vv(nc4,4))
                vc(nf2,5)=(cc4*vv(nc1,5)+cc1*vv(nc2,5)+cc2*vv(nc3,5)+cc3*vv(nc4,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 3                              --- C
! ------------------------------------------------------------------------ C
!
!$OMP SIMD
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
!
                if = 2*ic - 1
!
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf4 = nf1 + incf_j
                nf3 = nf4 + incf_i
!
!          if (ktrans.eq.3) then
! --     Interplotion Bilineaire Volumique
!             cix = ( vol(nf4) + .5*vol(nf3))/(vol(nf3)+vol(nf4))
!             ciy = ( vol(nf2) + .5*vol(nf3))/(vol(nf3)+vol(nf2))
!
!             cixm1 = 1-cix
!             ciym1 = 1-ciy
!
!             cc3 = cixm1*ciym1
!             cc4 = cix  *ciym1
!             cc1 = cix  *ciy
!             cc2 = cixm1*ciy
!          endif
! --------------------------------------------------
!
                vc(nf3,1)=(cc3*vv(nc1,1)+cc4*vv(nc2,1)+cc1*vv(nc3,1)+cc2*vv(nc4,1))
                vc(nf3,2)=(cc3*vv(nc1,2)+cc4*vv(nc2,2)+cc1*vv(nc3,2)+cc2*vv(nc4,2))
                vc(nf3,3)=(cc3*vv(nc1,3)+cc4*vv(nc2,3)+cc1*vv(nc3,3)+cc2*vv(nc4,3))
                vc(nf3,4)=(cc3*vv(nc1,4)+cc4*vv(nc2,4)+cc1*vv(nc3,4)+cc2*vv(nc4,4))
                vc(nf3,5)=(cc3*vv(nc1,5)+cc4*vv(nc2,5)+cc1*vv(nc3,5)+cc2*vv(nc4,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 4                              --- C
! ------------------------------------------------------------------------ C
!
!$OMP SIMD
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
!
                if = 2*ic - 1
!
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf4 = nf1 + incf_j
                nf3 = nf4 + incf_i
!
! --------------------------------------------------
!          if (ktrans.eq.3) then
! --     Interplotion Bilineaire Volumique
!             cix = (            .5*vol(nf4))/(vol(nf4)+vol(nf3))
!             ciy = ( vol(nf1) + .5*vol(nf4))/(vol(nf4)+vol(nf1))
!
!             cixm1 = 1-cix
!             ciym1 = 1-ciy
!
!             cc2 = cixm1*ciym1
!             cc3 = cix  *ciym1
!             cc4 = cix  *ciy
!             cc1 = cixm1*ciy
!          endif
! --------------------------------------------------
!
                vc(nf4,1)=(cc2*vv(nc1,1)+cc3*vv(nc2,1)+cc4*vv(nc3,1)+cc1*vv(nc4,1))
                vc(nf4,2)=(cc2*vv(nc1,2)+cc3*vv(nc2,2)+cc4*vv(nc3,2)+cc1*vv(nc4,2))
                vc(nf4,3)=(cc2*vv(nc1,3)+cc3*vv(nc2,3)+cc4*vv(nc3,3)+cc1*vv(nc4,3))
                vc(nf4,4)=(cc2*vv(nc1,4)+cc3*vv(nc2,4)+cc4*vv(nc3,4)+cc1*vv(nc4,4))
                vc(nf4,5)=(cc2*vv(nc1,5)+cc3*vv(nc2,5)+cc4*vv(nc3,5)+cc1*vv(nc4,5))
             enddo
!
! ------------------------------------------------------------------------ C
!
          enddo
       enddo
!
    enddo
!
    return
  contains
    function    indn_c(i,j,k)
      implicit none
      integer          ::      i,indn_c,     j,     k
      indn_c=n0nc+1+(i-id1(lmc))+(j-jd1(lmc))*nic+(k-kd1(lmc))*nijc
    end function indn_c
    function    indc_f(i,j,k)
      implicit none
      integer          ::      i,indc_f,     j,     k
      indc_f=n0cf+1+(i-id1(lmf))+(j-jd1(lmf))*nif+(k-kd1(lmf))*nijf
    end function indc_f
    function    inc_c(id,jd,kd)
      implicit none
      integer          ::    id,inc_c,   jd,   kd
      inc_c=id+jd*nic+kd*nijc
    end function inc_c
    function    inc_f(id,jd,kd)
      implicit none
      integer          ::    id,inc_f,   jd,   kd
      inc_f=id+jd*nif+kd*nijf
    end function inc_f
  end subroutine smg_cf_2d
end module mod_smg_cf_2d
