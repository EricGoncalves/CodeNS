module mod_smg_cf_3d
  implicit none
contains
  subroutine smg_cf_3d( &
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
    integer          ::      nc1,     nc2,     nc3,     nc4,     nc5
    integer          ::      nc6,     nc7,     nc8,     nf1,     nf2
    integer          ::      nf3,     nf4,     nf5,     nf6,     nf7
    integer          ::      nf8,     nic,     nif,    nijc,    nijf
    integer          ::      njc,     njf
    double precision ::           cc1,          cc2,          cc3,          cc4,          cc5
    double precision ::           cc6,          cc7,          cc8,          cix,        cixm1
    double precision ::           ciy,        ciym1,          ciz,        cizm1,          sss
    double precision :: vc(ip11,ip60),    vol(ip11),vv(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!     i,j+1,k          i+1,j+1,k     i,j+1,k+1        i+1,j+1,k+1
!        X----------------X              X----------------X
!        |        |       |              |        |       |
!        |    4   |   3   |              |   8    |   7   |
!        |        |       |              |        |       |
!        |--------|-------|              |--------|-------|
!        |        |       |              |        |       |
!        |    1   |   2   |              |   5    |   6   |
!        |        |       |              |        |       |
!        X----------------X              X----------------X
!       i,j,k           i+1,j,k        i,j,k+1        i+1,j,k+1
!
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
          sss =  64.
          cc1 = 27./sss
          cc2 =  9./sss
          cc3 =  3./sss
          cc4 =  9./sss
          cc5 =  9./sss
          cc6 =  3./sss
          cc7 =  1./sss
          cc8 =  3./sss
       elseif (ktrans.eq.1) then
          sss =  4./(3.*sqrt(3.)) + 4./sqrt(11.) + 2./sqrt(19.)
          cc1 =  1./ (sqrt(3.)*sss)
          cc2 =  1./ (sqrt(11.)*sss)
          cc3 =  1./ (sqrt(19.)*sss)
          cc4 =  1./ (sqrt(11.)*sss)
          cc5 =  1./ (sqrt(11.)*sss)
          cc6 =  1./ (sqrt(19.)*sss)
          cc7 =  1./ (3*sqrt(3.)*sss)
          cc8 =  1./ (sqrt(11.)*sss)
       endif
!
       do kc = k1c,k2c-1
          kf = 2*kc -1
          do jc = j1c,j2c-1
             jf = 2*jc -1
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 1                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = .5*vol(nf1)/(vol(nf1)+vol(nf2))
                   ciy = .5*vol(nf1)/(vol(nf1)+vol(nf4))
                   ciz = .5*vol(nf1)/(vol(nf1)+vol(nf5))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc1 = cixm1*ciym1*cizm1
                   cc2 = cix  *ciym1*cizm1
                   cc3 = cix  *ciy  *cizm1
                   cc4 = cixm1*ciy  *cizm1
                   cc5 = cixm1*ciym1*ciz
                   cc6 = cix  *ciym1*ciz
                   cc7 = cix  *ciy  *ciz
                   cc8 = cixm1*ciy  *ciz
                endif
!
                vc(nf1,1)=( cc1*vv(nc1,1) + cc2*vv(nc2,1) + cc3*vv(nc3,1) + &
                     cc4*vv(nc4,1) + cc5*vv(nc5,1) + cc6*vv(nc6,1) + &
                     cc7*vv(nc7,1) + cc8*vv(nc8,1))
                vc(nf1,2)=( cc1*vv(nc1,2) + cc2*vv(nc2,2) + cc3*vv(nc3,2) + &
                     cc4*vv(nc4,2) + cc5*vv(nc5,2) + cc6*vv(nc6,2) + &
                     cc7*vv(nc7,2) + cc8*vv(nc8,2))
                vc(nf1,3)=( cc1*vv(nc1,3) + cc2*vv(nc2,3) + cc3*vv(nc3,3) + &
                     cc4*vv(nc4,3) + cc5*vv(nc5,3) + cc6*vv(nc6,3) + &
                     cc7*vv(nc7,3) + cc8*vv(nc8,3))
                vc(nf1,4)=( cc1*vv(nc1,4) + cc2*vv(nc2,4) + cc3*vv(nc3,4) + &
                     cc4*vv(nc4,4) + cc5*vv(nc5,4) + cc6*vv(nc6,4) + &
                     cc7*vv(nc7,4) + cc8*vv(nc8,4))
                vc(nf1,5)=( cc1*vv(nc1,5) + cc2*vv(nc2,5) + cc3*vv(nc3,5) + &
                     cc4*vv(nc4,5) + cc5*vv(nc5,5) + cc6*vv(nc6,5) + &
                     cc7*vv(nc7,5) + cc8*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 2                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if(ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = ( vol(nf1) + .5*vol(nf2) )/(vol(nf2)+vol(nf1))
                   ciy = (            .5*vol(nf2) )/(vol(nf2)+vol(nf3))
                   ciz = (            .5*vol(nf2) )/(vol(nf2)+vol(nf6))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc2 = cixm1*ciym1*cizm1
                   cc1 = cix  *ciym1*cizm1
                   cc4 = cix  *ciy  *cizm1
                   cc3 = cixm1*ciy  *cizm1
                   cc6 = cixm1*ciym1*ciz
                   cc5 = cix  *ciym1*ciz
                   cc8 = cix  *ciy  *ciz
                   cc7 = cixm1*ciy  *ciz
                endif
!
                vc(nf2,1)=( cc2*vv(nc1,1) + cc1*vv(nc2,1) + cc4*vv(nc3,1) + &
                     cc3*vv(nc4,1) + cc6*vv(nc5,1) + cc5*vv(nc6,1) + &
                     cc8*vv(nc7,1) + cc7*vv(nc8,1))
                vc(nf2,2)=( cc2*vv(nc1,2) + cc1*vv(nc2,2) + cc4*vv(nc3,2) + &
                     cc3*vv(nc4,2) + cc6*vv(nc5,2) + cc5*vv(nc6,2) + &
                     cc8*vv(nc7,2) + cc7*vv(nc8,2))
                vc(nf2,3)=( cc2*vv(nc1,3) + cc1*vv(nc2,3) + cc4*vv(nc3,3) + &
                     cc3*vv(nc4,3) + cc6*vv(nc5,3) + cc5*vv(nc6,3) + &
                     cc8*vv(nc7,3) + cc7*vv(nc8,3))
                vc(nf2,4)=( cc2*vv(nc1,4) + cc1*vv(nc2,4) + cc4*vv(nc3,4) + &
                     cc3*vv(nc4,4) + cc6*vv(nc5,4) + cc5*vv(nc6,4) + &
                     cc8*vv(nc7,4) + cc7*vv(nc8,4))
                vc(nf2,5)=( cc2*vv(nc1,5) + cc1*vv(nc2,5) + cc4*vv(nc3,5) + &
                     cc3*vv(nc4,5) + cc6*vv(nc5,5) + cc5*vv(nc6,5) + &
                     cc8*vv(nc7,5) + cc7*vv(nc8,5))
!
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 3                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = ( vol(nf4) + .5*vol(nf3) )/(vol(nf3)+vol(nf4))
                   ciy = ( vol(nf2) + .5*vol(nf3) )/(vol(nf3)+vol(nf2))
                   ciz = (            .5*vol(nf3) )/(vol(nf3)+vol(nf7))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc3 = cixm1*ciym1*cizm1
                   cc4 = cix  *ciym1*cizm1
                   cc1 = cix  *ciy  *cizm1
                   cc2 = cixm1*ciy  *cizm1
                   cc7 = cixm1*ciym1*ciz
                   cc8 = cix  *ciym1*ciz
                   cc5 = cix  *ciy  *ciz
                   cc6 = cixm1*ciy  *ciz
                endif
!
                vc(nf3,1)=( cc3*vv(nc1,1) + cc4*vv(nc2,1) + cc1*vv(nc3,1) + &
                     cc2*vv(nc4,1) + cc7*vv(nc5,1) + cc8*vv(nc6,1) + &
                     cc5*vv(nc7,1) + cc6*vv(nc8,1))
                vc(nf3,2)=( cc3*vv(nc1,2) + cc4*vv(nc2,2) + cc1*vv(nc3,2) + &
                     cc2*vv(nc4,2) + cc7*vv(nc5,2) + cc8*vv(nc6,2) + &
                     cc5*vv(nc7,2) + cc6*vv(nc8,2))
                vc(nf3,3)=( cc3*vv(nc1,3) + cc4*vv(nc2,3) + cc1*vv(nc3,3) + &
                     cc2*vv(nc4,3) + cc7*vv(nc5,3) + cc8*vv(nc6,3) + &
                     cc5*vv(nc7,3) + cc6*vv(nc8,3))
                vc(nf3,4)=( cc3*vv(nc1,4) + cc4*vv(nc2,4) + cc1*vv(nc3,4) + &
                     cc2*vv(nc4,4) + cc7*vv(nc5,4) + cc8*vv(nc6,4) + &
                     cc5*vv(nc7,4) + cc6*vv(nc8,4))
                vc(nf3,5)=( cc3*vv(nc1,5) + cc4*vv(nc2,5) + cc1*vv(nc3,5) + &
                     cc2*vv(nc4,5) + cc7*vv(nc5,5) + cc8*vv(nc6,5) + &
                     cc5*vv(nc7,5) + cc6*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 4                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = (            .5*vol(nf4) )/(vol(nf4)+vol(nf3))
                   ciy = ( vol(nf1) + .5*vol(nf4) )/(vol(nf4)+vol(nf1))
                   ciz = (            .5*vol(nf4) )/(vol(nf4)+vol(nf8))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc4 = cixm1*ciym1*cizm1
                   cc3 = cix  *ciym1*cizm1
                   cc2 = cix  *ciy  *cizm1
                   cc1 = cixm1*ciy  *cizm1
                   cc8 = cixm1*ciym1*ciz
                   cc7 = cix  *ciym1*ciz
                   cc6 = cix  *ciy  *ciz
                   cc5 = cixm1*ciy  *ciz
                endif
!
                vc(nf4,1)=( cc4*vv(nc1,1) + cc3*vv(nc2,1) + cc2*vv(nc3,1) + &
                     cc1*vv(nc4,1) + cc8*vv(nc5,1) + cc7*vv(nc6,1) + &
                     cc6*vv(nc7,1) + cc5*vv(nc8,1))
                vc(nf4,2)=( cc4*vv(nc1,2) + cc3*vv(nc2,2) + cc2*vv(nc3,2) + &
                     cc1*vv(nc4,2) + cc8*vv(nc5,2) + cc7*vv(nc6,2) + &
                     cc6*vv(nc7,2) + cc5*vv(nc8,2))
                vc(nf4,3)=( cc4*vv(nc1,3) + cc3*vv(nc2,3) + cc2*vv(nc3,3) + &
                     cc1*vv(nc4,3) + cc8*vv(nc5,3) + cc7*vv(nc6,3) + &
                     cc6*vv(nc7,3) + cc5*vv(nc8,3))
                vc(nf4,4)=( cc4*vv(nc1,4) + cc3*vv(nc2,4) + cc2*vv(nc3,4) + &
                     cc1*vv(nc4,4) + cc8*vv(nc5,4) + cc7*vv(nc6,4) + &
                     cc6*vv(nc7,4) + cc5*vv(nc8,4))
                vc(nf4,5)=( cc4*vv(nc1,5) + cc3*vv(nc2,5) + cc2*vv(nc3,5) + &
                     cc1*vv(nc4,5) + cc8*vv(nc5,5) + cc7*vv(nc6,5) + &
                     cc6*vv(nc7,5) + cc5*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 5                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
!--Interplotion Trilineaire Volumique
                   cix = (            .5*vol(nf5) )/(vol(nf5)+vol(nf6))
                   ciy = (            .5*vol(nf5) )/(vol(nf5)+vol(nf8))
                   ciz = ( vol(nf1) + .5*vol(nf5) )/(vol(nf5)+vol(nf1))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc5 = cixm1*ciym1*cizm1
                   cc6 = cix  *ciym1*cizm1
                   cc7 = cix  *ciy  *cizm1
                   cc8 = cixm1*ciy  *cizm1
                   cc1 = cixm1*ciym1*ciz
                   cc2 = cix  *ciym1*ciz
                   cc3 = cix  *ciy  *ciz
                   cc4 = cixm1*ciy  *ciz
                endif
!
                vc(nf5,1)=( cc5*vv(nc1,1) + cc6*vv(nc2,1) + cc7*vv(nc3,1) + &
                     cc8*vv(nc4,1) + cc1*vv(nc5,1) + cc2*vv(nc6,1) + &
                     cc3*vv(nc7,1) + cc4*vv(nc8,1))
                vc(nf5,2)=( cc5*vv(nc1,2) + cc6*vv(nc2,2) + cc7*vv(nc3,2) + &
                     cc8*vv(nc4,2) + cc1*vv(nc5,2) + cc2*vv(nc6,2) + &
                     cc3*vv(nc7,2) + cc4*vv(nc8,2))
                vc(nf5,3)=( cc5*vv(nc1,3) + cc6*vv(nc2,3) + cc7*vv(nc3,3) + &
                     cc8*vv(nc4,3) + cc1*vv(nc5,3) + cc2*vv(nc6,3) + &
                     cc3*vv(nc7,3) + cc4*vv(nc8,3))
                vc(nf5,4)=( cc5*vv(nc1,4) + cc6*vv(nc2,4) + cc7*vv(nc3,4) + &
                     cc8*vv(nc4,4) + cc1*vv(nc5,4) + cc2*vv(nc6,4) + &
                     cc3*vv(nc7,4) + cc4*vv(nc8,4))
                vc(nf5,5)=( cc5*vv(nc1,5) + cc6*vv(nc2,5) + cc7*vv(nc3,5) + &
                     cc8*vv(nc4,5) + cc1*vv(nc5,5) + cc2*vv(nc6,5) + &
                     cc3*vv(nc7,5) + cc4*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 6                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interpolation Trilineaire Volumique
                   cix = ( vol(nf5) + .5*vol(nf6) )/(vol(nf6)+vol(nf5))
                   ciy = (            .5*vol(nf6) )/(vol(nf6)+vol(nf7))
                   ciz = ( vol(nf2) + .5*vol(nf6) )/(vol(nf6)+vol(nf2))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc6 = cixm1*ciym1*cizm1
                   cc5 = cix  *ciym1*cizm1
                   cc8 = cix  *ciy  *cizm1
                   cc7 = cixm1*ciy  *cizm1
                   cc2 = cixm1*ciym1*ciz
                   cc1 = cix  *ciym1*ciz
                   cc4 = cix  *ciy  *ciz
                   cc3 = cixm1*ciy  *ciz
                endif
!
                vc(nf6,1)=( cc6*vv(nc1,1) + cc5*vv(nc2,1) + cc8*vv(nc3,1) + &
                     cc7*vv(nc4,1) + cc2*vv(nc5,1) + cc1*vv(nc6,1) + &
                     cc4*vv(nc7,1) + cc3*vv(nc8,1))
                vc(nf6,2)=( cc6*vv(nc1,2) + cc5*vv(nc2,2) + cc8*vv(nc3,2) + &
                     cc7*vv(nc4,2) + cc2*vv(nc5,2) + cc1*vv(nc6,2) + &
                     cc4*vv(nc7,2) + cc3*vv(nc8,2))
                vc(nf6,3)=( cc6*vv(nc1,3) + cc5*vv(nc2,3) + cc8*vv(nc3,3) + &
                     cc7*vv(nc4,3) + cc2*vv(nc5,3) + cc1*vv(nc6,3) + &
                     cc4*vv(nc7,3) + cc3*vv(nc8,3))
                vc(nf6,4)=( cc6*vv(nc1,4) + cc5*vv(nc2,4) + cc8*vv(nc3,4) + &
                     cc7*vv(nc4,4) + cc2*vv(nc5,4) + cc1*vv(nc6,4) + &
                     cc4*vv(nc7,4) + cc3*vv(nc8,4))
                vc(nf6,5)=( cc6*vv(nc1,5) + cc5*vv(nc2,5) + cc8*vv(nc3,5) + &
                     cc7*vv(nc4,5) + cc2*vv(nc5,5) + cc1*vv(nc6,5) + &
                     cc4*vv(nc7,5) + cc3*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 7                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = ( vol(nf8) + .5*vol(nf7) )/(vol(nf7)+vol(nf8))
                   ciy = ( vol(nf6) + .5*vol(nf7) )/(vol(nf7)+vol(nf6))
                   ciz = ( vol(nf3) + .5*vol(nf7) )/(vol(nf7)+vol(nf3))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc7 = cixm1*ciym1*cizm1
                   cc8 = cix  *ciym1*cizm1
                   cc5 = cix  *ciy  *cizm1
                   cc6 = cixm1*ciy  *cizm1
                   cc3 = cixm1*ciym1*ciz
                   cc4 = cix  *ciym1*ciz
                   cc1 = cix  *ciy  *ciz
                   cc2 = cixm1*ciy  *ciz
                endif
!
                vc(nf7,1)=( cc7*vv(nc1,1) + cc8*vv(nc2,1) + cc5*vv(nc3,1) + &
                     cc6*vv(nc4,1) + cc3*vv(nc5,1) + cc4*vv(nc6,1) + &
                     cc1*vv(nc7,1) + cc2*vv(nc8,1))
                vc(nf7,2)=( cc7*vv(nc1,2) + cc8*vv(nc2,2) + cc5*vv(nc3,2) + &
                     cc6*vv(nc4,2) + cc3*vv(nc5,2) + cc4*vv(nc6,2) + &
                     cc1*vv(nc7,2) + cc2*vv(nc8,2))
                vc(nf7,3)=( cc7*vv(nc1,3) + cc8*vv(nc2,3) + cc5*vv(nc3,3) + &
                     cc6*vv(nc4,3) + cc3*vv(nc5,3) + cc4*vv(nc6,3) + &
                     cc1*vv(nc7,3) + cc2*vv(nc8,3))
                vc(nf7,4)=( cc7*vv(nc1,4) + cc8*vv(nc2,4) + cc5*vv(nc3,4) + &
                     cc6*vv(nc4,4) + cc3*vv(nc5,4) + cc4*vv(nc6,4) + &
                     cc1*vv(nc7,4) + cc2*vv(nc8,4))
                vc(nf7,5)=( cc7*vv(nc1,5) + cc8*vv(nc2,5) + cc5*vv(nc3,5) + &
                     cc6*vv(nc4,5) + cc3*vv(nc5,5) + cc4*vv(nc6,5) + &
                     cc1*vv(nc7,5) + cc2*vv(nc8,5))
             enddo
!
! ------------------------------------------------------------------------ C
! ---                           Cellule 8                              --- C
! ------------------------------------------------------------------------ C
!
             do ic = i1c,i2c-1
!
                nc1= indn_c(ic,jc,kc)
                nc2= nc1 + incc_i
                nc4= nc1 + incc_j
                nc3= nc4 + incc_i
                nc5= nc1 + incc_k
                nc6= nc2 + incc_k
                nc7= nc3 + incc_k
                nc8= nc4 + incc_k
!
                if = 2*ic - 1
                nf1 = indc_f(if,jf,kf)
                nf2 = nf1 + incf_i
                nf3 = nf1 + incf_ij
                nf4 = nf1 + incf_j
                nf5 = nf1 + incf_k
                nf6 = nf2 + incf_k
                nf7 = nf3 + incf_k
                nf8 = nf4 + incf_k
!
                if (ktrans.eq.3) then
! --Interplotion Trilineaire Volumique
                   cix = (            .5*vol(nf8) )/(vol(nf8)+vol(nf7))
                   ciy = ( vol(nf5) + .5*vol(nf8) )/(vol(nf8)+vol(nf5))
                   ciz = ( vol(nf4) + .5*vol(nf8) )/(vol(nf8)+vol(nf4))
!
                   cixm1 = 1-cix
                   ciym1 = 1-ciy
                   cizm1 = 1-ciz
!
                   cc8 = cixm1*ciym1*cizm1
                   cc7 = cix  *ciym1*cizm1
                   cc6 = cix  *ciy  *cizm1
                   cc5 = cixm1*ciy  *cizm1
                   cc4 = cixm1*ciym1*ciz
                   cc3 = cix  *ciym1*ciz
                   cc2 = cix  *ciy  *ciz
                   cc1 = cixm1*ciy  *ciz
                endif
!
                vc(nf8,1)=( cc8*vv(nc1,1) + cc7*vv(nc2,1) + cc6*vv(nc3,1) + &
                     cc5*vv(nc4,1) + cc4*vv(nc5,1) + cc3*vv(nc6,1) + &
                     cc2*vv(nc7,1) + cc1*vv(nc8,1))
                vc(nf8,2)=( cc8*vv(nc1,2) + cc7*vv(nc2,2) + cc6*vv(nc3,2) + &
                     cc5*vv(nc4,2) + cc4*vv(nc5,2) + cc3*vv(nc6,2) + &
                     cc2*vv(nc7,2) + cc1*vv(nc8,2))
                vc(nf8,3)=( cc8*vv(nc1,3) + cc7*vv(nc2,3) + cc6*vv(nc3,3) + &
                     cc5*vv(nc4,3) + cc4*vv(nc5,3) + cc3*vv(nc6,3) + &
                     cc2*vv(nc7,3) + cc1*vv(nc8,3))
                vc(nf8,4)=( cc8*vv(nc1,4) + cc7*vv(nc2,4) + cc6*vv(nc3,4) + &
                     cc5*vv(nc4,4) + cc4*vv(nc5,4) + cc3*vv(nc6,4) + &
                     cc2*vv(nc7,4) + cc1*vv(nc8,4))
                vc(nf8,5)=( cc8*vv(nc1,5) + cc7*vv(nc2,5) + cc6*vv(nc3,5) + &
                     cc5*vv(nc4,5) + cc4*vv(nc5,5) + cc3*vv(nc6,5) + &
                     cc2*vv(nc7,5) + cc1*vv(nc8,5))
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
  end subroutine smg_cf_3d
end module mod_smg_cf_3d
