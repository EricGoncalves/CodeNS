module mod_smg_fcv
  implicit none
contains
  subroutine smg_fcv( &
       imgf,imgc, &
       vol,v)
!
!***********************************************************************
!
!     ACT
!_A    Fine to Coarse  :  Transfer of Variables from Fine --> Coarse
!_A                       COARSE level : imgc
!_A                       FINE   level : imgf
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use schemanum
    implicit none
    integer          ::      i,   i1c,   i2c,    ic,    id
    integer          ::     if,  imgc,  imgf,incf_i,incf_j
    integer          :: incf_k,     j,   j1c,   j2c,    jc
    integer          ::     jd,    jf,     k,   k1c,   k2c
    integer          ::     kc,    kd,    kf,     l,   lmc
    integer          ::    lmf,   n0c,   n0f,    nc,   nf1
    integer          ::    nf2,   nf3,   nf4,   nf5,   nf6
    integer          ::    nf7,   nf8,   nic,   nif,  nijc
    integer          ::   nijf,   njc,   njf
    double precision :: v(ip11,ip60),   vol(ip11)
!
!-----------------------------------------------------------------------
!
!



!
    do l = 1,lzx
       lmc  = l+(imgc-1)*lz
       lmf  = l+(imgf-1)*lz
!
       n0c  = npc(lmc)
       nic  = id2(lmc)-id1(lmc)+1
       njc  = jd2(lmc)-jd1(lmc)+1
       nijc = nic*njc
!
       n0f  = npc(lmf)
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
       incf_i  = incf(1,0,0)
       incf_j  = incf(0,1,0)
       incf_k  = incf(0,0,1)
!
       if(equat(3:5).eq.'2dk') then
!
          do kc = k1c,k2c-1
             kf = kc
             do jc = j1c,j2c-1
                jf = 2*jc -1
                do ic = i1c,i2c-1
                   if = 2*ic - 1
                   nf1= indf(if,jf,kf)
                   nf2= nf1 + incf_i
                   nf4= nf1 + incf_j
                   nf3= nf4 + incf_i
                   nc = indc(ic,jc,kc)
!
                   v(nc,1) = (vol(nf1)*v(nf1,1) + vol(nf2)*v(nf2,1) &
                        +vol(nf3)*v(nf3,1) + vol(nf4)*v(nf4,1))/vol(nc)
                   v(nc,2) = (vol(nf1)*v(nf1,2) + vol(nf2)*v(nf2,2) &
                        +vol(nf3)*v(nf3,2) + vol(nf4)*v(nf4,2))/vol(nc)
                   v(nc,3) = (vol(nf1)*v(nf1,3) + vol(nf2)*v(nf2,3) &
                        +vol(nf3)*v(nf3,3) + vol(nf4)*v(nf4,3))/vol(nc)
                   v(nc,4) = (vol(nf1)*v(nf1,4) + vol(nf2)*v(nf2,4) &
                        +vol(nf3)*v(nf3,4) + vol(nf4)*v(nf4,4))/vol(nc)
                   v(nc,5) = (vol(nf1)*v(nf1,5) + vol(nf2)*v(nf2,5) &
                        +vol(nf3)*v(nf3,5) + vol(nf4)*v(nf4,5))/vol(nc)
                enddo
             enddo
          enddo
!
       else if(equat(3:5).eq.'2dj') then
          stop 'MG2dj'
!
       else if(equat(3:5).eq.'2di') then
          stop 'MG2di'
!
       else if(equat(3:4).eq.'3d') then
!
          do kc = k1c,k2c-1
             kf = 2*kc -1
             do jc = j1c,j2c-1
                jf = 2*jc -1
                do ic = i1c,i2c-1
                   if = 2*ic - 1
!
                   nf1= indf(if,jf,kf)
                   nf2= nf1 + incf_i
                   nf4= nf1 + incf_j
                   nf3= nf4 + incf_i
!
                   nf5= nf1 + incf_k
                   nf6= nf2 + incf_k
                   nf7= nf3 + incf_k
                   nf8= nf4 + incf_k
!
                   nc = indc(ic,jc,kc)
!
                   v(nc,1) = (vol(nf1)*v(nf1,1) + vol(nf2)*v(nf2,1) &
                        +vol(nf3)*v(nf3,1) + vol(nf4)*v(nf4,1) &
                        +vol(nf5)*v(nf5,1) + vol(nf6)*v(nf6,1) &
                        +vol(nf7)*v(nf7,1) + vol(nf8)*v(nf8,1))/vol(nc)
                   v(nc,2) = (vol(nf1)*v(nf1,2) + vol(nf2)*v(nf2,2) &
                        +vol(nf3)*v(nf3,2) + vol(nf4)*v(nf4,2) &
                        +vol(nf5)*v(nf5,2) + vol(nf6)*v(nf6,2) &
                        +vol(nf7)*v(nf7,2) + vol(nf8)*v(nf8,2))/vol(nc)
                   v(nc,3) = (vol(nf1)*v(nf1,3) + vol(nf2)*v(nf2,3) &
                        +vol(nf3)*v(nf3,3) + vol(nf4)*v(nf4,3) &
                        +vol(nf5)*v(nf5,3) + vol(nf6)*v(nf6,3) &
                        +vol(nf7)*v(nf7,3) + vol(nf8)*v(nf8,3))/vol(nc)
                   v(nc,4) = (vol(nf1)*v(nf1,4) + vol(nf2)*v(nf2,4) &
                        +vol(nf3)*v(nf3,4) + vol(nf4)*v(nf4,4) &
                        +vol(nf5)*v(nf5,4) + vol(nf6)*v(nf6,4) &
                        +vol(nf7)*v(nf7,4) + vol(nf8)*v(nf8,4))/vol(nc)
                   v(nc,5) = (vol(nf1)*v(nf1,5) + vol(nf2)*v(nf2,5) &
                        +vol(nf3)*v(nf3,5) + vol(nf4)*v(nf4,5) &
                        +vol(nf5)*v(nf5,5) + vol(nf6)*v(nf6,5) &
                        +vol(nf7)*v(nf7,5) + vol(nf8)*v(nf8,5))/vol(nc)
                enddo
             enddo
          enddo
!
       else
          stop 'smg_fcv'
       endif
!
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(lmc))+(j-jd1(lmc))*nic+(k-kd1(lmc))*nijc
    end function indc
    function    indf(i,j,k)
      implicit none
      integer          ::    i,indf,   j,   k
      indf=n0f+1+(i-id1(lmf))+(j-jd1(lmf))*nif+(k-kd1(lmf))*nijf
    end function indf
    function    incf(id,jd,kd)
      implicit none
      integer          ::   id,incf,  jd,  kd
      incf=id+jd*nif+kd*nijf
    end function incf
  end subroutine smg_fcv
end module mod_smg_fcv
