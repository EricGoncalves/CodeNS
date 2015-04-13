module mod_smg_fcr
  implicit none
contains
  subroutine smg_fcr( &
       imgf,imgc,df)
!
!***********************************************************************
!
!     ACT
!_A    Fine to Coarse  :  Residual Restriction from Fine --> Coarse
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
    double precision :: df(ip11,ip60)
!$OMP MASTER
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
!$OMP SIMD
                do ic = i1c,i2c-1
                   if = 2*ic - 1
                   nf1= indf(if,jf,kf)
                   nf2= nf1 + incf_i
                   nf4= nf1 + incf_j
                   nf3= nf4 + incf_i
!
                   nc = indc(ic,jc,kc)
!
                   df(nc,1) = df(nf1,1) + df(nf2,1) + df(nf3,1) + df(nf4,1)
                   df(nc,2) = df(nf1,2) + df(nf2,2) + df(nf3,2) + df(nf4,2)
                   df(nc,3) = df(nf1,3) + df(nf2,3) + df(nf3,3) + df(nf4,3)
                   df(nc,4) = df(nf1,4) + df(nf2,4) + df(nf3,4) + df(nf4,4)
                   df(nc,5) = df(nf1,5) + df(nf2,5) + df(nf3,5) + df(nf4,5)
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
!$OMP SIMD
                do ic = i1c,i2c-1
                   if = 2*ic - 1
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
                   df(nc,1) = df(nf1,1) + df(nf2,1) + df(nf3,1) + df(nf4,1) &
                        + df(nf5,1) + df(nf6,1) + df(nf7,1) + df(nf8,1)
                   df(nc,2) = df(nf1,2) + df(nf2,2) + df(nf3,2) + df(nf4,2) &
                        + df(nf5,2) + df(nf6,2) + df(nf7,2) + df(nf8,2)
                   df(nc,3) = df(nf1,3) + df(nf2,3) + df(nf3,3) + df(nf4,3) &
                        + df(nf5,3) + df(nf6,3) + df(nf7,3) + df(nf8,3)
                   df(nc,4) = df(nf1,4) + df(nf2,4) + df(nf3,4) + df(nf4,4) &
                        + df(nf5,4) + df(nf6,4) + df(nf7,4) + df(nf8,4)
                   df(nc,5) = df(nf1,5) + df(nf2,5) + df(nf3,5) + df(nf4,5) &
                        + df(nf5,5) + df(nf6,5) + df(nf7,5) + df(nf8,5)
!
                enddo
             enddo
          enddo
!
       else
          stop 'smg_fcr'
       endif
!
    enddo
!
!$OMP END MASTER
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
  end subroutine smg_fcr
end module mod_smg_fcr
