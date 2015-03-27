module mod_smg_fcs
  implicit none
contains
  subroutine smg_fcs( &
       imgf,imgc, &
       vol,v)
!
!***********************************************************************
!
!     ACT
!_A    Fine to Coarse  :  Transfer of scalar mut from Fine --> Coarse
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
    double precision ::   v(ip11),vol(ip11)
!
!-----------------------------------------------------------------------
!
!



!
    do l = 1,lzx
!
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
!!$OMP SIMD
                do ic = i1c,i2c-1
                   if = 2*ic - 1
!
                   nf1= indf(if,jf,kf)
                   nf2= nf1 + incf_i
                   nf4= nf1 + incf_j
                   nf3= nf4 + incf_i
                   nc = indc(ic,jc,kc)
!
                   v(nc) = (vol(nf1)*v(nf1) + vol(nf2)*v(nf2) &
                        +vol(nf3)*v(nf3) + vol(nf4)*v(nf4))/vol(nc)
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
!!$OMP SIMD
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
                   v(nc)=(vol(nf1)*v(nf1) + vol(nf2)*v(nf2) &
                        +vol(nf3)*v(nf3) + vol(nf4)*v(nf4) &
                        +vol(nf5)*v(nf5) + vol(nf6)*v(nf6) &
                        +vol(nf7)*v(nf7) + vol(nf8)*v(nf8))/vol(nc)
                enddo
             enddo
          enddo
!
       else
          stop 'smg_fcs'
       end if
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
  end subroutine smg_fcs
end module mod_smg_fcs
