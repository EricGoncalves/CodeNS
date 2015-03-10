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
integer :: incf
integer :: indc
integer :: indf
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: imgf
integer :: imgc
double precision :: vol
double precision :: v
integer :: i1c
integer :: i2c
integer :: ic
integer :: if
integer :: incf_i
integer :: incf_j
integer :: incf_k
integer :: j1c
integer :: j2c
integer :: jc
integer :: jf
integer :: k1c
integer :: k2c
integer :: kc
integer :: kf
integer :: l
integer :: lmc
integer :: lmf
integer :: n0c
integer :: n0f
integer :: nc
integer :: nf1
integer :: nf2
integer :: nf3
integer :: nf4
integer :: nf5
integer :: nf6
integer :: nf7
integer :: nf8
integer :: nic
integer :: nif
integer :: nijc
integer :: nijf
integer :: njc
integer :: njf
!
!-----------------------------------------------------------------------
!
      dimension vol(ip11),v(ip11)
!
      indc(i,j,k)=n0c+1+(i-id1(lmc))+(j-jd1(lmc))*nic+(k-kd1(lmc))*nijc
      indf(i,j,k)=n0f+1+(i-id1(lmf))+(j-jd1(lmf))*nif+(k-kd1(lmf))*nijf
      incf(id,jd,kd)=id+jd*nif+kd*nijf
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
!DEC$ IVDEP
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
!DEC$ IVDEP
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
      end
