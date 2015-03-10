      subroutine smg_fcm( &
                 imgf,imgc, &
                 l, &
                 x,y,z)
!
!***********************************************************************
!
!     ACT
!_A    Extraction of coarse mesh from the fine mesh
!_A                       COARSE level : imgc
!_A                       FINE   level : imgf
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use chainecarac
implicit none
integer :: indc
integer :: indf
integer :: i
integer :: j
integer :: k
integer :: imgf
integer :: imgc
integer :: l
double precision :: x
double precision :: y
double precision :: z
integer :: i1c
integer :: i2c
integer :: ic
integer :: if
integer :: j1c
integer :: j2c
integer :: jc
integer :: jf
integer :: k1c
integer :: k2c
integer :: kc
integer :: kf
integer :: lmc
integer :: lmf
integer :: n0c
integer :: n0f
integer :: nc
integer :: nf
integer :: nic
integer :: nif
integer :: nijc
integer :: nijf
integer :: njc
integer :: njf
!
!-----------------------------------------------------------------------
!
      dimension x(ip21),y(ip21),z(ip21)

      indc(i,j,k)=n0c+1+(i-id1(lmc))+(j-jd1(lmc))*nic+(k-kd1(lmc))*nijc
      indf(i,j,k)=n0f+1+(i-id1(lmf))+(j-jd1(lmf))*nif+(k-kd1(lmf))*nijf
!
      lmc  = l+(imgc-1)*lz
      lmf  = l+(imgf-1)*lz
!
      n0c  = npn(lmc)
      nic  = id2(lmc)-id1(lmc)+1
      njc  = jd2(lmc)-jd1(lmc)+1
      nijc = nic*njc
!
      n0f  = npn(lmf)
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
      do kc = k1c,k2c
       kf = 2*kc -1
       if(equat(3:5).eq.'2dk' .or. equat(3:5).eq.'2xk') kf = kc
!
       do jc = j1c,j2c
        jf = 2*jc -1
        if(equat(3:5).eq.'2dj') jf = jc
!
!DEC$ IVDEP
        do ic = i1c,i2c
         if = 2*ic - 1
         nf = indf(if,jf,kf)
         nc = indc(ic,jc,kc)
!
         x(nc) = x(nf)
         y(nc) = y(nf)
         z(nc) = z(nf)
        enddo
       enddo
      enddo
!
      return
      end
