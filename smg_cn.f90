module mod_smg_cn
implicit none
contains
      subroutine smg_cn( &
                 img, &
                 vol,volt, &
                 vc,vv)
!
!***********************************************************************
!
!     ACT
!_A     Projection of Cell-Centered values to Nodal values
!_A     Transfer of Corrections
!_A     Computation carried out for level -> img
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use chainecarac
implicit none
integer :: inc
integer :: indc
integer :: indn
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: img
double precision :: vol
double precision :: volt
double precision :: vc
double precision :: vv
double precision :: eps
integer :: i1
integer :: i2
integer :: iinf
integer :: inc_i
integer :: inc_ij
integer :: inc_ijk
integer :: inc_ik
integer :: inc_j
integer :: inc_jk
integer :: inc_k
integer :: isup
integer :: j1
integer :: j2
integer :: jinf
integer :: jsup
integer :: k1
integer :: k2
integer :: kinf
integer :: ksup
integer :: l
integer :: lm
integer :: n
integer :: n0c
integer :: n0n
integer :: nc
integer :: ndeb
integer :: nfin
integer :: ni
integer :: nij
integer :: nj
integer :: nn
double precision :: ts
double precision :: vols
!
!-----------------------------------------------------------------------
!
      dimension vc(ip11,ip60),vv(ip11,ip60)
      dimension vol(ip11),volt(ip11)
!

!
      indn(i,j,k)=n0n+1+(i-id1(lm))+(j-jd1(lm))*ni+(k-kd1(lm))*nij
      indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*ni+(k-kd1(lm))*nij
      inc(id,jd,kd)=id+jd*ni+kd*nij
      eps=0.001
!
      do l = 1,lzx
       lm = l+(img-1)*lz
       ndeb = npn(lm)+1
       nfin = npn(lm)+nnn(lm)
!
       do n = ndeb,nfin
        vv(n,1)=0.
        vv(n,2)=0.
        vv(n,3)=0.
        vv(n,4)=0.
        vv(n,5)=0.
        volt(n)=0.
       enddo
!
       n0n = npn(lm)
       n0c = npc(lm)
       ni  = id2(lm)-id1(lm)+1
       nj  = jd2(lm)-jd1(lm)+1
       nij = ni*nj
!
       i1 =ii1(lm)
       i2 =ii2(lm)
       j1 =jj1(lm)
       j2 =jj2(lm)
       k1 =kk1(lm)
       k2 =kk2(lm)
!
       inc_i  = inc(1,0,0)
       inc_j  = inc(0,1,0)
       inc_k  = inc(0,0,1)
       inc_ij = inc(1,1,0)
       inc_jk = inc(0,1,1)
       inc_ik = inc(1,0,1)
       inc_ijk= inc(1,1,1)
!
       if(equat(3:5).eq.'2dk') then
        kinf = k1
        ksup = k2-1
       else
        kinf = k1-1
        ksup = k2
       endif
!
       if(equat(3:5).eq.'2dj') then
        jinf = j1
        jsup = j2-1
       else
        jinf = j1-1
        jsup = j2
       endif
!
       iinf = i1-1
       isup = i2
!
       do k=kinf,ksup
        do j=jinf,jsup
!
         do i=iinf,isup
          nn = indn(i,j,k)
          nc = indc(i,j,k)
          vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
          vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
          vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
          vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
          vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
          volt(nn)  = volt(nn)+vol(nc)
         enddo
!
         do i=iinf,isup
          nn = indn(i,j,k)+inc_i
          nc = indc(i,j,k)
          vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
          vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
          vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
          vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
          vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
          volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_j
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_ij
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_k
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_ik
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_jk
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
        do i=iinf,isup
         nn = indn(i,j,k)+inc_ijk
         nc = indc(i,j,k)
         vv(nn,1)  = vv(nn,1) +vc(nc,1)*vol(nc)
         vv(nn,2)  = vv(nn,2) +vc(nc,2)*vol(nc)
         vv(nn,3)  = vv(nn,3) +vc(nc,3)*vol(nc)
         vv(nn,4)  = vv(nn,4) +vc(nc,4)*vol(nc)
         vv(nn,5)  = vv(nn,5) +vc(nc,5)*vol(nc)
         volt(nn)  = volt(nn)+vol(nc)
        enddo
!
       enddo
      enddo
!
      do n = ndeb,nfin
       ts=sign(0.5,-volt(n))
       vols = (0.5+ts)*eps+(0.5-ts)*volt(n)+tiny(1.)
       vv(n,1)=vv(n,1)/vols
       vv(n,2)=vv(n,2)/vols
       vv(n,3)=vv(n,3)/vols
       vv(n,4)=vv(n,4)/vols
       vv(n,5)=vv(n,5)/vols
      enddo
!
      enddo
!
      return
      end subroutine
end module
