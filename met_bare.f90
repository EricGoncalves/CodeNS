module mod_met_bare
  implicit none
contains
  subroutine met_bare( &
       l, &
       equat, &
       sn, &
       vol,s,mu,mut, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dtdx,dtdy,dtdz,bare, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!     calcul des termes bas Renolds de l'equation pour epsilon
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use mod_teq_grads
    implicit none
    integer          ::      i,    i1,    i2,  i2m1,  ind1
    integer          ::   ind2,     j,    j1,    j2,  j2m1
    integer          ::      k,    k1,    k2,  k2m1,     l
    integer          :: lgsnlt,     m,     n,   n0c,   nid
    integer          ::   nijd,   njd,  npsn
    double precision ::    bare(ip00),  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21)
    double precision ::   cmuk1(ip21),  cmuk2(ip21),   dtdx(ip00),   dtdy(ip00),   dtdz(ip00)
    double precision ::    dvxx(ip00),   dvxy(ip00),   dvxz(ip00),   dvyx(ip00),   dvyy(ip00)
    double precision ::    dvyz(ip00),   dvzx(ip00),   dvzy(ip00),   dvzz(ip00),     mu(ip12)
    double precision ::     mut(ip12), s(ip11,ip60),sn(ip31*ndir),    vol(ip11)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!

!
    n0c=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
!     initialisation
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=0.
          enddo
       enddo
    enddo
!--------------------------------------------------------------
!     calcul de vx
!------------------------------------------------------------
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvxx, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvxy, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvxz, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!--------------------------------------------------------------
!     calcul de vy
!---------------------------------------------------------------
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvyx, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvyy, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvyz, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!---------------------------------------------------------
!     calcul de vz
!--------------------------------------------------------
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvzx, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvzy, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call teq_grads( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         dvzz, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=bare(m)+dtdx(m)*dtdx(m) &
                  +dtdy(m)*dtdy(m) &
                  +dtdz(m)*dtdz(m)
          enddo
       enddo
    enddo
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0c
             bare(m)=2.*mut(n)*mu(n)*bare(m)/s(n,1)
          enddo
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_bare
end module mod_met_bare
