module mod_met_roe2o
  implicit none
contains
  subroutine met_roe2o( &
       l,t,d, &
       sn,lgsnlt, &
       vol,del6,del7)
!
!***********************************************************************
!
!_DA :  janvier 2002 - Eric Goncalves / SINUMEF
!
!     Schema decentre de Roe - ordre 2 - Methode Flux Splitting
!     Pas de correction de Harten
!     Stockage de la dissipation dans le tableau d
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    implicit none
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,  i2p1,    id,  ind1,  ind2
    integer          ::     is,     j,    j1,  j1m1,  j1p1
    integer          ::     j2,  j2m1,  j2p1,    jd,     k
    integer          ::     k1,  k1m1,  k1p1,    k2,  k2m1
    integer          ::   k2p1,    kd,     l,lgsnlt,     m
    integer          ::     m1,    m2,     n,    n0,    n1
    integer          ::    nci,   ncj,   nck,   nid,  nijd
    integer          ::   ninc,   njd
    double precision ::                    a,                   b,        d(ip11,ip60),                  dd,          del6(ip00)
    double precision ::           del7(ip00),                rlam,                 rro,                rro1,sn(lgsnlt,nind,ndir)
    double precision ::         t(ip11,ip60),                   u,                   v,           vol(ip11),                   w
    double precision ::                    x,                   y,                   z
!
!-----------------------------------------------------------------------
!
!


!



    n0 = npc(l)
    i1 = ii1(l)
    i2 = ii2(l)
    j1 = jj1(l)
    j2 = jj2(l)
    k1 = kk1(l)
    k2 = kk2(l)
!
    nid  = id2(l)-id1(l)+1
    njd  = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    i1m1 = i1-1
    j1m1 = j1-1
    k1m1 = k1-1
    i1p1 = i1+1
    j1p1 = j1+1
    k1p1 = k1+1
    i2m1 = i2-1
    j2m1 = j2-1
    k2m1 = k2-1
    i2p1 = i2+1
    j2p1 = j2+1
    k2p1 = k2+1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
!
! -----initialisation-----------------------------------------
!
    do k = k1m1,k2
       do j = j1m1,j2
          ind1 = ind(i1m1,j,k)+n0
          ind2 = ind(i2  ,j,k)+n0
!$OMP SIMD
          do n = ind1,ind2
             d(n,6)=0.
             d(n,7)=0.
          enddo
       enddo
    enddo
!
! ------------------------------------------------------------
!
    if(equat(3:5).ne.'2dk') then
!
!     (direction k)
!
       ninc=nck
!
       do k = k1,k2
          do j = j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n = m+n0
                n1 = n-ninc
                del6(m)=(t(n,6)-t(n1,6))
                del7(m)=(t(n,7)-t(n1,7))
             enddo
          enddo
       enddo
!
       is=-1
       do k = k1m1,k2p1,(k2p1-k1m1)
          is=-is
          do j = j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n = m+n0
                n1 = n-ninc
                del6(m)=del6(m+is*ninc)
                del7(m)=del7(m+is*ninc)
             enddo
          enddo
       enddo
!
       do k = k1,k2
          do j = j1,j2m1
             ind1 = ind(i1,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n=m+n0
                m1=m-ninc
                m2=m+ninc
                n1=n-ninc
                rro=2*t(n ,1)
                rro1=2*t(n1,1)
                u=t(n,2)/rro + t(n1,2)/rro1
                v=t(n,3)/rro + t(n1,3)/rro1
                w=t(n,4)/rro + t(n1,4)/rro1
                rlam=u*sn(m,3,1)+v*sn(m,3,2)+w*sn(m,3,3)
                dd=0.5*abs(rlam)*(del6(m)-fi1(del6(m1),del6(m),del6(m2)))
                d(n     ,6) = d(n     ,6)-dd
                d(n-ninc,6) = d(n-ninc,6)+dd
                dd=0.5*abs(rlam)*(del7(m)-fi1(del7(m1),del7(m),del7(m2)))
                d(n     ,7) = d(n     ,7)-dd
                d(n-ninc,7) = d(n-ninc,7)+dd
             enddo
          enddo
       enddo
!
    endif
!
! ------------------------------------------------------------
    if(equat(3:5).ne.'2dj') then
!
       ninc=ncj
!
!     (direction j)
!
       do k = k1,k2m1
          do j = j1,j2
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n       = m+n0
                n1      = n-ninc
                del6(m)=(t(n,6)-t(n1,6))
                del7(m)=(t(n,7)-t(n1,7))
             enddo
          enddo
       enddo
!
       do k = k1,k2m1
          is=-1
          do j = j1m1,j2p1,(j2p1-j1m1)
             is=-is
             ind1 = ind(i1,  j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                del6(m)=del6(m+is*ninc)
                del7(m)=del7(m+is*ninc)
             enddo
          enddo
       enddo
!
       do k = k1,k2m1
          do j = j1,j2
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n       = m+n0
                m1      = m-ninc
                m2      = m+ninc
                n1      = n-ninc
                rro     = 2*t(n ,1)
                rro1    = 2*t(n1,1)
                u       = t(n,2)/rro + t(n1,2)/rro1
                v       = t(n,3)/rro + t(n1,3)/rro1
                w       = t(n,4)/rro + t(n1,4)/rro1
                rlam=u*sn(m,2,1)+v*sn(m,2,2)+w*sn(m,2,3)
                dd=0.5*abs(rlam)*(del6(m)-fi1(del6(m1),del6(m),del6(m2)))
                d(n     ,6) = d(n     ,6)-dd
                d(n-ninc,6) = d(n-ninc,6)+dd
                dd=0.5*abs(rlam)*(del7(m)-fi1(del7(m1),del7(m),del7(m2)))
                d(n     ,7) = d(n     ,7)-dd
                d(n-ninc,7) = d(n-ninc,7)+dd
             enddo
          enddo
       enddo
!
    endif
!
! ------------------------------------------------------------
    if(equat(3:5).ne.'2di') then
!
       ninc=nci
!
!     (direction i)
!
       do k = k1,k2m1
          do j = j1,j2
             ind1 = ind(i1,j,k)
             ind2 = ind(i2,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n       = m+n0
                n1      = n-ninc
                del6(m)=(t(n,6)-t(n1,6))
                del7(m)=(t(n,7)-t(n1,7))
             enddo
          enddo
       enddo
!
       is=-1
       do i = i1m1,i2p1,(i2p1-i1m1)
          is=-is
          do k = k1,k2m1
             ind1 = ind(i,j1  ,k)
             ind2 = ind(i,j2m1,k)
!$OMP SIMD
             do m = ind1,ind2,ncj
                del6(m)=del6(m+is*ninc)
                del7(m)=del7(m+is*ninc)
             enddo
          enddo
       enddo
!
       do k = k1,k2m1
          do j = j1,j2m1
             ind1 = ind(i1,j,k)
             ind2 = ind(i2,j,k)
!$OMP SIMD
             do m = ind1,ind2
                n       = m+n0
                m1      = m-ninc
                m2      = m+ninc
                n1      = n-ninc
                rro     = 2*t(n ,1)
                rro1    = 2*t(n1,1)
                u       = t(n,2)/rro + t(n1,2)/rro1
                v       = t(n,3)/rro + t(n1,3)/rro1
                w       = t(n,4)/rro + t(n1,4)/rro1
                rlam=u*sn(m,1,1)+v*sn(m,1,2)+w*sn(m,1,3)
                dd=0.5*abs(rlam)*(del6(m)-fi1(del6(m1),del6(m),del6(m2)))
                d(n     ,6) = d(n     ,6)-dd
                d(n-ninc,6) = d(n-ninc,6)+dd
                dd=0.5*abs(rlam)*(del7(m)-fi1(del7(m1),del7(m),del7(m2)))
                d(n     ,7) = d(n     ,7)-dd
                d(n-ninc,7) = d(n-ninc,7)+dd
             enddo
          enddo
       enddo
!
    endif

    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
    function    amimd(a,b)
      implicit none
      double precision ::     a,amimd,    b
      amimd=sign(1.,a)*max(0.,min(abs(a),b*sign(1.,a)))
    end function amimd
    function    fi1(x,y,z)
      implicit none
      double precision :: fi1,  x,  y,  z
      fi1=amimd(x,amimd(y,z))
    end function fi1
  end subroutine met_roe2o
end module mod_met_roe2o
