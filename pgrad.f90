module mod_pgrad
  implicit none
contains
  subroutine pgrad(  &
       sn,vol,lgsnlt,l, &
       ps,dpdx,dpdy,dpdz)
!
!***********************************************************************
!
!_DA  DATE_C : mai 1999 -- AUTEUR : DMAE - Eric Goncalves
!
!     ACT
!_A    Calcul du gradient de la pression.
!
!     OUT
!_O    dpdx        : arg real(ip00      ) ; composante x du gradient de pression
!_O    dpdy        : arg real(ip00      ) ; composante y du gradient de pression
!_O    dpdz        : arg real(ip00      ) ; composante z du gradient de pression
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use chainecarac
    implicit none
    integer          ::      i,    i1,  i1m1,  i1p1,    i2
    integer          ::   i2m1,  i2p1,    id,  ind1,  ind2
    integer          ::      j,    j1,  j1m1,  j1p1,    j2
    integer          ::   j2m1,  j2p1,    jd,     k,    k1
    integer          ::   k1m1,  k1p1,    k2,  k2m1,  k2p1
    integer          ::     kd,     l,lgsnlt,     m,     n
    integer          ::     n0,   nci,   ncj,   nck,   nid
    integer          ::   nijd,   njd
    double precision ::                  cc0,          dpdx(ip00),          dpdy(ip00),          dpdz(ip00),                 eps
    double precision ::             ps(ip11),sn(lgsnlt,nind,ndir),                  ts,           vol(ip11),                vols
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!

!


!
    eps=0.001
    n0=npc(l)
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
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    i2p1=i2+1
    j1m1=j1-1
    j2p1=j2+1
    k1m1=k1-1
    k2p1=k2+1
!
    nci=inc(1,0,0)
    ncj=inc(0,1,0)
    nck=inc(0,0,1)
!
!     initialisation
    ind1=ind(id1(l),jd1(l),kd1(l))-n0
    ind2=ind(id2(l),jd2(l),kd2(l))-n0
!$OMP SIMD
    do m=ind1,ind2
       dpdx(m)=0.
       dpdy(m)=0.
       dpdz(m)=0.
    enddo
!
!----Calcul du gradient de pression----------------------------
!
!     direction k (a travers les facettes k=cste)
    if (equat(3:5).ne.'2dk') then
!
       do k=k1p1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                dpdx(m    )=dpdx(m    )-sn(m,3,1)*(ps(n)+ps(n-nck))
                dpdy(m    )=dpdy(m    )-sn(m,3,2)*(ps(n)+ps(n-nck))
                dpdz(m    )=dpdz(m    )-sn(m,3,3)*(ps(n)+ps(n-nck))
                dpdx(m-nck)=dpdx(m-nck)+sn(m,3,1)*(ps(n)+ps(n-nck))
                dpdy(m-nck)=dpdy(m-nck)+sn(m,3,2)*(ps(n)+ps(n-nck))
                dpdz(m-nck)=dpdz(m-nck)+sn(m,3,3)*(ps(n)+ps(n-nck))
             enddo
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = ind(i1  ,j,k1)
          ind2 = ind(i2m1,j,k1)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0
             dpdx(m)=dpdx(m)-sn(m,3,1)*2*ps(n-nck)
             dpdy(m)=dpdy(m)-sn(m,3,2)*2*ps(n-nck)
             dpdz(m)=dpdz(m)-sn(m,3,3)*2*ps(n-nck)
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = ind(i1  ,j,k2)
          ind2 = ind(i2m1,j,k2)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0
             dpdx(m-nck)=dpdx(m-nck)+sn(m,3,1)*2*ps(n)
             dpdy(m-nck)=dpdy(m-nck)+sn(m,3,2)*2*ps(n)
             dpdz(m-nck)=dpdz(m-nck)+sn(m,3,3)*2*ps(n)
          enddo
       enddo
!
    endif
!
!     direction j (a travers les facettes j=cste)
    if (equat(3:5).ne.'2dj') then
!
       do j=j1p1,j2m1
          do k=k1,k2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                dpdx(m    )=dpdx(m    )-sn(m,2,1)*(ps(n)+ps(n-ncj))
                dpdy(m    )=dpdy(m    )-sn(m,2,2)*(ps(n)+ps(n-ncj))
                dpdz(m    )=dpdz(m    )-sn(m,2,3)*(ps(n)+ps(n-ncj))
                dpdx(m-ncj)=dpdx(m-ncj)+sn(m,2,1)*(ps(n)+ps(n-ncj))
                dpdy(m-ncj)=dpdy(m-ncj)+sn(m,2,2)*(ps(n)+ps(n-ncj))
                dpdz(m-ncj)=dpdz(m-ncj)+sn(m,2,3)*(ps(n)+ps(n-ncj))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1  ,j1,k)
          ind2 = ind(i2m1,j1,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0
             dpdx(m)=dpdx(m)-sn(m,2,1)*2*ps(n-ncj)
             dpdy(m)=dpdy(m)-sn(m,2,2)*2*ps(n-ncj)
             dpdz(m)=dpdz(m)-sn(m,2,3)*2*ps(n-ncj)
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1  ,j2,k)
          ind2 = ind(i2m1,j2,k)
!$OMP SIMD
          do n=ind1,ind2
             m=n-n0
             dpdx(m-ncj)=dpdx(m-ncj)+sn(m,2,1)*2*ps(n)
             dpdy(m-ncj)=dpdy(m-ncj)+sn(m,2,2)*2*ps(n)
             dpdz(m-ncj)=dpdz(m-ncj)+sn(m,2,3)*2*ps(n)
          enddo
       enddo
!
    end if
!
!     direction i (a travers les facettes i=cste)
    if (equat(3:5).ne.'2di') then
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1p1,j,k)
             ind2 = ind(i2m1,j,k)
!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                dpdx(m    )=dpdx(m    )-sn(m,1,1)*(ps(n)+ps(n-nci))
                dpdy(m    )=dpdy(m    )-sn(m,1,2)*(ps(n)+ps(n-nci))
                dpdz(m    )=dpdz(m    )-sn(m,1,3)*(ps(n)+ps(n-nci))
                dpdx(m-nci)=dpdx(m-nci)+sn(m,1,1)*(ps(n)+ps(n-nci))
                dpdy(m-nci)=dpdy(m-nci)+sn(m,1,2)*(ps(n)+ps(n-nci))
                dpdz(m-nci)=dpdz(m-nci)+sn(m,1,3)*(ps(n)+ps(n-nci))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1,j1  ,k)
          ind2 = ind(i1,j2m1,k)
          do n=ind1,ind2,ncj
             m=n-n0
             dpdx(m)=dpdx(m)-sn(m,1,1)*2*ps(n-nci)
             dpdy(m)=dpdy(m)-sn(m,1,2)*2*ps(n-nci)
             dpdz(m)=dpdz(m)-sn(m,1,3)*2*ps(n-nci)
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i2,j1  ,k)
          ind2 = ind(i2,j2m1,k)
          do  n=ind1,ind2,ncj
             m=n-n0
             dpdx(m-nci)=dpdx(m-nci)+sn(m,1,1)*2*ps(n)
             dpdy(m-nci)=dpdy(m-nci)+sn(m,1,2)*2*ps(n)
             dpdz(m-nci)=dpdz(m-nci)+sn(m,1,3)*2*ps(n)
          enddo
       enddo
!
    endif
!
!-----calcul du gradient :
!
    ind1 = ind(i1  ,j1  ,k1  )
    ind2 = ind(i2m1,j2m1,k2m1)
!$OMP SIMD
    do n=ind1,ind2
       m=n-n0
!       le coefficient 1/2 provient de la moyenne de ps
       ts=sign(0.5,-vol(n))
       vols=(0.5+ts)*eps+(0.5-ts)*vol(n)
       cc0=0.5/vols
!
       dpdx(m)=dpdx(m)*cc0
       dpdy(m)=dpdy(m)*cc0
       dpdz(m)=dpdz(m)*cc0
    enddo
!
!$OMP END MASTER
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine pgrad
end module mod_pgrad
