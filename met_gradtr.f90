module mod_met_gradtr
  implicit none
contains
  subroutine met_gradtr( &
       l, &
       sn, &
       vol,s,mu,mut, &
       tt,dtdx,dtdy,dtdz, &
       fd5x,fd5y,fd5z,fd6x,fd6y,fd6z, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_A   Calcul des gradients de k et de la seconde variable de transport.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    s          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!
!     OUT
!_O    fd5x       : arg real(ip12     )  ; comp x grad(k)
!_O    fd5y       : arg real(ip12     )  ; comp y grad(k)
!_O    fd5z       : arg real(ip12     )  ; comp z grad(k)
!_O    fd6x       : arg real(ip12     )  ; comp x grad(e)
!_O    fd6y       : arg real(ip12     )  ; comp y grad(e)
!_O    fd6z       : arg real(ip12     )  ; comp z grad(e)
!
!     LOC
!_L    tt         : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  tt,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  tt,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  tt,z
!
!***********************************************************************
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use mod_teq_grads
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2m1
    integer          ::   imax,  imin,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,    j2,  j2m1,  jmax
    integer          ::   jmin,     k,    k1,  k1m1,    k2
    integer          ::   k2m1,  kmax,  kmin,     l,lgsnlt
    integer          ::      m,     n,    n0,   nid,  nijd
    integer          ::    njd,  npsn
    double precision :: cmui1,cmui2,cmuj1,cmuj2,cmuk1
    double precision :: cmuk2, dtdx, dtdy, dtdz, fd5x
    double precision ::  fd5y, fd5z, fd6x, fd6y, fd6z
    double precision ::    mu,  mut,    s,   sn,   tt
    double precision ::   vol
!
!-----------------------------------------------------------------------
!
    dimension sn(ip31*ndir), &
         vol(ip11)
    dimension s(ip11,ip60)
    dimension mut(ip12),mu (ip12), &
         fd5x(ip12),fd5y(ip12),fd5z(ip12), &
         fd6x(ip12),fd6y(ip12),fd6z(ip12)
    dimension tt(ip00),dtdx(ip00),dtdy(ip00),dtdz(ip00)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
!


    n0=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    imin=i1m1
    imax=i2
    jmin=j1m1
    jmax=j2
    kmin=k1m1
    kmax=k2
!
    if (equat(3:5).eq.'2di') then
       imin=i1
       imax=i2m1
    endif
    if (equat(3:5).eq.'2dj') then
       jmin=j1
       jmax=j2m1
    endif
    if (equat(3:5).eq.'2dk') then
       kmin=k1
       kmax=k2m1
    endif
!---------------------------------------------------------
!     grad(k)
!--------------------------------------------------------
    do k=kmin,kmax
       do j=jmin,jmax
          ind1=indc(imin,j,k)
          ind2=indc(imax,j,k)
          do n=ind1,ind2
             m=n-n0
             tt(m)= s(n,6)/s(n,1)
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
         tt, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    ind1=indc(i1  ,j1  ,k1  )
    ind2=indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       m=n-n0
       fd5x(n)=dtdx(m)
       fd5y(n)=dtdy(m)
       fd5z(n)=dtdz(m)
    enddo
!--------------------------------------------------------
!     grad(epsilon)
!--------------------------------------------------------
    do k=kmin,kmax
       do j=jmin,jmax
          ind1=indc(imin,j,k)
          ind2=indc(imax,j,k)
          do n=ind1,ind2
             m=n-n0
             tt(m)= s(n,7)/s(n,1)
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
         tt, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)

    ind1=indc(i1  ,j1  ,k1  )
    ind2=indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       m=n-n0
       fd6x(n)=dtdx(m)
       fd6y(n)=dtdy(m)
       fd6z(n)=dtdz(m)
    enddo

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_gradtr
end module mod_met_gradtr
