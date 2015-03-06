      subroutine met_dual2( &
                 icycle,l,u,v, &
                 vol,ptdual)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 2002 : Eric GONCALVES - SINUMEF
!
!     ACT
!_A    Calcul du residu instationnaire R* pour les equations
!-A    de transport de la turbulence - ordre 3
!
!-----------------------------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use schemanum
!
!-----------------------------------------------------------------------
!
      integer icycle,l
      real dti,fact
      dimension v(ip11,ip60),u(ip11,ip60),ptdual(ip11,ip60)
      dimension vol(ip11)
!
      indc(i,j,k)=1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      fact = 11./6.
      if(icycle.eq.1) fact = 1.
      dti = 1./dt1min
!
      n0c  = npc(l)
      i1   = ii1(l)
      i2   = ii2(l)
      j1   = jj1(l)
      j2   = jj2(l)
      k1   = kk1(l)
      k2   = kk2(l)
      nid  = id2(l)-id1(l)+1
      njd  = jd2(l)-jd1(l)+1
      nijd = nid*njd
      i2m1 = i2-1
      j2m1 = j2-1
      k2m1 = k2-1
!
      do k=k1,k2m1
       do j=j1,j2m1
        ind1=indc(i1,j,k)
        ind2=indc(i2m1,j,k)
        do n=ind1,ind2
         c0=vol(n)*dti
         u(n,6)=u(n,6) + c0*(fact*v(n,6)+ptdual(n,6))
         u(n,7)=u(n,7) + c0*(fact*v(n,7)+ptdual(n,7))
         enddo
        enddo
      enddo
!
      return
      end
