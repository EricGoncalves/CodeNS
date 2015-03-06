      subroutine met_fludcsa( &
                 l, &
                 s,mu, &
                 fd5x,fd5y,fd5z,fd6x,fd6y,fd6z)
!
!***********************************************************************
!
!     ACT
!_A   modele Spalart Allmaras
!_A   Calcul des densites de flux dissipatifs a partir de grad(nutilde)
!_A   mise a zero de fd6* pour resoudre 0=0
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    s          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!     OUT
!
!     I/O
!_/    fd5x       : arg real(ip12     )  ; comp x grad(nu tilde) puis flux diss.
!_/    fd5y       : arg real(ip12     )  ; comp y grad(nu tilde) puis flux diss.
!_/    fd5z       : arg real(ip12     )  ; comp z grad(nu tilde) puis flux diss.
!_/    fd6x       : arg real(ip12     )  ; comp x grad(e) puis flux diss.
!_/    fd6y       : arg real(ip12     )  ; comp y grad(e) puis flux diss.
!_/    fd6z       : arg real(ip12     )  ; comp z grad(e) puis flux diss.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use chainecarac
!
!-----------------------------------------------------------------------
!
      real mu
!
      dimension s(ip11,ip60)
      dimension mu(ip12), &
                fd5x(ip12),fd5y(ip12),fd5z(ip12), &
                fd6x(ip12),fd6y(ip12),fd6z(ip12)
!
      indc(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
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
!
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
!
      sigma1=1./sigma
!
      ind1=indc(i1  ,j1  ,k1  )
      ind2=indc(i2m1,j2m1,k2m1)
      do n=ind1,ind2
!        m=n-n0
         smutot=sigma1*(mu(n)+s(n,6))
         fd5x(n)=smutot*fd5x(n)
         fd5y(n)=smutot*fd5y(n)
         fd5z(n)=smutot*fd5z(n)
         fd6x(n)=0.
         fd6y(n)=0.
         fd6z(n)=0.
      enddo
!
      return
      end
