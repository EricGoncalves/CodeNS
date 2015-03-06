      subroutine met_brko( &
                 l, &
                 mut,v, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
!***********************************************************************
!
!     ACT
!_A   initialisation des grandeurs k et omega pour le modele de Menter
!_A   a partir de mu_t (u'v'/k=0.3)
!
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
      real mut,mutmx
      dimension v(ip11,ip60)
      dimension mut(ip12)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00)
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
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      if(equatt(1:3).eq.'2KO') then
         cmu1=1./sqrt(betas)
      else
         cmu1=1./sqrt(betae)
      endif
!
      do k=k1,k2m1
       do j=j1,j2m1
        ind1=indc(i1,j,k)
        ind2=indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0
         rota = sqrt((dvzy(m)-dvyz(m))**2 &
                    +(dvxz(m)-dvzx(m))**2 &
                    +(dvyx(m)-dvxy(m))**2)
         mutmx=max(mut(n),1.e-12)
         v(n,6)=cmu1*mutmx*rota
         v(n,7)=v(n,1)*v(n,6)/mutmx
        enddo
       enddo
      enddo
!
      return
      end
