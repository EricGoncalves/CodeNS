      subroutine met_smkesas( &
                 l, &
                 sn, &
                 vol,mu,s,cfke, &
                 tprod,bark,bare,tsv6,tsv7, &
                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                 t,dtdx,tn1,tn2,tn3, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT  calcul du terme source modèle k-e SAS      
!     Auteur: Jean Decaix 07/2010
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
      real mu,nlaplu,psas
      dimension s(ip11,ip60)
      dimension mu(ip12),tsv6(ip12),tsv7(ip12)
      dimension vol(ip11)
      dimension tprod(ip00),bark(ip00),bare(ip00)
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00), &
                dvyz(ip00),dvzx(ip00),dvzy(ip00),dvzz(ip00), &
                tn1(ip00),tn2(ip00),tn3(ip00),t(ip00),dtdx(ip00)
      dimension cfke(ip13)
      dimension sn(ip31*ndir)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
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
      i1p1=i1+1
      j1p1=j1+1
      k1p1=k1+1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
      i2p1=i2+1
      j2p1=j2+1
      k2p1=k2+1
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      nci  = inc(1,0,0)
      ncj  = inc(0,1,0)
      nck  = inc(0,0,1)
!

      npsn  =ndir*npfb(l)+1
      lgsnlt=nnn(l)
!     calcul du laplacien de la vitesse
!     premiere composante
       call met_laplaciens( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvxx,dvxy,dvxz, &
                 tn1,tn2,tn3, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!     deuxieme composante
        call met_laplaciens( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvyx,dvyy,dvyz, &
                 tn2,tn3,dtdx, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
       if(equat(3:4).eq.'3d') then
!        troisieme composante
          call met_laplaciens( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 s, &
                 dvzx,dvzy,dvzz, &
                 tn3,dtdx,t, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
        endif
!
!     constantes du modele     
      f1=1.
      qsi2=1.47   ! constante Menter
      xkappa=0.41    
      cc43=4./3.
!
      do k=k1,k2m1
       do j=j1,j2m1
        n=indc(i1m1,j,k)
        do i=i1,i2m1
          n=n+nci
          m=n-n0c
          retur=(s(n,6)**2)/(s(n,7)*mu(n))
          arg  =-retur**2
          f2   =(1.-0.3*exp(arg))
          esk  =s(n,7)/s(n,6)
          roe2sk=(s(n,7)*s(n,7))/s(n,6)
!         mise a zero des termes bas Reynolds pour lois de paroi
          if(lparoi.eq.1) then
            bark(m)=0.
            bare(m)=0.
            f2=1.
          endif
          ss=sqrt(cc43*(dvxx(m)**2+dvyy(m)**2+dvzz(m)**2) &
               + (dvzy(m)+dvyz(m))**2 &
               + (dvxz(m)+dvzx(m))**2 &
               + (dvyx(m)+dvxy(m))**2)
!         nlaplu=sqrt(tn1(m)**2+tn2(m)**2+tn3(m)**2)    !3D
          nlaplu=sqrt(tn1(m)**2+tn2(m)**2)
          xlvk=xkappa*ss/nlaplu
          xlvk2=xlvk*xlvk
          xl=((s(n,6)/s(n,1))**0.5)/(esk)
          xl2=xl*xl
          psas=xl2/xlvk2
          tsv6(n)=tprod(m)-s(n,7) + bark(m)
          tsv7(n)=cke1*f1*tprod(m)*esk+qsi2*psas*tprod(m)*esk-cke2*f2*roe2sk+bare(m)
!
!--------rayon spectral ke
!
          c1f1=cke1*f1
          c2f2=cke2*f2
          d=-bark(m)
          x=c1f1*tprod(m)-c2f2*s(n,7)
          g=amax1(0.,d*(d-2.*(x-c2f2*s(n,7))))
          rdelta=sqrt((x-c2f2*s(n,7))**2+4*x*s(n,7)+g)
          b=x-c2f2*s(n,7)-d
          xl1=0.5*(b+rdelta)/s(n,6)
          xl2=0.5*(b-rdelta)/s(n,6)
          cfke(n)=amax1(abs(xl1),abs(xl2))
         enddo
       enddo
      enddo
!
      return
      end
