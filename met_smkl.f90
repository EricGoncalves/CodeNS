module mod_met_smkl
  implicit none
contains
  subroutine met_smkl( &
       l, &
       sn, &
       vol,v,mu,mut,dist,mnpar, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       tprod,cfke, &
       t,dtdx,dtdy,dtdz,bark, &
       qcxts5,qcyts6, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A   Calcul du terme source des equations pour k-l
!_A   Modele de Smith
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    tprod      : arg real(ip00     )  ; production de k
!_I    dvxx       : arg real(ip00     )  ; grad(V)  vx,x
!_I    dvxy       : arg real(ip00     )  ; grad(V)  vx,y
!_I    dvxz       : arg real(ip00     )  ; grad(V)  vx,z
!_I    dvyx       : arg real(ip00     )  ; grad(V)  vy,x
!_I    dvyy       : arg real(ip00     )  ; grad(V)  vy,y
!_I    dvyz       : arg real(ip00     )  ; grad(V)  vy,z
!_I    dvzx       : arg real(ip00     )  ; grad(V)  vz,x
!_I    dvzy       : arg real(ip00     )  ; grad(V)  vz,y
!_I    dvzz       : arg real(ip00     )  ; grad(V)  vz,z
!
!     I/O
!_/    txxf5x     : arg real(ip12     )  ; comp x grad(k) puis
!_/                                        tenseur visqueux
!_/    txyf5y     : arg real(ip12     )  ; comp y grad(k) puis
!_/                                        tenseur visqueux
!_/    txzf5z     : arg real(ip12     )  ; comp z grad(k) puis
!_/                                        tenseur visqueux
!_/    tyyf6x     : arg real(ip12     )  ; comp x grad(e) puis
!_/                                        tenseur visqueux
!_/    tyzf6y     : arg real(ip12     )  ; comp y grad(e) puis
!_/                                        tenseur visqueux
!_/    tzzf6z     : arg real(ip12     )  ; comp z grad(e) puis
!_/                                        tenseur visqueux
!_/    qcxts5     : arg real(ip12     )  ; terme source equation pour k puis
!_/                                        vecteur flux de chaleur
!_/    qcyts6     : arg real(ip12     )  ; terme source seconde equation puis
!
!     LOC
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  t,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  t,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  t,z
!_L    bark       : arg real(ip00     )  ; terme faible reynolds equation k
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use chainecarac
    use mod_met_bark
    implicit none
    integer          ::           i,         i1,       i1m1,         i2,       i2m1
    integer          ::        imax,       imin,       ind1,       ind2
    integer          ::           j,         j1,       j1m1,         j2,       j2m1
    integer          ::        jmax,       jmin,          k,         k1
    integer          ::        k1m1,         k2,       k2m1,       kmax
    integer          ::        kmin,          l,     lgsnlt,          m,mnpar(ip12)
    integer          ::           n,        n0c,        nci,        nid,       nijd
    integer          ::         njd,       npsn
    double precision ::    bark(ip00),         c132,   cfke(ip13),  cmui1(ip21),  cmui2(ip21)
    double precision ::   cmuj1(ip21),  cmuj2(ip21),  cmuk1(ip21),  cmuk2(ip21),         csk2
    double precision ::          csl1,        csl12,   dist(ip12),         divv,        dqkdk
    double precision ::         dqkdl,        dqldk,        dqldl,       dssigl,   dtdx(ip00)
    double precision ::    dtdy(ip00),   dtdz(ip00),   dvxx(ip00),   dvxy(ip00),   dvxz(ip00)
    double precision ::    dvyx(ip00),   dvyy(ip00),   dvyz(ip00),   dvzx(ip00),   dvzy(ip00)
    double precision ::    dvzz(ip00),        gkgro,         glgk,         glgl,        glgro
    double precision ::      mu(ip12),    mut(ip12), qcxts5(ip12), qcyts6(ip12),         rack
    double precision ::         racrk,        racro,       rdelta,          sk2,          sl1
    double precision ::           sl2,          sl3,          sl4,sn(ip31*ndir),      t(ip00)
    double precision ::   tprod(ip00), txxf5x(ip12), txyf5y(ip12), txzf5z(ip12), tyyf6x(ip12)
    double precision ::  tyzf6y(ip12), tzzf6z(ip12), v(ip11,ip60),    vol(ip11),       xdelta
    double precision ::            xk,      xkapad2,           xl,          xl1,          xl2
    double precision ::       xlskap2
    logical          :: impli
!
!-----------------------------------------------------------------------
!
!
    impli=.true.
!      impli=.false.
!
!     ----------------------------------------------------------
!com  met_bark --> calcul de -2*mu*||grad(racine(k)||**2 (terme bas-reynolds)
!
    call met_bark( &
         l, &
         equat, &
         sn, &
         vol,v,mu, &
         t,dtdx,dtdy,dtdz,bark, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!     --------------------------------------------------------------
!     notation pour les constantes
!
!     B1 <-> cklb1
!     E2 <-> ckle2
!     sigma_k <-> sigmak
!     sigma_l <-> sigmal
!
    csk2  =2.*sqrt(2.)/cklb1
    csl1  =sqrt(2.)*(2.-ckle2)/cklb1
    dssigl=2./sigmal
!
    c132  =1.5*csk2
    csl12 =2.*csl1
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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    if(impli) then
!       calcul grad( rho ) pour implicitation
!       dtdx<->d(rho)/dx  dtdy<->d(rho)/dy  dtdz<->d(rho)/dz
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
!
       do k=kmin,kmax
          do j=jmin,jmax
             ind1=indc(imin,j,k)
             ind2=indc(imax,j,k)
             do n=ind1,ind2
                m=n-n0c
                t(m)=v(n,1)
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
            t, &
            dtdx,dtdy,dtdz, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
    endif
!
    nci=inc(1,0,0)
    do k=k1,k2m1
       do j=j1,j2m1
          n=indc(i1m1,j,k)
          do i=i1,i2m1
             n=n+nci
             m=n-n0c
             xk     =v(n,6)/v(n,1)
             xl     =v(n,7)/v(n,1)
             xkapad2=(xkappa*dist(n))**2
             xlskap2=xl*xl/xkapad2
             rack   =sqrt(xk)
             racrk  =sqrt(v(n,6))
             racro  =sqrt(v(n,1))
!           divv= div(V)
!           glgl= grad(l) . grad(l)
!           glgk= grad(l) . grad(k)/k
             divv =dvxx(m)+dvyy(m)+dvzz(m)
             glgl =tyyf6x(n)**2+tyzf6y(n)**2+tzzf6z(n)**2
             glgk =(tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+ &
                  tzzf6z(n)*txzf5z(n))/xk
!
             sk2  =csk2*rack*v(n,6)/xl
             sl1  =csl1*v(n,1)*rack*(1.-xlskap2)
             sl2  =v(n,7)*divv
             sl3  =mut(n)*glgl*xlskap2/(sigmal*xl)
             sl4  =dssigl*mut(n)*glgk
!
             qcxts5(n)=tprod(m) - sk2 + bark(m)
             qcyts6(n)=sl1 + sl2 - sl3 + sl4
!
!-----rayon spectral matrice jacobienne terme source pour implicitation
!        gkgro= ( grad(k)/k ) * ( grad(rho)/rho )
!        glgro= grad(l) . grad(rho)/rho
             gkgro=(txxf5x(n)*dtdx(m)+txyf5y(n)*dtdy(m)+ &
                  txzf5z(n)*dtdz(m))/v(n,6)
             glgro=(tyyf6x(n)*dtdx(m)+tyzf6y(n)*dtdy(m)+ &
                  tzzf6z(n)*dtdz(m))/v(n,1)
             dqkdk=-c132*racro*racrk/v(n,7)-bark(m)/v(n,6)+ &
                  mu (n)*gkgro/v(n,1)
             dqkdl=sk2/v(n,7)
             dqldl=-csl12*rack*xl/xkapad2+divv-mut(n)* &
                  (glgl-2*xl*glgro)/(sigmal*v(n,1)*xkapad2)
             dqldk=dssigl*mut(n)*(-glgk+glgro)/v(n,6)
             xdelta=dqkdk**2+dqldl**2+abs( -2.*dqkdk*dqldl+ &
                  4.*dqkdl*dqldk)
             rdelta=sqrt(xdelta)
             xl1=(dqkdk+dqldl+rdelta)*0.5
             xl2=(dqkdk+dqldl-rdelta)*0.5
             cfke(n)=max(abs(xl1),abs(xl2))
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
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_smkl
end module mod_met_smkl
