module mod_met_smdes
  implicit none
contains
  subroutine met_smdes( &
       l,x,y,z,delta, &
       sn, &
       vol,v,mu,dist, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       txxf5x,txyf5y,txzf5z,cfke, &
       t,dtdx,dtdy,dtdz,sdif, &
       qcxts5,qcyts6, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_H   DATE_C : juin 2002      - AUTEUR : Eric Goncalves / SINUMEF
!
!     ACT
!_A   Modele DES de Spalart.
!_A   Calcul du second  membre de l'equation pour nu-tilde.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    dist        : arg real(ip12      ) ; distance
!_L    dvxx       : arg real(ip00     )  ; grad(V)  vx,x
!_L    dvxy       : arg real(ip00     )  ; grad(V)  vx,y
!_L    dvxz       : arg real(ip00     )  ; grad(V)  vx,z
!_L    dvyx       : arg real(ip00     )  ; grad(V)  vy,x
!_L    dvyy       : arg real(ip00     )  ; grad(V)  vy,y
!_L    dvyz       : arg real(ip00     )  ; grad(V)  vy,z
!_L    dvzx       : arg real(ip00     )  ; grad(V)  vz,x
!_L    dvzy       : arg real(ip00     )  ; grad(V)  vz,y
!_L    dvzz       : arg real(ip00     )  ; grad(V)  vz,z
!_L    txxf5x     : arg real(ip12     )  ; comp x grad (v6)
!_L    txyf5y     : arg real(ip12     )  ; comp y grad (v6)
!_L    txzf5z     : arg real(ip12     )  ; comp z grad (v6)
!
!_/    qcxts5     : arg real(ip12     )  ; terme source equation pour nu-tilde
!_/    qcyts6     : arg real(ip12     )  ; terme source seconde equation=0
!-O
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  t,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  t,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  t,z
!_L    sdif       : arg real(ip00     )  ; grad(nu_tilde).grad(ro*nu_tulde)
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use chainecarac
    use mod_met_dist
    use mod_met_difsa
    implicit none
    integer          ::     i1,  i1m1,    i2,  i2m1
    integer          ::   imax,  imin,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,    j2,  j2m1,  jmax
    integer          ::   jmin,     k,    k1,  k1m1,    k2
    integer          ::   k2m1,  kmax,  kmin,     l,lgsnlt
    integer          ::      m,     n,   n0c,   nid,  nijd
    integer          ::    njd,  npsn
    double precision ::        cb2sig,   cfke(ip13),  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21)
    double precision ::   cmuj2(ip21),  cmuk1(ip21),  cmuk2(ip21),         ct42,        ctdes
    double precision ::          cv13,        cv133,         cw36,  delta(ip00),         dft2
    double precision ::          dfv1,         dfv2,          dfw,   dist(ip12),        dist2
    double precision ::           dpr,        dsdif,         dsm1,         dsm2,          dst
    double precision ::    dtdx(ip00),   dtdy(ip00),   dtdz(ip00),   dvxx(ip00),   dvxy(ip00)
    double precision ::    dvxz(ip00),   dvyx(ip00),   dvyy(ip00),   dvyz(ip00),   dvzx(ip00)
    double precision ::    dvzy(ip00),   dvzz(ip00),          dxg,          ft2,          fv1
    double precision ::           fv2,           fw,        kapd2,       kappa2,     mu(ip12)
    double precision ::       nutilde, qcxts5(ip12), qcyts6(ip12),        rtil6,       rtilde
    double precision ::    sdif(ip00),          sm1,          sm2,sn(ip31*ndir),       stilde
    double precision ::       t(ip00), txxf5x(ip12), txyf5y(ip12), txzf5z(ip12),         uns3
    double precision ::          uns6, v(ip11,ip60),    vol(ip11),         vort,      x(ip21)
    double precision ::            xg,          xg6,         xkhi,        xkhi2,        xkhi3
    double precision ::         xkhi4,      y(ip21),      z(ip21)
!
!-----------------------------------------------------------------------
!
!
!     ---------------------------------------------------------------
!com  sdif --> grad(ro nu_tilde).grad(nu_tilde) * cb2/sigma
!
    call met_difsa( &
         l, &
         sn, &
         vol,v, &
         t,dtdx,dtdy,dtdz, &
         txxf5x,txyf5y,txzf5z,sdif, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!     ---------------------------------------------------------------
!com  calcul de la dimension max d'une cellule de calcul
!
    call met_dist(l,x,y,z,delta)
!
!     --------------------------------------------------------------
!com  calcul des termes sources
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
    if(equat(3:5).eq.'2di') then
       imin=i1
       imax=i2m1
    endif
    if(equat(3:5).eq.'2dj') then
       jmin=j1
       jmax=j2m1
    endif
    if(equat(3:5).eq.'2dk') then
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
         t   , &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    kappa2=kappa**2
    uns6  =1./6.
    uns3  =-1./3.
    cv13=cv1**3
    cw36=cw3**6
!
    cv133 =3.*cv13
    ct42  =2.*ct4
    cb2sig=cb2/sigma
!
!     constante DES
      ctdes=0.65 !valeur init
!    ctdes=0.95
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             vort  =sqrt((dvzy(m)-dvyz(m))**2+(dvxz(m)-dvzx(m))**2 &
                  +(dvyx(m)-dvxy(m))**2)
!
             xkhi  =v(n,6)/mu(n)
             xkhi2 =xkhi**2
             xkhi3 =xkhi2*xkhi
!            dist2 =dist(n)**2
             dist2 =min(dist(n),ctdes*delta(m))**2
             nutilde=v(n,6)/v(n,1)
             kapd2  =kappa2*dist2
!
             fv1   =xkhi3/(xkhi3+cv13)
             fv2   =abs(1.-(xkhi/(1+xkhi*fv1)))
!            fv2   =(1.+xkhi/5.)**uns3
!            fv3 =(1.+xkhi*fv1)*(1.-fv2)/xkhi
             stilde=vort+nutilde*fv2/kapd2
!            stilde=fv3*vort+nutilde*fv2/kapd2
             rtilde=nutilde/(stilde*kapd2)
             rtilde=min(rtilde,1.)
             rtil6 =rtilde**6
             ft2   =ct3*exp(-ct4*xkhi2)
             sm1   =cb1*(1.-ft2)*stilde*v(n,6)
!
             xg   =rtilde+cw2*(rtil6-rtilde)
             xg6  =xg**6
             fw   =xg*((1.+cw36)/(xg6+cw36))**uns6
             sm2  =-(cw1*fw-cb1*ft2/kappa2)*v(n,6)*nutilde/dist2
!
             qcxts5(n)=sm1+sm2+sdif(m)
             qcyts6(n)=0.
!
!           rayon spectral matrice jacobienne terme source
!
             xkhi4=xkhi3*xkhi
             dfv1 =cv133*fv1**2/(mu(n)*xkhi4)
             dfv2 =(xkhi2*dfv1-1./mu(n))/((1.+xkhi*fv1)**2)
             dst  =(fv2/v(n,1)+nutilde*dfv2)/kapd2
             dpr  =(1./v(n,1)-nutilde*dst/stilde)/(stilde*kapd2)
             dxg  =dpr*(1.+cw2*(6.*rtil6/rtilde-1.))
             dfw  =(1.-xg6/(xg6+cw36))*fw*dxg/xg
             dft2 =-ct42*xkhi*ft2
             dsm1 =cb1*((1.-ft2)*(stilde+dst*v(n,6))- &
                  stilde*v(n,6)*dft2)
             dsm2 = 2.*sm2/v(n,6)-(cw1*dfw-cb1*dft2/kappa2)* &
                  nutilde*v(n,6) /dist2
             dsdif=-cb2sig*(txxf5x(m)*dtdx(m)+txyf5y(m)*dtdy(m)+ &
                  txzf5z(m)*dtdz(m))/(v(n,1)**2)
             cfke(n)=abs(dsm1+dsm2+dsdif)
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
  end subroutine met_smdes
end module mod_met_smdes
