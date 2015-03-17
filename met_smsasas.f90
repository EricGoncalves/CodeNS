module mod_met_smsasas
  implicit none
contains
  subroutine met_smsasas( &
       l, &
       sn, &
       vol,v,mu,dist, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       txxf5x,txyf5y,txzf5z,cfke, &
       t,dtdx,dtdy,dtdz,sdif, &
       qcxts5,qcyts6, &
       tn1,tn2,tn3,tn10, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     DATE 
!      juillet 2010 - Jean DECAIX
!     ACT
!_A    Modele de Spalart-Allmaras avec terme SAS.
!_A    Calcul du second  membre de l'equation pour nu-tilde.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    dist       : arg real(ip12      ) ; distance
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
    use mod_met_laplaciens
    use mod_met_difsa
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2m1
    integer          ::   imax,  imin,  ind1,  ind2,     j
    integer          ::     j1,  j1m1,    j2,  j2m1,  jmax
    integer          ::   jmin,     k,    k1,  k1m1,    k2
    integer          ::   k2m1,  kmax,  kmin,     l,lgsnlt
    integer          ::      m,     n,   n0c,   nid,  nijd
    integer          ::    njd,  npsn
    double precision ::  cb2sig,   cc43,   cfke,  cmui1,  cmui2
    double precision ::   cmuj1,  cmuj2,  cmuk1,  cmuk2,   csas
    double precision ::    ct42,   cv13,  cv133,   cw36,   dft2
    double precision ::    dfv1,   dfv2,    dfw,   dist,  dist2
    double precision ::     dpr,  dsdif,   dsm1,   dsm2,    dst
    double precision ::    dtdx,   dtdy,   dtdz,   dvxx,   dvxy
    double precision ::    dvxz,   dvyx,   dvyy,   dvyz,   dvzx
    double precision ::    dvzy,   dvzz,    dxg,     fd,    ft2
    double precision ::     fv1,    fv2,     fw,  kapd2, kappa2
    double precision ::      mu, nlaplu,nutilde, qcxts5, qcyts6
    double precision ::   rtil3,  rtil6, rtilde,   sdif,    sm1
    double precision ::     sm2,     sn,     ss, stilde,      t
    double precision ::     tn1,   tn10,    tn2,    tn3, txxf5x
    double precision ::  txyf5y, txzf5z,   uns6,      v,    vol
    double precision ::    vort,     xg,    xg6,   xkhi,  xkhi2
    double precision ::   xkhi3,  xkhi4,   xlvk,  xlvk2
!
!-----------------------------------------------------------------------
!
!
    dimension v(ip11,ip60)
    dimension mu(ip12),dist(ip12), &
         txxf5x(ip12),txyf5y(ip12),txzf5z(ip12), &
         qcxts5(ip12),qcyts6(ip12)
    dimension cfke(ip13)
    dimension sn(ip31*ndir)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00), &
         dtdx(ip00),dtdy(ip00),dtdz(ip00), &
         t(ip00),tn1(ip00),tn2(ip00),tn3(ip00),tn10(ip00)
    dimension sdif(ip00),vol(ip11)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
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
!     calcul grad( rho ) pour implicitation
!     dtdx<->d(rho)/dx  dtdy<->d(rho)/dy  dtdz<->d(rho)/dz
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
         t, &
         dtdx,dtdy,dtdz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!      calcul du laplacien de la vitesse
!       première composante
    call met_laplaciens( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         v, &
         dvxx,dvxy,dvxz, &
         tn1,tn2,tn3, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!       deuxième composante
    call met_laplaciens( &
         l, &
         equat, &
         sn(npsn),lgsnlt, &
         vol, &
         v, &
         dvyx,dvyy,dvyz, &
         tn2,tn3,tn10, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!       troisième composante
    if(equat(3:4).eq.'3d') then
       call met_laplaciens( &
            l, &
            equat, &
            sn(npsn),lgsnlt, &
            vol, &
            v, &
            dvzx,dvzy,dvzz, &
            tn3,tn10,t, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
    endif
!
!     constantes du modele de Spalart
    kappa2=kappa**2
    uns6=1./6.
    cv13=cv1**3
    cw36=cw3**6
    cv133=3.*cv13
    ct42=2.*ct4
    cb2sig=cb2/sigma
    cc43=4./3. 
!
    csas=1.   ! constante modèle SAS
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1  ,j,k)
          ind2=indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             vort=sqrt((dvzy(m)-dvyz(m))**2+(dvxz(m)-dvzx(m))**2 &
                  +(dvyx(m)-dvxy(m))**2)
!
             xkhi=v(n,6)/mu(n)
             xkhi2=xkhi**2
             xkhi3=xkhi2*xkhi
             dist2=dist(n)**2
             nutilde=v(n,6)/v(n,1)
             kapd2=kappa2*dist2
!
             fv1=xkhi3/(xkhi3+cv13)
!         fv2=1.-(xkhi/(1+xkhi*fv1))
             fv2=abs(1.-(xkhi/(1.+xkhi*fv1)))
             stilde=vort+nutilde*fv2/kapd2
             rtilde=nutilde/(stilde*kapd2)
!          rtilde=min(rtilde,1.)
             rtil6=rtilde**6
             rtil3=(8.*rtilde)**3          
             fd=1.-tanh(rtil3)            ! fonction fd DDES
             ft2=ct3*exp(-ct4*xkhi2)
             sm1=cb1*(1.-ft2)*stilde*v(n,6)
             xg=rtilde+cw2*(rtil6-rtilde)
             xg6=xg**6
             fw=xg*((1.+cw36)/(xg6+cw36))**uns6
             sm2=-(cw1*fw-cb1*ft2/kappa2)*v(n,6)*nutilde/dist2   !sans SAS

             if(fd.gt.0.8)then    ! activation terme SAS avec fonction fd
                ss=sqrt(cc43*(dvxx(m)**2+dvyy(m)**2+dvzz(m)**2) &
                     + (dvzy(m)+dvyz(m))**2 &
                     + (dvxz(m)+dvzx(m))**2+(dvyx(m)+dvxy(m))**2)
!           nlaplu=sqrt(tn1(m)**2+tn2(m)**2+tn3(m)**2) !3D
                nlaplu=sqrt(tn1(m)**2+tn2(m)**2)
                xlvk=kappa*ss/nlaplu
                xlvk2=xlvk**2
                sm2=-cw1*fw*v(n,6)*csas*nutilde/xlvk2 &
                     + cb1*ft2*v(n,6)*nutilde/(kappa2*dist2)
             endif
             qcxts5(n)=sm1+sm2+sdif(m)
             qcyts6(n)=0.

!         rayon spectral matrice jacobienne terme source
             xkhi4=xkhi3*xkhi
             dfv1=cv133*fv1**2/(mu(n)*xkhi4)
             dfv2=(xkhi2*dfv1-1./mu(n))/((1.+xkhi*fv1)**2)
             dst=(fv2/v(n,1)+nutilde*dfv2)/kapd2
             dpr=(1./v(n,1)-nutilde*dst/stilde)/(stilde*kapd2)
             dxg=dpr*(1.+cw2*(6.*rtil6/rtilde-1.))
             dfw=(1.-xg6/(xg6+cw36))*fw*dxg/xg
             dft2=-ct42*xkhi*ft2
             dsm1=cb1*((1.-ft2)*(stilde+dst*v(n,6))- &
                  stilde*v(n,6)*dft2)
             dsm2=2.*sm2/v(n,6)-(cw1*dfw-cb1*dft2/kappa2)* &
                  nutilde*v(n,6)/dist2
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
  end subroutine met_smsasas
end module mod_met_smsasas
