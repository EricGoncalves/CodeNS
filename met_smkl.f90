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
implicit none
integer :: inc
integer :: indc
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: sn
double precision :: vol
double precision :: v
double precision :: dist
integer :: mnpar
double precision :: dvxx
double precision :: dvxy
double precision :: dvxz
double precision :: dvyx
double precision :: dvyy
double precision :: dvyz
double precision :: dvzx
double precision :: dvzy
double precision :: dvzz
double precision :: txxf5x
double precision :: txyf5y
double precision :: txzf5z
double precision :: tyyf6x
double precision :: tyzf6y
double precision :: tzzf6z
double precision :: tprod
double precision :: cfke
double precision :: t
double precision :: dtdx
double precision :: dtdy
double precision :: dtdz
double precision :: bark
double precision :: qcxts5
double precision :: qcyts6
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
double precision :: c132
double precision :: csk2
double precision :: csl1
double precision :: csl12
double precision :: divv
double precision :: dqkdk
double precision :: dqkdl
double precision :: dqldk
double precision :: dqldl
double precision :: dssigl
double precision :: gkgro
double precision :: glgk
double precision :: glgl
double precision :: glgro
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: imax
integer :: imin
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: jmax
integer :: jmin
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: kmax
integer :: kmin
integer :: lgsnlt
integer :: m
integer :: n
integer :: n0c
integer :: nci
integer :: nid
integer :: nijd
integer :: njd
integer :: npsn
double precision :: rack
double precision :: racrk
double precision :: racro
double precision :: rdelta
double precision :: sk2
double precision :: sl1
double precision :: sl2
double precision :: sl3
double precision :: sl4
double precision :: xdelta
double precision :: xk
double precision :: xkapad2
double precision :: xl
double precision :: xl1
double precision :: xl2
double precision :: xlskap2
!
!-----------------------------------------------------------------------
!
      logical impli
      real mu,mut
      dimension v(ip11,ip60)
      dimension mut(ip12),mu(ip12),dist(ip12),mnpar(ip12), &
                txxf5x(ip12),txyf5y(ip12),txzf5z(ip12), &
                tyyf6x(ip12),tyzf6y(ip12),tzzf6z(ip12), &
                qcxts5(ip12),qcyts6(ip12)
      dimension cfke(ip13),vol(ip11)
      dimension sn(ip31*ndir)              
      dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
                dvyx(ip00),dvyy(ip00),dvyz(ip00), &
                dvzx(ip00),dvzy(ip00),dvzz(ip00), &
                dtdx(ip00),dtdy(ip00),dtdz(ip00), &
                t(ip00),tprod(ip00),bark(ip00)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd
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
      end
