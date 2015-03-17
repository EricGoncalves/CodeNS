module mod_met_smmtr
  implicit none
contains
  subroutine met_smmtr( &
       l, &
       v,mu,mut,dist,mnpar,ncin, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       tprod,cfke,frac, &
       qcxts5,qcyts6)
!
!***********************************************************************
!
!_H   DATE_C : septembre 2002  - AUTEUR: Eric Goncalves / SINUMEF
!
!     ACT
!_A   Calcul du second membre des equations pour k-omega de Menter.
!_A   La fonction de raccord des modeles est aussi calculee.
!_A   Conditions de realisabilite de Durbin.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant courant
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_L    tprod      : arg real(ip00     )  ; production de k
!
!_O    frac       : arg real(ip12     )  ; fonction de raccord des modeles
!
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
!_L    gkgo       ; grad(k) * grad (omega)
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use mod_met_mtcorf1
    implicit none
    integer          ::     i,   i1, i1m1,   i2, i2m1
    integer          ::  ind1, ind2,    j,   j1, j1m1
    integer          ::    j2, j2m1,    k,   k1, k1m1
    integer          ::    k2, k2m1,    l,    m,mnpar
    integer          ::     n,  n0c, ncin,  nid, nijd
    integer          ::   njd
    double precision ::  alpha,  cfke, coef1, coef2, coef3
    double precision ::  coef4,  dist, dist2, dkomg,  dvxx
    double precision ::   dvxy,  dvxz,  dvyx,  dvyy,  dvyz
    double precision ::   dvzx,  dvzy,  dvzz, exp2x,   fm1
    double precision ::    fm2,  frac,  gkgo,    mu,   mut
    double precision ::   omeg,qcxts5,qcyts6,  rcmu,   sl3
    double precision ::     ss, tprod,txxf5x,txyf5y,txzf5z
    double precision :: tyyf6x,tyzf6y,tzzf6z,     v,w4sig2
    double precision ::   wgam,  wsig, xgam1, xgam2,  zeta
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension mut(ip12),mu(ip12),dist(ip12),mnpar(ip12), &
         txxf5x(ip12),txyf5y (ip12),txzf5z(ip12), &
         tyyf6x(ip12),tyzf6y (ip12),tzzf6z(ip12), &
         qcxts5(ip12),qcyts6 (ip12),frac(ip12)
    dimension tprod(ip00)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00)
    dimension ncin(ip41),cfke(ip13)
!

!
    xgam1=beta1/betae-sigma1*okappa**2/sqrt(betae)
    xgam2=beta2/betae-sigma2*okappa**2/sqrt(betae)
    w4sig2=4.*wsig2
    rcmu=1./sqrt(0.09)
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
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
! ----------------------------------------------------------
!com  calcul de la fonction de raccord entre modeles
! ----------------------------------------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             omeg=v(n,7)/v(n,1)
             gkgo=tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+ &
                  tzzf6z(n)*txzf5z(n)
             coef1=sqrt(v(n,6)/v(n,1))/(betae*omeg*dist(n))
             dist2=dist(n)**2
             coef2=500.*mu(n)/(v(n,7)*dist2)
             coef4=wsig2*v(n,1)*gkgo/omeg
             dkomg=max(coef4,1.e-20)
             coef3=w4sig2*v(n,6)/(dist2*dkomg)
             zeta =min( max(coef1,coef2),coef3)
             if(zeta.le.2.5) then
                exp2x=exp(2.*zeta**4)
                frac(n)=(exp2x-1.)/(exp2x+1.)
             else
                frac(n)=1.
             endif
          enddo
       enddo
    enddo
!
! --------------------------------------------------------------------------
! verification que la fonction de raccord reste a 1 dans la region de paroi
! --------------------------------------------------------------------------
!
    call met_mtcorf1( &
         l,ncin, &
         dist,mnpar,frac)
!
!  -----------------------------------------------------------------------
!com  calcul des termes sources et du rayon spectral pour implicitation
! ------------------------------------------------------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             fm1=frac(n)
             fm2=1.-fm1
             wgam=fm1*xgam1+fm2*xgam2
             wsig=fm1*wsig1+fm2*wsig2
             beta=fm1*beta1+fm2*beta2
             gkgo=tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+ &
                  tzzf6z(n)*txzf5z(n)
             omeg=v(n,7)/v(n,1)
             sl3=2.*wsig*v(n,1)*gkgo/omeg
             ss=sqrt(4.*(dvxx(m)**2+dvyy(m)**2+dvzz(m)**2)/3. &
                  + (dvzy(m)+dvyz(m))**2 + (dvxz(m)+dvzx(m))**2 &
                  + (dvyx(m)+dvxy(m))**2)/omeg
             alpha=min(1.,rcmu/ss)
             qcxts5(n)=tprod(m)- betae*v(n,6)*omeg
             qcyts6(n)=tprod(m)*wgam*v(n,1)/mut(n)/alpha &
                  - beta*v(n,7)*omeg + sl3
             cfke(n)=0.18*omeg
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
  end subroutine met_smmtr
end module mod_met_smmtr
