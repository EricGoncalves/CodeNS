module mod_met_smko
  implicit none
contains
  subroutine met_smko( &
       l,v, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       tprod,cfke, &
       qcxts5,qcyts6,cson)
!
!***********************************************************************
!
!_DA   DATE_C : mars 2010     - Eric Goncalves / LEGI
!
!     ACT
!_A   Calcul du second membre des equations pour k-omega
!_A   Modele compressible de Wilcox
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
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
!_L    gkgo       : arg real(ip00     )  ; grad(k).grad(omega)
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    implicit none
    integer          ::    i,  i1,  i2,i2m1,ind1
    integer          :: ind2,   j,  j1,  j2,j2m1
    integer          ::    k,  k1,  k2,k2m1,   l
    integer          ::    m,   n, n0c, nid,nijd
    integer          ::  njd
    double precision ::     aa, betac,betasc,    cd,  cfke
    double precision ::   cson,  gkgo,  omeg,qcxts5,qcyts6
    double precision ::  tprod,txxf5x,txyf5y,txzf5z,tyyf6x
    double precision :: tyzf6y,tzzf6z,     v,   xmt,  xmt0
    double precision ::     xw
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension cson(ip11)
    dimension txxf5x(ip12),txyf5y(ip12),txzf5z(ip12), &
         tyyf6x(ip12),tyzf6y(ip12),tzzf6z(ip12), &
         qcxts5(ip12),qcyts6(ip12)
    dimension tprod(ip00)
    dimension cfke(ip13)
!

!

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
!
!      xw=beta/betas-sigmw*0.41**2/sqrt(betas)
    xw=5./9.
!       xmt0=0.8  !coupure Mach turbulent
!       xmt0=0.5  !coupure Mach turbulent
!       xmt0=0.25  !coupure Mach turbulent
    xmt0=0.1  !coupure Mach turbulent
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
!         gkgo=tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+
!     &        tzzf6z(n)*txzf5z(n)
             omeg=v(n,7)/v(n,1)
!         cd=sigmd*v(n,1)*max(gkgo,0.)/omeg
             cd=0.
             xmt=sqrt(2.*v(n,6)/v(n,1))/cson(n)
             betasc=betas*(1.+1.5*fmt(xmt))
             betac=beta-betas*1.5*fmt(xmt)
             qcxts5(n)=tprod(m)-betasc*v(n,6)*omeg
             qcyts6(n)=xw*tprod(m)*v(n,7)/v(n,6)-betac*v(n,7)*omeg+cd
             cfke(n)=0.18*omeg   !=2*betas*omega
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
    function    fmt(aa)
      implicit none
      double precision ::  aa,fmt
      fmt=0.5*(1.+sign(1.,aa-xmt0))*(aa**2-xmt0**2)
    end function fmt
  end subroutine met_smko
end module mod_met_smko
