module mod_met_smkor
  implicit none
contains
  subroutine met_smkor( &
       l,v, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       tprod,cfke, &
       qcxts5,qcyts6)
!
!***********************************************************************
!
!DA         avril 2002     - Eric Goncalves / SINUMEF
!_H
!     ACT
!_A   Calcul du second membre des equations pour k-omega
!_A   Modele de Kok avec condition realisabilite de Durbin
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
    integer          :: ind2,indc,   j,  j1,  j2
    integer          :: j2m1,   k,  k1,  k2,k2m1
    integer          ::    l,   m,   n, n0c, nid
    integer          :: nijd, njd
    double precision ::  alpha,    cd,  cfke,  dvxx,  dvxy
    double precision ::   dvxz,  dvyx,  dvyy,  dvyz,  dvzx
    double precision ::   dvzy,  dvzz,  gkgo,  omeg,qcxts5
    double precision :: qcyts6,  rcmu,    ss, tprod,txxf5x
    double precision :: txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z
    double precision ::      v,    xw
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension txxf5x(ip12),txyf5y(ip12),txzf5z(ip12), &
         tyyf6x(ip12),tyzf6y(ip12),tzzf6z(ip12), &
         qcxts5(ip12),qcyts6(ip12)
    dimension tprod(ip00)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00)
    dimension cfke(ip13)
!
    indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
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
    xw=beta/betas-sigmw*0.41**2/sqrt(betas)
    rcmu=1./sqrt(0.09)
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             gkgo=tyyf6x(n)*txxf5x(n)+tyzf6y(n)*txyf5y(n)+ &
                  tzzf6z(n)*txzf5z(n)
             omeg=v(n,7)/v(n,1)
             cd=sigmd*v(n,1)*max(gkgo,0.)/omeg
             ss=sqrt(4.*(dvxx(m)**2+dvyy(m)**2+dvzz(m)**2)/3. &
                  + (dvzy(m)+dvyz(m))**2 + (dvxz(m)+dvzx(m))**2 &
                  + (dvyx(m)+dvxy(m))**2)/omeg
             alpha=min(1.,rcmu/ss)
             qcxts5(n)=tprod(m)-betas*v(n,6)*omeg
             qcyts6(n)=xw*tprod(m)*v(n,7)/v(n,6)/alpha-beta*v(n,7)*omeg+cd
             cfke(n)=0.18*omeg
          enddo
       enddo
    enddo
!
    return
  end subroutine met_smkor
end module mod_met_smkor
