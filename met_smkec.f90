module mod_met_smkec
  implicit none
contains
  subroutine met_smkec( &
       l, &
       mu,s,cfke, &
       cson, &
       tprod,bark,bare,tsv6,tsv7)
!
!***********************************************************************
!
!   DATE  juillet 2010 - Jean DECAIX  - LEGI
!
!   ACT  
!   calcul du terme source pour modele k-e compressible type Sarkar
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    implicit none
    integer          ::    i,  i1,i1m1,i1p1,  i2
    integer          :: i2m1,  id,   j,  j1,j1p1
    integer          ::   j2,j2m1,  jd,   k,  k1
    integer          :: k1p1,  k2,k2m1,  kd,   l
    integer          ::    m,   n, n0c, nci, ncj
    integer          ::  nck, nid,nijd, njd
    double precision ::   alf1,  alf2,  alf3,   arg,     b
    double precision ::   bare,  bark,  c1f1,  c2f2,  cfke
    double precision ::   cson,     d,  epsc,   esk,    f1
    double precision ::     f2,     g,    mu,   pid,rdelta
    double precision ::  retur,roe2sk,     s, tprod,  tsv6
    double precision ::   tsv7,     x,   xl1,   xl2,   xmt
!
!-----------------------------------------------------------------------
!
    dimension s(ip11,ip60)
    dimension cson(ip11)
    dimension mu(ip12),tsv6(ip12),tsv7(ip12)
    dimension tprod(ip00),bark(ip00),bare(ip00)
    dimension cfke(ip13)
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
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
!
    nci  = inc(1,0,0)
    ncj  = inc(0,1,0)
    nck  = inc(0,0,1)
!
!     constantes du modele
    f1=1.
    alf2=0.15
    alf3=0.2  !valeur initiale
!      alf3=0.05
!      alf3=0.025
!      alf3=0.01
!      alf3=0.001
!      alf1=0.5  !valeur initiale
    alf1=0.1
!      alf1=0.05
!      alf1=0.025
!      alf1=0.005
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
!         calcul du Mach turbulent dans la cellule
             xmt=sqrt(2.*s(n,6)/s(n,1))/cson(n)
!         calcul du terme de correlation pression-dilatation
!          pid=-alf2*tprod(m)*xmt+alf3*s(n,7)*xmt**2
             pid=0.
!         calcul de la composante compressible de la dissipation
             epsc=alf1*s(n,7)*xmt**2
!          epsc=0.
!         mise a zero des termes bas Reynolds pour lois de paroi
             if(lparoi.eq.1) then
                bark(m)=0.
                bare(m)=0.
                f2=1.
             endif
             tsv6(n)=tprod(m)+pid-s(n,7)-epsc+bark(m)
             tsv7(n)=cke1*f1*tprod(m)*esk -cke2*f2*roe2sk+bare(m)
!
!---------rayon spectral ke
!
             c1f1=cke1*f1
             c2f2=cke2*f2
             d=-bark(m)
             x=c1f1*tprod(m)-c2f2*s(n,7)
             g=max(0.,d*(d-2.*(x-c2f2*s(n,7))))
             rdelta=sqrt((x-c2f2*s(n,7))**2+4*x*s(n,7)+g)
             b=x-c2f2*s(n,7)-d
             xl1=0.5*(b+rdelta)/s(n,6)
             xl2=0.5*(b-rdelta)/s(n,6)
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
  end subroutine met_smkec
end module mod_met_smkec
