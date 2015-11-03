module mod_met_ke2mut
  implicit none
contains
  subroutine met_ke2mut( &
       l,ncyc, &
       v,mu,mut,dist,mnpar,ncin, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       frac)
!
!***********************************************************************
!
!     ACT
!_A   Calcul de mu_t avec la fonction correctrice de Jones Launder
!_A   ou Launder Sharma dans les regions de couches limites proches des
!_A   parois et la fonction de Smith pour les regions externes et les
!_A   sillages. Utilisation de la fonction de raccord de Menter.
!_A   Appel si equatt(1:4) egal a "2LS2" ou "2JL2"
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    txxf5x     : arg real(ip12     )  ; comp x grad(k)
!_I    txyf5y     : arg real(ip12     )  ; comp y grad(k)
!_I    txzf5z     : arg real(ip12     )  ; comp z grad(k)
!_I    tyyf6x     : arg real(ip12     )  ; comp x grad(e)
!_I    tyzf6y     : arg real(ip12     )  ; comp y grad(e)
!_I    tzzf6z     : arg real(ip12     )  ; comp z grad(e)
!
!     OUT
!_O    frac       : arg real(ip12     )  ; fonction de raccord des modeles
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use chainecarac
    implicit none
    integer          ::     i,   i1, i1m1,   i2, i2m1
    integer          ::    id,    j,   j1, j1m1,   j2
    integer          ::  j2m1,   jd,    k,   k1, k1m1
    integer          ::    k2, k2m1,   kd,    l,    m
    integer          :: mnpar,    n,  n0c,  nci, ncin
    integer          ::   ncj,  nck, ncyc,  nid, nijd
    integer          ::   njd
    double precision ::    c14,   c22,   c50,   chi,  chi2
    double precision ::   chi4,  cmb1,cmu500, coef1, coef2
    double precision ::  coef3, coef4,  dist, dist2,   eps
    double precision ::  exp2x,    f1,   fmu, fmujl, fmusm
    double precision ::   frac,    mu,   mut,  mut0, retur
    double precision ::    s24,txxf5x,txyf5y,txzf5z,tyyf6x
    double precision :: tyzf6y,tzzf6z,     v,    x1,    xk
    double precision ::     xl,  zeta
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension mut(ip12),mu(ip12),dist(ip12),mnpar(ip12), &
         txxf5x(ip12),txyf5y (ip12),txzf5z(ip12), &
         tyyf6x(ip12),tyzf6y (ip12),tzzf6z(ip12),frac(ip12)
    dimension ncin(ip41)
!
    cmu500=500.*cmukl
    c14=25.5**4
    c22=2.**2
    cmb1=cmukl* cklb1**(1./3.) /sqrt(2.)
    c50=50./xkappa**2
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
    nci  = inc(1,0,0)
    ncj  = inc(0,1,0)
    nck  = inc(0,0,1)
!     ----------------------------------------------------------
!com  calcul de la fonction de raccord entre modeles : F1 de Menter
!
!                       n : pointeur cellule tableaux tous domaines
!                       m : pointeur cellule tableaux un   domaine
!
!
    do k=k1,k2m1
       do j=j1,j2m1
          n=indc(i1m1,j,k)
          do i=i1,i2m1
             n=n+nci
             m=n-n0c
             xk   =v(n,6)/v(n,1)
             eps  =v(n,7)/v(n,1)
             coef1=sqrt(xk)*xk/(eps*dist(n))
             dist2=dist(n)**2
             coef2=cmu500*mu(n)*xk/(v(n,7)*dist2)
             coef4=txxf5x(n)*(tyyf6x(n)/eps-txxf5x(n)/xk)+ &
                  txyf5y(n)*(tyzf6y(n)/eps-txyf5y(n)/xk)+ &
                  txzf5z(n)*(tzzf6z(n)/eps-txzf5z(n)/xk)
             coef4=max(coef4,1.e-20)
             coef3=4.*xk/(dist2*coef4)
             zeta =min( max(coef1,coef2),coef3)
             if(zeta.le.2.5) then
                exp2x=exp(2.*zeta**4)
                frac(n)=(exp2x-1.)/(exp2x+1.)
             else
                frac(n)=1.
             end if
          enddo
       enddo
    enddo
!
!     --------------------------------------------------------------
!
    if(equatt(1:4).eq.'2JL2') then
       do k=k1,k2m1
          do j=j1,j2m1
             n=indc(i1-1,j,k)
             do i=i1,i2m1
                n=n+nci
!
!         fonction correctrice f_mu de Jones Launder
!
                retur=(v(n,6)**2)/(v(n,7)*mu(n))
                mut0=cmu*retur*mu(n)
                fmujl=exp(-2.5/(1.+retur/50.))
!
!         fonction correctrice f_mu de Smith
!
                chi=cmukl*retur
                chi2=chi**2
                chi4=chi2**2
                xk  =v(n,6)/v(n,1)
                eps =v(n,7)/v(n,1)
                xl  =cmb1*xk**1.5/eps
                f1  =exp(-c50*(xl/dist(n))**2)
                s24 =c22*chi2+chi4
                fmusm=((c14*f1+s24)/(c14+s24))**.25
                fmu=frac(n)*fmujl+(1.-frac(n))*fmusm
                mut(n)=fmu*mut0
             enddo
          enddo
       enddo
!
    elseif(equatt(1:3).eq.'2LS2') then
       do k=k1,k2m1
          do j=j1,j2m1
             n=indc(i1-1,j,k)
             do i=i1,i2m1
                n=n+nci
                retur=(v(n,6)**2)/(v(n,7)*mu(n))
                mut0=cmu*retur*mu(n)
                fmu=exp(-3.4/((1.+retur/50.)**2) )
                mut(n)=fmu*mut0
             enddo
          enddo
       enddo
    endif
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
  end subroutine met_ke2mut
end module mod_met_ke2mut
