module mod_zvismo
  implicit none
contains
  subroutine zvismo(l,mu,s,temp)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 2004 - Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul du coefficient de viscosite moleculaire.
!
!
!_I    l          : arg int              ; numero de domaine
!_I    s          : arg real(ip11,ip60 ) ; variables de calcul
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    tnz        : com real             ; etat pour adimensionnement,
!_I                                        temperature
!_I    reynz      : com real             ; nombre de Reynolds calcule avec
!_I                                        les grandeurs d'adimensionnement,
!_I                                        pour definir la loi de Sutherland
!
!     OUT
!_O    mu         : arg real(ip12      ) ; viscosite moleculaire
!
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use definition
    use proprieteflu
    use sortiefichier
    implicit none
    integer          ::     i1,    i2,  i2m1,iarret
    integer          ::   ind1,  ind2,     j,    j1,    j2
    integer          ::   j2m1,     k,    k1,    k2,  k2m1
    integer          ::      l,     n,    n0,   nid,  nijd
    integer          ::    njd
    double precision ::   mu(ip12),s(ip11,ip60),temp(ip11)
    double precision ::   a, bl, usrey
!
!-----------------------------------------------------------------------
!
    iarret=0
!
    n0=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd

    if(iflu.eq.1) then    !air  
!   constante sutherland pour l'air (=110.4K)
      a=110.4/tnz
    elseif(iflu.eq.2) then  !eau froide
!     constante loi viscosite liquide: mu=A*exp(B/T)
!     pour l'eau  : A=1.24e-6 Pa.s   et B=1968K
!     pour le R114: A=10.336e-6 Pa.s et B=976.9738K
       bl=1968./tnz
    elseif(iflu.eq.3) then  !freon R114
       bl=976.9738/tnz
    endif
    usrey=1./reynz

    if(iflu.eq.1) then
     do k=k1,k2m1
      do j=j1,j2m1
       ind1=ind(i1  ,j,k)
       ind2=ind(i2m1,j,k)
       do n=ind1,ind2
         if(temp(n).le.0.) then
           iarret=iarret+1
           mu(n)=0.
         else
           mu(n)=usrey*temp(n)*sqrt(temp(n))*(1.+a)/(temp(n)+a)
         endif
       enddo
      enddo
     enddo
     elseif(iflu.ge.2) then
      do k=k1,k2m1
       do j=j1,j2m1
        ind1=ind(i1  ,j,k)
        ind2=ind(i2m1,j,k)
        do n=ind1,ind2
         if(temp(n).le.0.) then
          iarret=iarret+1
          mu(n)=0.
         else
!          mu(n)=usrey*exp(bl*(1./temp(n)-1.))
          mu(n)=usrey
         endif
        enddo
       enddo
      enddo
     endif
!
    if(iarret.ne.0) then
       write(imp,'(/,"!!!zvismo: temperature negative en",i8," cellules domaine= ",i4)')iarret,l
       call abort
       stop
    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine zvismo
end module mod_zvismo
