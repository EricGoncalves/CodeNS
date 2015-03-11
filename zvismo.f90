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
integer :: ind
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: s
double precision :: temp
integer :: i1
integer :: i2
integer :: i2m1
integer :: iarret
integer :: ind1
integer :: ind2
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
double precision :: usrey
!
!-----------------------------------------------------------------------
!
      real mu,bl,a
      dimension mu(ip12),temp(ip11)
      dimension s(ip11,ip60)
!
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
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

!     constante sutherland pour l'air (=110.4K)
      a=110.4/tnz
      usrey=1./reynz

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
!
      if(iarret.ne.0) then
        write(imp,'(/,"!!!zvismo: temperature negative en",i8," cellules domaine= ",i4)')iarret,l
          stop
      endif
!
      return
      end subroutine
end module
