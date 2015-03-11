module mod_vervol
implicit none
contains
      subroutine vervol(lm,vol)
!
!***********************************************************************
!
!     ACT
!_A    Detection de mailles a volumes negatifs.
!_A    En cas de volumes negatifs impression des informations correspondantes
!_A    et arret du calcul.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    imp        : com int              ; unite logiq, sorties de controle
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
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use sortiefichier 
implicit none
integer :: ind
integer :: i
integer :: j
integer :: k
integer :: lm
double precision :: vol
integer :: i1
integer :: i1m1
integer :: i2
integer :: img
integer :: j1
integer :: j1m1
integer :: j2
integer :: k1
integer :: k1m1
integer :: k2
integer :: kneg
integer :: l
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
!
!-----------------------------------------------------------------------
!
      dimension vol(ip11)
!
      ind(i,j,k)=n0+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
!
      n0=npc(lm)
      i1=ii1(lm)
      i2=ii2(lm)
      j1=jj1(lm)
      j2=jj2(lm)
      k1=kk1(lm)
      k2=kk2(lm)
!
      nid = id2(lm)-id1(lm)+1
      njd = jd2(lm)-jd1(lm)+1
      nijd = nid*njd
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      kneg=0
      l=mod(lm,lz)
      if (l.eq.0) l=lz
      img=(lm-l)/lz+1
      write(imp,1000) l,img
      do k=k1,k2-1
      do j=j1,j2-1
      do i=i1,i2-1
         n = ind(i,j,k)
!         write(imp,1001) i,j,n,vol(n)
!         if(vol(n).lt.0.) then
         if(vol(n).le.0.) then
         kneg=kneg+1
         write(imp,1001) i,j,k,vol(n)
         end if
       enddo
       enddo
       enddo
         write(imp,1002) kneg
         if(kneg.gt.0) then
         write(imp,1003)
         stop '  !! volumes negatifs !! '
         end if
!
 1000 format(//,5x,"volumes negatifs sur domaine no : ",i3/, &
             //,5x,"                 sur grille  no : ",i3/)
 1001 format(5x,' i=',i3,'        j=',i3,'        n=',i5,'        vol=',e12.4)
 1002 format(//,5x,' nombre de volumes negatifs :',i6)
 1003 format(///,5x,'arret du programme car mailles negatives')
!
      return
      end subroutine
end module
