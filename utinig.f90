module mod_utinig
  implicit none
contains
      subroutine utinig(l,x,y,z)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme utilisateur d' initialisation des coordonnees
!_A    aux points de discretisation.
!_A    Il doit remplir les tableaux x, y et z.
!_A
!_A     A partir d'un maillage monodomaine bi-dimensionnel, initialisation
!_A     d'une configuration bi-domaine par coupure en deux parties dans
!_A     la direction des indices i a la valeur i0 lue.
!_A     Dans le premier domaine deplacement de 20% d'une ligne k sur deux
!_A     jusqu'a la valeur kmil lue en donnee.
!_A     Rotation du domaine 2 par rapport au domaine 1 de l'angle protat.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    ptrans     : arg real             ; distance pour periodicite
!_I    protat     : arg real             ; angle(rad) pour periodicite
!_I    king       : arg int              ; cle initialisation maillage
!_I    lzx        : com int              ; nbr total de domaines
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
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
!_I    inig1      : com int              ; unite logiq, donnees utilisateur
!_I                                        initialisation geometrique 1
!
!     OUT
!_O    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_O    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_O    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!
!     LOC
!_L    tnte1      : arg real(ip11,ip60 ) ; tableau de travail
!
!***********************************************************************
!
      use para_var
      use para_fige
      use maillage
      use sortiefichier
    implicit none
    integer          ::    i,  i1,  i2,   j,  j1
    integer          ::   j2,   k,  k1,  k2, kdg
    integer          ::    l,   n, nid,nijd, njd
    integer          ::  npi, npj, npk, nkd
    integer          ::  n0 ,ni, nj, nk, numero
    double precision :: x,y,z
    logical          :: ecri
!     parameters locaux
!      parameter (idim=100,jdim=100,kdim=100)
!-----------------------------------------------------------------------
!
      character *80 titre
      character *20 fmt1
      dimension x(ip21),y(ip21),z(ip21)
!
      n0=npn(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
      ni=i2-i1+1
      nj=j2-j1+1
      nk=k2-k1+1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nkd=  kd2(l)-kd1(l)+1
      nijd = nid*njd
!
      if(l.eq.1) then
      close(inig1)
      open(inig1,file='finig1',form='formatted')
      end if
!
      print*,'finig1'
      read(inig1,2001) titre,numero
      print*,titre,numero
      read(inig1,2002) npi,npj,npk,fmt1
      print*,fmt1,npi,nid,npj,njd,npk,nkd,ip21

!!!!!!!!!!!!!!!!!
! 
!!!!!!!!!!!!!!!!!!!!!!
      if (numero.eq.1) then
      read(inig1,2003)(((x(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2003)(((y(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2003)(((z(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      elseif ((numero.ge.2).and.(numero.le.5)) then
      read(inig1,2004)(((x(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2004)(((y(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2004)(((z(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      else
      read(inig1,2005)(((x(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2005)(((y(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      read(inig1,2005)(((z(ind(i,j,k)),i=1,npi),j=1,npj),k=1,npk)
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maillage TIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 2003 format(977976ES16.9)
 2004 format(1361496ES16.9)
 2005 format(1030301ES16.9)
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maillage sphere canal_5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2003 format(977976ES16.9)
! 2004 format(1361496ES16.9)
! 2005 format(1030301ES16.9)
!!!!!!!!!!!!!!!!!!!!!!
! maillage sphere canal_5_plus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2003 format(1625151ES16.9)
! 2004 format(1625151ES16.9)
! 2005 format(1478741ES16.9)
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maillage sphere canal_plus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2003 format(1030301ES16.9)
! 2004 format(1030301ES16.9)
! 2005 format(1030301ES16.9)
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maillage sphere dans sphere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2003 format(977976ES16.9)
! 2004 format(1361496ES16.9)
! 2005 format(1030301ES16.9)
!!!!!!!!!!!!!!!!!!!!!!


!
!      open(unit=52,file='meshl.dat',form='formatted')
!      write(52,*) 'ZONE F=POINT, I=', npi, ',J=',npj  ,',K=',npk

!      print*,npi,npj,npk,l
!      do k=1,npk
!       do j=1,npj
!        do i=1,npi
!         write(52,*) x(ind(i,j,k)),y(ind(i,j,k)),z(ind(i,j,k))
!        enddo
!       enddo
!      enddo
!
 2001 format(a,/,i5)
 2002 format(3i5,a)

    return
    contains
    function    ind(i,j,k)
      implicit none
      integer          ::    i,ind,   j,   k
      ind=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    end subroutine utinig
    end module mod_utinig
