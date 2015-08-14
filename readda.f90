module mod_readda
  implicit none
contains
  subroutine readda( &
       l,kda, &
       v,mut,utau, &
       vdual,vdual1,vdual2)
!
!***********************************************************************
!
!     ACT
!_A    Lecture des variables (v) et de la viscosite turbulente (mut)
!_A    en toute cellule (non fictive) d'un domaine structure.
!
!_I    l          : arg int              ; numero de domaine
!_I    kda        : arg int              ; unite logique, variables
!_I    equat      : arg char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
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
!     OUT
!_O    v          : arg real(ip11,ip60 ) ; variables de calcul
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    use chainecarac
    use modeleturb
    implicit none
    integer          ::        i,      i1,      i2,    i2m1,       j
    integer          ::       j1,      j2,    j2m1,       k,      k1
    integer          ::       k2,    k2m1,     kda,       l,       m
    integer          :: mdimtnxl,       n,      n0,     nid,    nijd
    integer          ::      njd,    resu,err
    double precision ::         mut(ip12),     v(ip11,ip60), vdual(ip11,ip60),vdual1(ip11,ip60)
    double precision :: vdual2(ip11,ip60)
    double precision,allocatable :: utau(:)
!
!-----------------------------------------------------------------------
!
!
!

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
    nijd = nid*njd
!
    read(kda)((( v(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    read(kda)((( v(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    read(kda)((( v(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    read(kda)((( v(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    read(kda)((( v(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    if(equat(1:2).eq.'ns') then   !calcul Navier-stokes
       read(kda)(((mut(ind(i,j,k)),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    endif
!
!
!cEG d  mise en commentaire du test pour reprise en multi-domaine
!       if(keinit.eq.0) then
!       pb pour relire un fichier multi-domaines ecri avec k-e mais
!       calcul initialise avec keinit=1. keinit ne devrait etre
!       utilise que pour initialiser le calcul de mut, pas pour gerer la lecture

    if(equat(6:7).eq.'ke') then  !calcul avec equations de transport
!       if((equat(6:7).eq.'ke' .and. keinit.eq.0) .or. &
!          (equat(2:4).eq.'res' .and. ip60.eq.7 )) then

!         reprise calcul avec equations de transport
       read(kda)((( v(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
       read(kda)((( v(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!
       if(kutau.eq.1) then
!         modeles de Chien ou k-omega bas Reynolds ou
!         k-omega de Wilcox ou Menter avec rugosite
          read(kda,iostat=resu)mdimtnxl
          deallocate(utau)
          ip42=mdimtnxl
          allocate(utau(mdimtnxl))
          if(resu .EQ. 0) then
             read(kda) utau(1:mdimtnxl)
          else
             do m=1,ip42
                utau(m)=1.e-04
             enddo
          endif
       endif

!        elseif(equat(6:7).eq.'ke' .and. keinit.ge.1) then    !nouveau calcul avec equations de transport
!          do k=k1,k2m1
!           do j=j1,j2m1
!            do i=i1,i2m1
!             v(ind(i,j,k),6)=0.
!             v(ind(i,j,k),7)=0.
!            enddo
!           enddo
!          enddo
    endif         !fin test keinit
!cEG d     mise en commentaire du test pour reprise en multi-domaine

!     call mpi_open(200,file='testlec.tec',form='formatted',status='unknown')
!     do j=1,j2m1
!      do i=1,i2m1
!       do k=1,1
!        write(1,'(8(e20.12),3i5)')  &
!           v(ind(i,j,k),1) ,v(ind(i,j,k),2),v(ind(i,j,k),3),v(ind(i,j,k),4),v(ind(i,j,k),5), &
!     v(ind(i,j,k),6),v(ind(i,j,k),7),mut(ind(i,j,k)),i,j,k
!        enddo
!       enddo
!     enddo
!     call mpi_close(200)

    if(kfmg.eq.3) then
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n = ind(i,j,k)
                vdual(n,1) = v(n,1)
                vdual(n,2) = v(n,2)
                vdual(n,3) = v(n,3)
                vdual(n,4) = v(n,4)
                vdual(n,5) = v(n,5)
                vdual(n,6) = v(n,6)
                vdual(n,7) = v(n,7)

                vdual1(n,1) = v(n,1)
                vdual1(n,2) = v(n,2)
                vdual1(n,3) = v(n,3)
                vdual1(n,4) = v(n,4)
                vdual1(n,5) = v(n,5)
                vdual1(n,6) = v(n,6)
                vdual1(n,7) = v(n,7)

                vdual2(n,1) = v(n,1)
                vdual2(n,2) = v(n,2)
                vdual2(n,3) = v(n,3)
                vdual2(n,4) = v(n,4)
                vdual2(n,5) = v(n,5)
                vdual2(n,6) = v(n,6)
                vdual2(n,7) = v(n,7)
             enddo
          enddo
       enddo

!   initialisation pour calcul ordre 3 en temps
!       write(*,*)'Readda'
!       call mpi_open(UNIT=98,FILE='facdual',FORM='unformatted',STATUS='unknown')
!       read(98) &
!           (((vdual(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       read(98) &
!           (((vdual(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!
!       READ(98) &
!           (((vdual1(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual1(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual1(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual1(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual1(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!            (((vdual1(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!            (((vdual1(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)

!       READ(98) &
!           (((vdual2(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual2(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual2(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual2(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!           (((vdual2(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!            (((vdual2(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       READ(98) &
!            (((vdual2(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!       call mpi_close(98)
    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine readda
end module mod_readda
