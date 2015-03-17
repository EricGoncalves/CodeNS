module mod_utreadav
  implicit none
contains
  subroutine utreadav( &
       l,inia1,equat, &
       u,v,mut,keinit)
!
!***********************************************************************
!
!     ACT
!_A    Lecture des variables (v) et de la viscosite turbulente (mut)
!_A    en toute cellule (non fictive) d'un domaine structure.
!
!     INP
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
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2m1
    integer          ::  inia1,     j,    j1,  j1m1,    j2
    integer          ::   j2m1,     k,    k1,  k1m1,    k2
    integer          ::   k2m1,keinit,     l,     n,    n0
    integer          ::    nci,   ncj,   nck,   nid,  nijd
    integer          ::    njd
    double precision :: mut,  u,  v
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension u(ip11,ip60),v(ip11,ip60),mut(ip12)
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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    nci=1
    ncj=nid
    nck=nijd
!
    if(l.eq.1) then
       print*,'INITIALISATION A PARTIR CHAMP C.V'
    endif
!
    read(inia1)((( u(ind(i,j,k),1),i=i1,i2),j=j1,j2),k=k1,k2)
    read(inia1)((( u(ind(i,j,k),2),i=i1,i2),j=j1,j2),k=k1,k2)
    read(inia1)((( u(ind(i,j,k),3),i=i1,i2),j=j1,j2),k=k1,k2)
    read(inia1)((( u(ind(i,j,k),4),i=i1,i2),j=j1,j2),k=k1,k2)
    read(inia1)((( u(ind(i,j,k),5),i=i1,i2),j=j1,j2),k=k1,k2)
!
!    Centre des cellules
!
    do k=k1,k2m1
       do j=j1,j2m1
          do i=i1,i2m1
             n=ind(i,j,k)
             v(n,1)=.125*(u(n        ,1)+u(n+        nck,1) &
                  +u(n+nci    ,1)+u(n+nci    +nck,1) &
                  +u(n    +ncj,1)+u(n    +ncj+nck,1) &
                  +u(n+nci+ncj,1)+u(n+nci+ncj+nck,1))
             v(n,2)=.125*(u(n        ,2)+u(n+        nck,2) &
                  +u(n+nci    ,2)+u(n+nci    +nck,2) &
                  +u(n    +ncj,2)+u(n    +ncj+nck,2) &
                  +u(n+nci+ncj,2)+u(n+nci+ncj+nck,2))
             v(n,3)=.125*(u(n        ,3)+u(n+        nck,3) &
                  +u(n+nci    ,3)+u(n+nci    +nck,3) &
                  +u(n    +ncj,3)+u(n    +ncj+nck,3) &
                  +u(n+nci+ncj,3)+u(n+nci+ncj+nck,3))
             v(n,4)=.125*(u(n        ,4)+u(n+        nck,4) &
                  +u(n+nci    ,4)+u(n+nci    +nck,4) &
                  +u(n    +ncj,4)+u(n    +ncj+nck,4) &
                  +u(n+nci+ncj,4)+u(n+nci+ncj+nck,4))
             v(n,5)=.125*(u(n        ,5)+u(n+        nck,5) &
                  +u(n+nci    ,5)+u(n+nci    +nck,5) &
                  +u(n    +ncj,5)+u(n    +ncj+nck,5) &
                  +u(n+nci+ncj,5)+u(n+nci+ncj+nck,5))
          enddo
       enddo
    enddo
!
    if(equat(1:2).eq.'ns') then
       read(inia1)((( u(ind(i,j,k),1),i=i1,i2),j=j1,j2),k=k1,k2)
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=ind(i,j,k)
                mut(n)=.125*(u(n        ,1)+u(n+        nck,1) &
                     +u(n+nci    ,1)+u(n+nci    +nck,1) &
                     +u(n    +ncj,1)+u(n    +ncj+nck,1) &
                     +u(n+nci+ncj,1)+u(n+nci+ncj+nck,1))
             enddo
          enddo
       enddo
    end if
!
    if(equat(6:7).eq.'ke'.and.keinit.eq.0) then
       read(inia1)((( u(ind(i,j,k),6),i=i1,i2),j=j1,j2),k=k1,k2)
       read(inia1)((( u(ind(i,j,k),7),i=i1,i2),j=j1,j2),k=k1,k2)
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=ind(i,j,k)
                v(n,6)=.125*(u(n        ,6)+u(n+        nck,6) &
                     +u(n+nci    ,6)+u(n+nci    +nck,6) &
                     +u(n    +ncj,6)+u(n    +ncj+nck,6) &
                     +u(n+nci+ncj,6)+u(n+nci+ncj+nck,6))
                v(n,7)=.125*(u(n        ,7)+u(n+        nck,7) &
                     +u(n+nci    ,7)+u(n+nci    +nck,7) &
                     +u(n    +ncj,7)+u(n    +ncj+nck,7) &
                     +u(n+nci+ncj,7)+u(n+nci+ncj+nck,7))
             enddo
          enddo
       enddo
!
    elseif(equat(6:7).eq.'ke') then
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=ind(i,j,k)
                v(n,6)=0.
                v(n,7)=0.
             enddo
          enddo
       enddo
    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine utreadav
end module mod_utreadav
