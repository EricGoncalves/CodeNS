module mod_metrics
  implicit none
contains
  subroutine metrics( &
       l,x,y,z,r,exs1,exs2, &
       sn,vol, &
       tn1,tn2,tn3)
!
!***********************************************************************
!
!     ACT
!_A    Calcul de la metrique pour un domaine structure:
!_A    distance a l'axe x(r), surfaces orientees des
!_A    facettes (sn) et volumes des cellules (vol).
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    exs1       : arg real             ; premier coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    exs2       : arg real             ; deuxieme coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
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
!_I    kvn        : com int              ; cle controle mailles deversees
!
!     OUT
!_O    r          : arg real(ip11      ) ; distance a l'axe
!_O    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!
!     LOC
!_L    tn1        : arg real(ip00      ) ; tableau de travail
!_L    tn2        : arg real(ip00      ) ; tableau de travail
!_L    tn3        : arg real(ip00      ) ; tableau de travail
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use sortiefichier
    use mod_vervol
    use mod_snorm
    use mod_svol
    use mod_extmhg
    implicit none
    integer          ::      i,    id,     j,    jd,     k
    integer          ::     kd,     l,lgsnlt,    nc,   nci
    integer          ::   ncij, ncijk,  ncik,   ncj,  ncjk
    integer          ::    nck,   nid,   njd,    nn,  npsn
    double precision ::          exs1,         exs2,      r(ip11),sn(ip31*ndir),    tn1(ip00)
    double precision ::     tn2(ip00),    tn3(ip00),    vol(ip11),      x(ip21),      y(ip21)
    double precision ::            yc,      z(ip21),           zc
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!



!
    nid=id2(l)-id1(l)+1
    njd=jd2(l)-jd1(l)+1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
    ncij = inc(1,1,0)
    ncik = inc(1,0,1)
    ncjk = inc(0,1,1)
    ncijk= inc(1,1,1)
!
    do k=kk1(l),kk2(l)-1
       do j=jj1(l),jj2(l)-1
!$OMP SIMD
          do i=ii1(l),ii2(l)-1
             nn=indn(i,j,k)
             nc=indc(i,j,k)
             yc=0.125*( y(nn     )+y(nn+nci  ) &
                  +y(nn+ncj )+y(nn+ncij ) &
                  +y(nn+nck )+y(nn+ncik ) &
                  +y(nn+ncjk)+y(nn+ncijk))
             zc=0.125*( z(nn     )+z(nn+nci  ) &
                  +z(nn+ncj )+z(nn+ncij ) &
                  +z(nn+nck )+z(nn+ncik ) &
                  +z(nn+ncjk)+z(nn+ncijk))
             r(nc)=sqrt(yc**2+zc**2)
          enddo
       enddo
    enddo
!
    call extmhg(l,x,y,z,exs1,exs2)
!
    npsn  =ndir*npfb(l)+1
    lgsnlt=nnn(l)
!
    call snorm( &
         l,x,y,z, &
         sn(npsn),lgsnlt)
    call svol( &
         l,x,y,z, &
         sn(npsn),lgsnlt,vol, &
         tn1,tn2,tn3)
    if(kvn.eq.1) then
       call vervol(l,vol)
    endif
!
!$OMP END MASTER
    return
  contains
    function    indn(i,j,k)
      implicit none
      integer          ::    i,indn,   j,   k
      indn=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nid*njd
    end function indn
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=npc(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nid*njd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nid*njd
    end function inc
  end subroutine metrics
end module mod_metrics
