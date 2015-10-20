module mod_met_dist
  implicit none
contains
  subroutine met_dist(l,x,y,z,delta)
!
!***********************************************************************
!
!_DA  DATE_C : juin 2002 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul de la dimension max d'une cellule de calcul par direction.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!
!_OUT Delta       : arg real(ip00      ) ; max(dx,dy,dz)
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    implicit none
    integer          ::     i,   i1,   i2, i2m1
    integer          ::     j,   j1,   j2, j2m1
    integer          ::     k,   k1,   k2, k2m1
    integer          ::     l,    m,    n,  n0c,  nci
    integer          ::  ncij,ncijk, ncik,  ncj, ncjk
    integer          ::   nck,  nid,  njd
    double precision :: delta(ip00),        dx1,        dx2,        dx3,        dx4
    double precision ::         dy1,        dy2,        dy3,        dy4,        dz1
    double precision ::         dz2,        dz3,        dz4,    x(ip21),    y(ip21)
    double precision ::     z(ip21)
!
!-----------------------------------------------------------------------
!
!


!
    nid=id2(l)-id1(l)+1
    njd=jd2(l)-jd1(l)+1
!
    nci  = inc(1,0,0)
    ncj  = inc(0,1,0)
    nck  = inc(0,0,1)
    ncij = inc(1,1,0)
    ncik = inc(1,0,1)
    ncjk = inc(0,1,1)
    ncijk= inc(1,1,1)
!
    n0c=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
!
!******************************************************************************
    if(equat(3:5).eq.'2dk') then
!******************************************************************************
       do j=j1,j2m1
          do i=i1,i2m1
             n=indn(i,j,k1)
             m=n-n0c
             dx1=sqrt((x(n)-x(n+nci))**2+(y(n)-y(n+nci))**2+ &
                  (z(n)-z(n+nci))**2)
             dx2=sqrt((x(n+ncj)-x(n+ncij))**2+(y(n+ncj)-y(n+ncij))**2+ &
                  (z(n+ncj)-z(n+ncij))**2)
!
             dy1=sqrt((x(n)-x(n+ncj))**2+(y(n)-y(n+ncj))**2+ &
                  (z(n)-z(n+ncj))**2)
             dy2=sqrt((x(n+nci)-x(n+ncij))**2+(y(n+nci)-y(n+ncij))**2+ &
                  (z(n+nci)-z(n+ncij))**2)
!
             Delta(m)=max(dx1,dx2,dy1,dy2)
          enddo
       enddo
!
!******************************************************************************
    elseif(equat(3:4).eq.'3d') then
!******************************************************************************
!
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=indn(i,j,k)
                m=n-n0c
                dx1=sqrt((x(n)-x(n+nci))**2+(y(n)-y(n+nci))**2+ &
                     (z(n)-z(n+nci))**2)
                dx2=sqrt((x(n+ncj)-x(n+ncij))**2+(y(n+ncj)-y(n+ncij))**2+ &
                     (z(n+ncj)-z(n+ncij))**2)
                dx3=sqrt((x(n+nck)-x(n+ncik))**2+(y(n+nck)-y(n+ncik))**2+ &
                     (z(n+nck)-z(n+ncik))**2)
                dx4=sqrt((x(n+ncjk)-x(n+ncijk))**2+(y(n+ncjk)-y(n+ncijk))**2+ &
                     (z(n+ncjk)-z(n+ncijk))**2)
!
                dy1=sqrt((x(n)-x(n+ncj))**2+(y(n)-y(n+ncj))**2+ &
                     (z(n)-z(n+ncj))**2)
                dy2=sqrt((x(n+nci)-x(n+ncij))**2+(y(n+nci)-y(n+ncij))**2+ &
                     (z(n+nci)-z(n+ncij))**2)
                dy3=sqrt((x(n+nck)-x(n+ncjk))**2+(y(n+nck)-y(n+ncjk))**2+ &
                     (z(n+nck)-z(n+ncjk))**2)
                dy4=sqrt((x(n+ncik)-x(n+ncijk))**2+(y(n+ncik)-y(n+ncijk))**2+ &
                     (z(n+ncik)-z(n+ncijk))**2)
!
                dz1=sqrt((x(n)-x(n+nck))**2+(y(n)-y(n+nck))**2+ &
                     (z(n)-z(n+nck))**2)
                dz2=sqrt((x(n+nci)-x(n+ncik))**2+(y(n+nci)-y(n+ncik))**2+ &
                     (z(n+nci)-z(n+ncik))**2)
                dz3=sqrt((x(n+ncj)-x(n+ncjk))**2+(y(n+ncj)-y(n+ncjk))**2+ &
                     (z(n+ncj)-z(n+ncjk))**2)
                dz4=sqrt((x(n+ncij)-x(n+ncijk))**2+(y(n+ncij)-y(n+ncijk))**2+ &
                     (z(n+ncij)-z(n+ncijk))**2)
!
                Delta(m)=max(dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4,dz1,dz2,dz3,dz4)
             enddo
          enddo
       enddo
!******************************************************************
    endif
!
    return
  contains
    function    indn(i,j,k)
      implicit none
      integer          ::    i,indn,   j,   k
      indn=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nid*njd
    end function indn
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nid*njd
    end function inc
  end subroutine met_dist
end module mod_met_dist
