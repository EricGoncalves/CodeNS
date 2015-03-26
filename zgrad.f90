module mod_zgrad
  implicit none
contains
  subroutine zgrad( &
       l, &
       equat, &
       sn,lgsnlt, &
       vol, &
       s,temp, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dtx,dty,dtz)
!
!***********************************************************************
!
!     ACT
!_A    Calculs des gradients de vitesse et de temperature.
!
!
!_I    l          : arg int              ; numero de domaine
!_I    equat      : arg char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    sn         : arg real(lgsnlt,
!_I                          nind,ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    vx         : arg real(ip00      ) ; composante en x de la vitesse
!_I    vy         : arg real(ip00      ) ; composante en y de la vitesse
!_I    vz         : arg real(ip00      ) ; composante en z de la vitesse
!_I    temp       : arg real(ip00      ) ; temperature
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
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!
!     OUT
!_O    dvxx       : arg real(ip00      ) ; composante xx du gradient de vitesse
!_O    dvxy       : arg real(ip00      ) ; composante xy du gradient de vitesse
!_O    dvxz       : arg real(ip00      ) ; composante xz du gradient de vitesse
!_O    dvyx       : arg real(ip00      ) ; composante yx du gradient de vitesse
!_O    dvyy       : arg real(ip00      ) ; composante yy du gradient de vitesse
!_O    dvyz       : arg real(ip00      ) ; composante yz du gradient de vitesse
!_O    dvzx       : arg real(ip00      ) ; composante zx du gradient de vitesse
!_O    dvzy       : arg real(ip00      ) ; composante zy du gradient de vitesse
!_O    dvzz       : arg real(ip00      ) ; composante zz du gradient de vitesse
!_O    dtx        : arg real(ip00      ) ; composante x du gradient de temp
!_O    dty        : arg real(ip00      ) ; composante y du gradient de temp
!_O    dtz        : arg real(ip00      ) ; composante z du gradient de temp
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
  integer          ::      i,    i1,  i1m1,  i1p1,    i2
  integer          ::   i2m1,    id,  imax,  imin,  ind1
  integer          ::   ind2,     j,    j1,  j1m1,  j1p1
  integer          ::     j2,  j2m1,    jd,  jmax,  jmin
  integer          ::      k,    k1,  k1m1,  k1p1,    k2
  integer          ::   k2m1,    kd,  kmax,  kmin,     l
  integer          :: lgsnlt,     m,     n,    n0,   nci
  integer          ::    ncj,   nck,   nid,  nijd,   njd
  double precision ::                   c0,           dtx(ip00),           dty(ip00),           dtz(ip00),          dvxx(ip00)
  double precision ::           dvxy(ip00),          dvxz(ip00),          dvyx(ip00),          dvyy(ip00),          dvyz(ip00)
  double precision ::           dvzx(ip00),          dvzy(ip00),          dvzz(ip00),                 eps,        s(ip11,ip60)
  double precision :: sn(lgsnlt,nind,ndir),          temp(ip11),                  ts,           vol(ip11),                vols
  double precision,allocatable :: vx(:),vy(:),vz(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat





    eps=0.001
    ALLOCATE(vx(ip00),vy(ip00),vz(ip00))
!
    n0=npc(l)
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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nci  = inc(1,0,0)
    ncj  = inc(0,1,0)
    nck  = inc(0,0,1)

!     composantes de la vitesse
!
    imin=i1m1
    imax=i2
    jmin=j1m1
    jmax=j2
    kmin=k1m1
    kmax=k2
    if (equat(3:5).eq.'2di') then
       imin=i1
       imax=i2m1
    endif
    if (equat(3:5).eq.'2dj') then
       jmin=j1
       jmax=j2m1
    endif
    if (equat(3:5).eq.'2dk') then
       kmin=k1
       kmax=k2m1
    endif
!
    do k=kmin,kmax
       do j=jmin,jmax
          ind1=ind(imin,j,k)
          ind2=ind(imax,j,k)
          do n=ind1,ind2
             m=n-n0
             vx(m)=s(n,2)/s(n,1)
             vy(m)=s(n,3)/s(n,1)
             vz(m)=s(n,4)/s(n,1)
          enddo
       enddo
    enddo
!
!     initialisation
    ind1=ind(id1(l),jd1(l),kd1(l))-n0
    ind2=ind(id2(l),jd2(l),kd2(l))-n0
    do m=ind1,ind2
       dvxx(m)=0.
       dvxy(m)=0.
       dvxz(m)=0.
       dvyx(m)=0.
       dvyy(m)=0.
       dvyz(m)=0.
       dvzx(m)=0.
       dvzy(m)=0.
       dvzz(m)=0.
       dtx(m)=0.
       dty(m)=0.
       dtz(m)=0.
    enddo
!
!     direction k (a travers les facettes k=cste)
!
    if (equat(3:5).ne.'2dk') then
!
       do k=k1p1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                dvxx(m)=dvxx(m)-sn(m,3,1)*(vx(m)+vx(m-nck))
                dvxy(m)=dvxy(m)-sn(m,3,2)*(vx(m)+vx(m-nck))
                dvxz(m)=dvxz(m)-sn(m,3,3)*(vx(m)+vx(m-nck))
                dvyx(m)=dvyx(m)-sn(m,3,1)*(vy(m)+vy(m-nck))
                dvyy(m)=dvyy(m)-sn(m,3,2)*(vy(m)+vy(m-nck))
                dvyz(m)=dvyz(m)-sn(m,3,3)*(vy(m)+vy(m-nck))
                dvzx(m)=dvzx(m)-sn(m,3,1)*(vz(m)+vz(m-nck))
                dvzy(m)=dvzy(m)-sn(m,3,2)*(vz(m)+vz(m-nck))
                dvzz(m)=dvzz(m)-sn(m,3,3)*(vz(m)+vz(m-nck))
                dtx(m)=dtx(m)-sn(m,3,1)*(temp(n)+temp(n-nck))
                dty(m)=dty(m)-sn(m,3,2)*(temp(n)+temp(n-nck))
                dtz(m)=dtz(m)-sn(m,3,3)*(temp(n)+temp(n-nck))
             enddo
!!$OMP SIMD
             do n=ind1,ind2
                m=n-n0
                dvxx(m-nck)=dvxx(m-nck)+sn(m,3,1)*(vx(m)+vx(m-nck))
                dvxy(m-nck)=dvxy(m-nck)+sn(m,3,2)*(vx(m)+vx(m-nck))
                dvxz(m-nck)=dvxz(m-nck)+sn(m,3,3)*(vx(m)+vx(m-nck))
                dvyx(m-nck)=dvyx(m-nck)+sn(m,3,1)*(vy(m)+vy(m-nck))
                dvyy(m-nck)=dvyy(m-nck)+sn(m,3,2)*(vy(m)+vy(m-nck))
                dvyz(m-nck)=dvyz(m-nck)+sn(m,3,3)*(vy(m)+vy(m-nck))
                dvzx(m-nck)=dvzx(m-nck)+sn(m,3,1)*(vz(m)+vz(m-nck))
                dvzy(m-nck)=dvzy(m-nck)+sn(m,3,2)*(vz(m)+vz(m-nck))
                dvzz(m-nck)=dvzz(m-nck)+sn(m,3,3)*(vz(m)+vz(m-nck))
                dtx(m-nck)=dtx(m-nck)+sn(m,3,1)*(temp(n)+temp(n-nck))
                dty(m-nck)=dty(m-nck)+sn(m,3,2)*(temp(n)+temp(n-nck))
                dtz(m-nck)=dtz(m-nck)+sn(m,3,3)*(temp(n)+temp(n-nck))
             enddo
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = ind(i1  ,j,k1)
          ind2 = ind(i2m1,j,k1)
          do n=ind1,ind2
             m=n-n0
             dvxx(m)=dvxx(m)-sn(m,3,1)*2*vx(m-nck)
             dvxy(m)=dvxy(m)-sn(m,3,2)*2*vx(m-nck)
             dvxz(m)=dvxz(m)-sn(m,3,3)*2*vx(m-nck)
             dvyx(m)=dvyx(m)-sn(m,3,1)*2*vy(m-nck)
             dvyy(m)=dvyy(m)-sn(m,3,2)*2*vy(m-nck)
             dvyz(m)=dvyz(m)-sn(m,3,3)*2*vy(m-nck)
             dvzx(m)=dvzx(m)-sn(m,3,1)*2*vz(m-nck)
             dvzy(m)=dvzy(m)-sn(m,3,2)*2*vz(m-nck)
             dvzz(m)=dvzz(m)-sn(m,3,3)*2*vz(m-nck)
             dtx(m)=dtx(m)-sn(m,3,1)*2*temp(n-nck)
             dty(m)=dty(m)-sn(m,3,2)*2*temp(n-nck)
             dtz(m)=dtz(m)-sn(m,3,3)*2*temp(n-nck)
          enddo
       enddo
!
       do j=j1,j2m1
          ind1 = ind(i1  ,j,k2)
          ind2 = ind(i2m1,j,k2)
          do n=ind1,ind2
             m=n-n0
             dvxx(m-nck)=dvxx(m-nck)+sn(m,3,1)*2*vx(m)
             dvxy(m-nck)=dvxy(m-nck)+sn(m,3,2)*2*vx(m)
             dvxz(m-nck)=dvxz(m-nck)+sn(m,3,3)*2*vx(m)
             dvyx(m-nck)=dvyx(m-nck)+sn(m,3,1)*2*vy(m)
             dvyy(m-nck)=dvyy(m-nck)+sn(m,3,2)*2*vy(m)
             dvyz(m-nck)=dvyz(m-nck)+sn(m,3,3)*2*vy(m)
             dvzx(m-nck)=dvzx(m-nck)+sn(m,3,1)*2*vz(m)
             dvzy(m-nck)=dvzy(m-nck)+sn(m,3,2)*2*vz(m)
             dvzz(m-nck)=dvzz(m-nck)+sn(m,3,3)*2*vz(m)
             dtx(m-nck)=dtx(m-nck)+sn(m,3,1)*2*temp(n)
             dty(m-nck)=dty(m-nck)+sn(m,3,2)*2*temp(n)
             dtz(m-nck)=dtz(m-nck)+sn(m,3,3)*2*temp(n)
          enddo
       enddo
!
    endif
!
!     direction j (a travers les facettes j=cste)
!
    if (equat(3:5).ne.'2dj') then
!
       do j=j1p1,j2m1
          do k=k1,k2m1
             ind1 = ind(i1  ,j,k)
             ind2 = ind(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0
                dvxx(m)=dvxx(m)-sn(m,2,1)*(vx(m)+vx(m-ncj))
                dvxy(m)=dvxy(m)-sn(m,2,2)*(vx(m)+vx(m-ncj))
                dvxz(m)=dvxz(m)-sn(m,2,3)*(vx(m)+vx(m-ncj))
                dvyx(m)=dvyx(m)-sn(m,2,1)*(vy(m)+vy(m-ncj))
                dvyy(m)=dvyy(m)-sn(m,2,2)*(vy(m)+vy(m-ncj))
                dvyz(m)=dvyz(m)-sn(m,2,3)*(vy(m)+vy(m-ncj))
                dvzx(m)=dvzx(m)-sn(m,2,1)*(vz(m)+vz(m-ncj))
                dvzy(m)=dvzy(m)-sn(m,2,2)*(vz(m)+vz(m-ncj))
                dvzz(m)=dvzz(m)-sn(m,2,3)*(vz(m)+vz(m-ncj))
                dtx(m)=dtx(m)-sn(m,2,1)*(temp(n)+temp(n-ncj))
                dty(m)=dty(m)-sn(m,2,2)*(temp(n)+temp(n-ncj))
                dtz(m)=dtz(m)-sn(m,2,3)*(temp(n)+temp(n-ncj))
             enddo
             do n=ind1,ind2
                m=n-n0
                dvxx(m-ncj)=dvxx(m-ncj)+sn(m,2,1)*(vx(m)+vx(m-ncj))
                dvxy(m-ncj)=dvxy(m-ncj)+sn(m,2,2)*(vx(m)+vx(m-ncj))
                dvxz(m-ncj)=dvxz(m-ncj)+sn(m,2,3)*(vx(m)+vx(m-ncj))
                dvyx(m-ncj)=dvyx(m-ncj)+sn(m,2,1)*(vy(m)+vy(m-ncj))
                dvyy(m-ncj)=dvyy(m-ncj)+sn(m,2,2)*(vy(m)+vy(m-ncj))
                dvyz(m-ncj)=dvyz(m-ncj)+sn(m,2,3)*(vy(m)+vy(m-ncj))
                dvzx(m-ncj)=dvzx(m-ncj)+sn(m,2,1)*(vz(m)+vz(m-ncj))
                dvzy(m-ncj)=dvzy(m-ncj)+sn(m,2,2)*(vz(m)+vz(m-ncj))
                dvzz(m-ncj)=dvzz(m-ncj)+sn(m,2,3)*(vz(m)+vz(m-ncj))
                dtx(m-ncj)=dtx(m-ncj)+sn(m,2,1)*(temp(n)+temp(n-ncj))
                dty(m-ncj)=dty(m-ncj)+sn(m,2,2)*(temp(n)+temp(n-ncj))
                dtz(m-ncj)=dtz(m-ncj)+sn(m,2,3)*(temp(n)+temp(n-ncj))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1  ,j1,k)
          ind2 = ind(i2m1,j1,k)
          do n=ind1,ind2
             m=n-n0
             dvxx(m)=dvxx(m)-sn(m,2,1)*2*vx(m-ncj)
             dvxy(m)=dvxy(m)-sn(m,2,2)*2*vx(m-ncj)
             dvxz(m)=dvxz(m)-sn(m,2,3)*2*vx(m-ncj)
             dvyx(m)=dvyx(m)-sn(m,2,1)*2*vy(m-ncj)
             dvyy(m)=dvyy(m)-sn(m,2,2)*2*vy(m-ncj)
             dvyz(m)=dvyz(m)-sn(m,2,3)*2*vy(m-ncj)
             dvzx(m)=dvzx(m)-sn(m,2,1)*2*vz(m-ncj)
             dvzy(m)=dvzy(m)-sn(m,2,2)*2*vz(m-ncj)
             dvzz(m)=dvzz(m)-sn(m,2,3)*2*vz(m-ncj)
             dtx(m)=dtx(m)-sn(m,2,1)*2*temp(n-ncj)
             dty(m)=dty(m)-sn(m,2,2)*2*temp(n-ncj)
             dtz(m)=dtz(m)-sn(m,2,3)*2*temp(n-ncj)
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1  ,j2,k)
          ind2 = ind(i2m1,j2,k)
          do n=ind1,ind2
             m=n-n0
             dvxx(m-ncj)=dvxx(m-ncj)+sn(m,2,1)*2*vx(m)
             dvxy(m-ncj)=dvxy(m-ncj)+sn(m,2,2)*2*vx(m)
             dvxz(m-ncj)=dvxz(m-ncj)+sn(m,2,3)*2*vx(m)
             dvyx(m-ncj)=dvyx(m-ncj)+sn(m,2,1)*2*vy(m)
             dvyy(m-ncj)=dvyy(m-ncj)+sn(m,2,2)*2*vy(m)
             dvyz(m-ncj)=dvyz(m-ncj)+sn(m,2,3)*2*vy(m)
             dvzx(m-ncj)=dvzx(m-ncj)+sn(m,2,1)*2*vz(m)
             dvzy(m-ncj)=dvzy(m-ncj)+sn(m,2,2)*2*vz(m)
             dvzz(m-ncj)=dvzz(m-ncj)+sn(m,2,3)*2*vz(m)
             dtx(m-ncj)=dtx(m-ncj)+sn(m,2,1)*2*temp(n)
             dty(m-ncj)=dty(m-ncj)+sn(m,2,2)*2*temp(n)
             dtz(m-ncj)=dtz(m-ncj)+sn(m,2,3)*2*temp(n)
          enddo
       enddo
!
    endif
!
!     direction i (a travers les facettes i=cste)
!
    if (equat(3:5).ne.'2di') then
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = ind(i1p1,j,k)
             ind2 = ind(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0
                dvxx(m)=dvxx(m)-sn(m,1,1)*(vx(m)+vx(m-nci))
                dvxy(m)=dvxy(m)-sn(m,1,2)*(vx(m)+vx(m-nci))
                dvxz(m)=dvxz(m)-sn(m,1,3)*(vx(m)+vx(m-nci))
                dvyx(m)=dvyx(m)-sn(m,1,1)*(vy(m)+vy(m-nci))
                dvyy(m)=dvyy(m)-sn(m,1,2)*(vy(m)+vy(m-nci))
                dvyz(m)=dvyz(m)-sn(m,1,3)*(vy(m)+vy(m-nci))
                dvzx(m)=dvzx(m)-sn(m,1,1)*(vz(m)+vz(m-nci))
                dvzy(m)=dvzy(m)-sn(m,1,2)*(vz(m)+vz(m-nci))
                dvzz(m)=dvzz(m)-sn(m,1,3)*(vz(m)+vz(m-nci))
                dtx(m)=dtx(m)-sn(m,1,1)*(temp(n)+temp(n-nci))
                dty(m)=dty(m)-sn(m,1,2)*(temp(n)+temp(n-nci))
                dtz(m)=dtz(m)-sn(m,1,3)*(temp(n)+temp(n-nci))
             enddo
             do n=ind1,ind2
                m=n-n0
                dvxx(m-nci)=dvxx(m-nci)+sn(m,1,1)*(vx(m)+vx(m-nci))
                dvxy(m-nci)=dvxy(m-nci)+sn(m,1,2)*(vx(m)+vx(m-nci))
                dvxz(m-nci)=dvxz(m-nci)+sn(m,1,3)*(vx(m)+vx(m-nci))
                dvyx(m-nci)=dvyx(m-nci)+sn(m,1,1)*(vy(m)+vy(m-nci))
                dvyy(m-nci)=dvyy(m-nci)+sn(m,1,2)*(vy(m)+vy(m-nci))
                dvyz(m-nci)=dvyz(m-nci)+sn(m,1,3)*(vy(m)+vy(m-nci))
                dvzx(m-nci)=dvzx(m-nci)+sn(m,1,1)*(vz(m)+vz(m-nci))
                dvzy(m-nci)=dvzy(m-nci)+sn(m,1,2)*(vz(m)+vz(m-nci))
                dvzz(m-nci)=dvzz(m-nci)+sn(m,1,3)*(vz(m)+vz(m-nci))
                dtx(m-nci)=dtx(m-nci)+sn(m,1,1)*(temp(n)+temp(n-nci))
                dty(m-nci)=dty(m-nci)+sn(m,1,2)*(temp(n)+temp(n-nci))
                dtz(m-nci)=dtz(m-nci)+sn(m,1,3)*(temp(n)+temp(n-nci))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i1,j1  ,k)
          ind2 = ind(i1,j2m1,k)
          do n=ind1,ind2,ncj
             m=n-n0
             dvxx(m)=dvxx(m)-sn(m,1,1)*2*vx(m-nci)
             dvxy(m)=dvxy(m)-sn(m,1,2)*2*vx(m-nci)
             dvxz(m)=dvxz(m)-sn(m,1,3)*2*vx(m-nci)
             dvyx(m)=dvyx(m)-sn(m,1,1)*2*vy(m-nci)
             dvyy(m)=dvyy(m)-sn(m,1,2)*2*vy(m-nci)
             dvyz(m)=dvyz(m)-sn(m,1,3)*2*vy(m-nci)
             dvzx(m)=dvzx(m)-sn(m,1,1)*2*vz(m-nci)
             dvzy(m)=dvzy(m)-sn(m,1,2)*2*vz(m-nci)
             dvzz(m)=dvzz(m)-sn(m,1,3)*2*vz(m-nci)
             dtx(m)=dtx(m)-sn(m,1,1)*2*temp(n-nci)
             dty(m)=dty(m)-sn(m,1,2)*2*temp(n-nci)
             dtz(m)=dtz(m)-sn(m,1,3)*2*temp(n-nci)
          enddo
       enddo
!
       do k=k1,k2m1
          ind1 = ind(i2,j1  ,k)
          ind2 = ind(i2,j2m1,k)
          do n=ind1,ind2,ncj
             m=n-n0
             dvxx(m-nci)=dvxx(m-nci)+sn(m,1,1)*2*vx(m)
             dvxy(m-nci)=dvxy(m-nci)+sn(m,1,2)*2*vx(m)
             dvxz(m-nci)=dvxz(m-nci)+sn(m,1,3)*2*vx(m)
             dvyx(m-nci)=dvyx(m-nci)+sn(m,1,1)*2*vy(m)
             dvyy(m-nci)=dvyy(m-nci)+sn(m,1,2)*2*vy(m)
             dvyz(m-nci)=dvyz(m-nci)+sn(m,1,3)*2*vy(m)
             dvzx(m-nci)=dvzx(m-nci)+sn(m,1,1)*2*vz(m)
             dvzy(m-nci)=dvzy(m-nci)+sn(m,1,2)*2*vz(m)
             dvzz(m-nci)=dvzz(m-nci)+sn(m,1,3)*2*vz(m)
             dtx(m-nci)=dtx(m-nci)+sn(m,1,1)*2*temp(n)
             dty(m-nci)=dty(m-nci)+sn(m,1,2)*2*temp(n)
             dtz(m-nci)=dtz(m-nci)+sn(m,1,3)*2*temp(n)
          enddo
       enddo
!
    endif
!
!-----calcul du gradient :
!
    ind1 = ind(i1  ,j1  ,k1  )
    ind2 = ind(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       m=n-n0
!      le coefficient 1/2 provient de la moyenne de vx,vy,vz ou t
       ts=sign(0.5,-vol(n))
       vols=(0.5+ts)*eps+(0.5-ts)*vol(n)
       c0=0.5/vols
!
       dvxx(m)=dvxx(m)*c0
       dvxy(m)=dvxy(m)*c0
       dvxz(m)=dvxz(m)*c0
       dvyx(m)=dvyx(m)*c0
       dvyy(m)=dvyy(m)*c0
       dvyz(m)=dvyz(m)*c0
       dvzx(m)=dvzx(m)*c0
       dvzy(m)=dvzy(m)*c0
       dvzz(m)=dvzz(m)*c0
       dtx(m)=dtx(m)*c0
       dty(m)=dty(m)*c0
       dtz(m)=dtz(m)*c0
    enddo

    DEALLOCATE(vx,vy,vz)

    return
  contains
    function    ind(i,j,k)
      implicit none
  integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
  integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine zgrad
end module mod_zgrad
