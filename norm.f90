module mod_norm
  implicit none
contains
  subroutine norm( &
       l,x,y,z, &
       equat,ndt, &
       imin,imax,jmin,jmax,kmin,kmax, &
       tnix,tniy,tniz, &
       tnjx,tnjy,tnjz, &
       tnkx,tnky,tnkz)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des normales aux facettes, multipliees par l'aire de la facette,
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
!_I    ndt        : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!
!     OUT
!_O    tnix       : arg real(ndt       ) ; composante en x du vect normal a une
!_O                                        facette i et de norme egale a la
!_O                                        surface de celle-ci
!_O    tniy       : arg real(ndt       ) ; composante en y du vect normal a une
!_O                                        facette i et de norme egale a la
!_O                                        surface de celle-ci
!_O    tniz       : arg real(ndt       ) ; composante en z du vect normal a une
!_O                                        facette i et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnjx       : arg real(ndt       ) ; composante en x du vect normal a une
!_O                                        facette j et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnjy       : arg real(ndt       ) ; composante en y du vect normal a une
!_O                                        facette j et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnjz       : arg real(ndt       ) ; composante en z du vect normal a une
!_O                                        facette j et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnkx       : arg real(ndt       ) ; composante en x du vect normal a une
!_O                                        facette k et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnky       : arg real(ndt       ) ; composante en y du vect normal a une
!_O                                        facette k et de norme egale a la
!_O                                        surface de celle-ci
!_O    tnkz       : arg real(ndt       ) ; composante en z du vect normal a une
!_O                                        facette k et de norme egale a la
!_O                                        surface de celle-ci
!
!_C    Le vecteur normal est le demi-produit vectoriel des diagonales
!_C    de la facette.
!_C    Chaque vecteur normal est determine sur une facette d'indice
!_C    i,j ou k constant et est oriente dans le sens de cet indice croissant.
!_C    Pour s'assurer de l'orientation on effectue le produit scalaire
!_C    de la normale avec le vecteur d'indice choisi croissant.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
  integer          ::      i,    i1,  i1p1,    id,  imax
  integer          :: imaxm1,  imin,   ior,     j,    j1
  integer          ::   j1p1,    jd,  jmax,jmaxm1,  jmin
  integer          ::    jor,     k,    k1,  k1p1,    kd
  integer          ::   kmax,kmaxm1,  kmin,   kor,     l
  integer          ::      m,   mor,     n,    n0,   ndt
  integer          ::    nid,  nijd,   njd,   nor
  double precision ::       dx1,      dx2,      dy1,      dy2,      dz1
  double precision ::       dz2,      psi,      psj,      psk,tnix(ndt)
  double precision :: tniy(ndt),tniz(ndt),tnjx(ndt),tnjy(ndt),tnjz(ndt)
  double precision :: tnkx(ndt),tnky(ndt),tnkz(ndt),       vx,       vy
  double precision ::        vz,  x(ip21),  y(ip21),  z(ip21)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!


!
    n0=npn(l)
    i1=ii1(l)
    j1=jj1(l)
    k1=kk1(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
!
    imaxm1=imax-1
    jmaxm1=jmax-1
    kmaxm1=kmax-1
!
!-----indice pour test d'orientation de nds
!
    ior=i1p1
    jor=j1p1
    kor=k1p1
    if(equat(3:5).eq.'2di') ior=i1
    if(equat(3:5).eq.'2dj') jor=j1
    if(equat(3:5).eq.'2dk') kor=k1
    if(equat(3:5).eq.'2xk') kor=k1
    nor=ind(ior,jor,kor)
    mor=nor-n0
!
!-----facettes de type k=constante
!
!     test d'orientation de nds
    dx1=x(nor+inc(1,1,0))-x(nor)
    dy1=y(nor+inc(1,1,0))-y(nor)
    dz1=z(nor+inc(1,1,0))-z(nor)
    dx2=x(nor+inc(0,1,0))-x(nor+inc(1,0,0))
    dy2=y(nor+inc(0,1,0))-y(nor+inc(1,0,0))
    dz2=z(nor+inc(0,1,0))-z(nor+inc(1,0,0))
    tnkx(mor)=dy1*dz2-dz1*dy2
    tnky(mor)=dz1*dx2-dx1*dz2
    tnkz(mor)=dx1*dy2-dy1*dx2
    vx = x(nor+inc(0,0,1))-x(nor) + &
         x(nor+inc(1,0,1))-x(nor+inc(1,0,0)) + &
         x(nor+inc(0,1,1))-x(nor+inc(0,1,0)) + &
         x(nor+inc(1,1,1))-x(nor+inc(1,1,0))
    vy = y(nor+inc(0,0,1))-y(nor) + &
         y(nor+inc(1,0,1))-y(nor+inc(1,0,0)) + &
         y(nor+inc(0,1,1))-y(nor+inc(0,1,0)) + &
         y(nor+inc(1,1,1))-y(nor+inc(1,1,0))
    vz = z(nor+inc(0,0,1))-z(nor) + &
         z(nor+inc(1,0,1))-z(nor+inc(1,0,0)) + &
         z(nor+inc(0,1,1))-z(nor+inc(0,1,0)) + &
         z(nor+inc(1,1,1))-z(nor+inc(1,1,0))
    psk = tnkx(mor)*vx + tnky(mor)*vy + tnkz(mor)*vz
!
    do k=kmin,kmax
       do j=jmin,jmaxm1
          do i=imin,imaxm1
             n=ind(i,j,k)
             m=n-n0
             dx1=x(n+inc(1,1,0))-x(n)
             dy1=y(n+inc(1,1,0))-y(n)
             dz1=z(n+inc(1,1,0))-z(n)
             dx2=x(n+inc(0,1,0))-x(n+inc(1,0,0))
             dy2=y(n+inc(0,1,0))-y(n+inc(1,0,0))
             dz2=z(n+inc(0,1,0))-z(n+inc(1,0,0))
             tnkx(m)=dy1*dz2-dz1*dy2
             tnky(m)=dz1*dx2-dx1*dz2
             tnkz(m)=dx1*dy2-dy1*dx2
             tnkx(m)=sign(.5,psk)*tnkx(m)
             tnky(m)=sign(.5,psk)*tnky(m)
             tnkz(m)=sign(.5,psk)*tnkz(m)
          enddo
       enddo
    enddo
!
!-----facettes de type j=constante
!
!     test d'orientation de nds
    dx1=x(nor+inc(1,0,1))-x(nor)
    dy1=y(nor+inc(1,0,1))-y(nor)
    dz1=z(nor+inc(1,0,1))-z(nor)
    dx2=x(nor+inc(0,0,1))-x(nor+inc(1,0,0))
    dy2=y(nor+inc(0,0,1))-y(nor+inc(1,0,0))
    dz2=z(nor+inc(0,0,1))-z(nor+inc(1,0,0))
    tnjx(mor)=dy1*dz2-dz1*dy2
    tnjy(mor)=dz1*dx2-dx1*dz2
    tnjz(mor)=dx1*dy2-dy1*dx2
    vx = x(nor+inc(0,1,0))-x(nor) + &
         x(nor+inc(1,1,0))-x(nor+inc(1,0,0)) + &
         x(nor+inc(0,1,1))-x(nor+inc(0,0,1)) + &
         x(nor+inc(1,1,1))-x(nor+inc(1,0,1))
    vy = y(nor+inc(0,1,0))-y(nor) + &
         y(nor+inc(1,1,0))-y(nor+inc(1,0,0)) + &
         y(nor+inc(0,1,1))-y(nor+inc(0,0,1)) + &
         y(nor+inc(1,1,1))-y(nor+inc(1,0,1))
    vz = z(nor+inc(0,1,0))-z(nor) + &
         z(nor+inc(1,1,0))-z(nor+inc(1,0,0)) + &
         z(nor+inc(0,1,1))-z(nor+inc(0,0,1)) + &
         z(nor+inc(1,1,1))-z(nor+inc(1,0,1))
    psj = tnjx(mor)*vx + tnjy(mor)*vy + tnjz(mor)*vz
!
    do j=jmin,jmax
       do k=kmin,kmaxm1
          do i=imin,imaxm1
             n=ind(i,j,k)
             m=n-n0
             dx1=x(n+inc(1,0,1))-x(n)
             dy1=y(n+inc(1,0,1))-y(n)
             dz1=z(n+inc(1,0,1))-z(n)
             dx2=x(n+inc(0,0,1))-x(n+inc(1,0,0))
             dy2=y(n+inc(0,0,1))-y(n+inc(1,0,0))
             dz2=z(n+inc(0,0,1))-z(n+inc(1,0,0))
             tnjx(m)=dy1*dz2-dz1*dy2
             tnjy(m)=dz1*dx2-dx1*dz2
             tnjz(m)=dx1*dy2-dy1*dx2
             tnjx(m)=sign(.5,psj)*tnjx(m)
             tnjy(m)=sign(.5,psj)*tnjy(m)
             tnjz(m)=sign(.5,psj)*tnjz(m)
          enddo
       enddo
    enddo
!
!-----facettes de type i=constante
!
!     test d'orientation de nds
    dx1=x(nor+inc(0,1,1))-x(nor)
    dy1=y(nor+inc(0,1,1))-y(nor)
    dz1=z(nor+inc(0,1,1))-z(nor)
    dx2=x(nor+inc(0,0,1))-x(nor+inc(0,1,0))
    dy2=y(nor+inc(0,0,1))-y(nor+inc(0,1,0))
    dz2=z(nor+inc(0,0,1))-z(nor+inc(0,1,0))
    tnix(mor)=dy1*dz2-dz1*dy2
    tniy(mor)=dz1*dx2-dx1*dz2
    tniz(mor)=dx1*dy2-dy1*dx2
    vx = x(nor+inc(1,0,0))-x(nor) + &
         x(nor+inc(1,1,0))-x(nor+inc(0,1,0)) + &
         x(nor+inc(1,0,1))-x(nor+inc(0,0,1)) + &
         x(nor+inc(1,1,1))-x(nor+inc(0,1,1))
    vy = y(nor+inc(1,0,0))-y(nor) + &
         y(nor+inc(1,1,0))-y(nor+inc(0,1,0)) + &
         y(nor+inc(1,0,1))-y(nor+inc(0,0,1)) + &
         y(nor+inc(1,1,1))-y(nor+inc(0,1,1))
    vz = z(nor+inc(1,0,0))-z(nor) + &
         z(nor+inc(1,1,0))-z(nor+inc(0,1,0)) + &
         z(nor+inc(1,0,1))-z(nor+inc(0,0,1)) + &
         z(nor+inc(1,1,1))-z(nor+inc(0,1,1))
    psi = tnix(mor)*vx + tniy(mor)*vy + tniz(mor)*vz
!
    do i=imin,imax
       do k=kmin,kmaxm1
          do j=jmin,jmaxm1
             n=ind(i,j,k)
             m=n-n0
             dx1=x(n+inc(0,1,1))-x(n)
             dy1=y(n+inc(0,1,1))-y(n)
             dz1=z(n+inc(0,1,1))-z(n)
             dx2=x(n+inc(0,0,1))-x(n+inc(0,1,0))
             dy2=y(n+inc(0,0,1))-y(n+inc(0,1,0))
             dz2=z(n+inc(0,0,1))-z(n+inc(0,1,0))
             tnix(m)=dy1*dz2-dz1*dy2
             tniy(m)=dz1*dx2-dx1*dz2
             tniz(m)=dx1*dy2-dy1*dx2
             tnix(m)=sign(.5,psi)*tnix(m)
             tniy(m)=sign(.5,psi)*tniy(m)
             tniz(m)=sign(.5,psi)*tniz(m)
          enddo
       enddo
    enddo
!
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
  end subroutine norm
end module mod_norm
