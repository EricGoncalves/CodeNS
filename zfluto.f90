module mod_zfluto
  implicit none
contains
  subroutine zfluto( &
       l,mu,mut,toxx,toxy,toxz,toyy,toyz,tozz, &
       qcx,qcy,qcz, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       dtx,dty,dtz)
!
!***********************************************************************
!
!     ACT
!_A    Calcul du tenseur des contraintes visqueuses et du vecteur flux de chaleur.
!
!_I    l          : arg int              ; numero de domaine
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    dvxx       : arg real(ip00      ) ; composante xx du gradient de vitesse
!_I    dvxy       : arg real(ip00      ) ; composante xy du gradient de vitesse
!_I    dvxz       : arg real(ip00      ) ; composante xz du gradient de vitesse
!_I    dvyx       : arg real(ip00      ) ; composante yx du gradient de vitesse
!_I    dvyy       : arg real(ip00      ) ; composante yy du gradient de vitesse
!_I    dvyz       : arg real(ip00      ) ; composante yz du gradient de vitesse
!_I    dvzx       : arg real(ip00      ) ; composante zx du gradient de vitesse
!_I    dvzy       : arg real(ip00      ) ; composante zy du gradient de vitesse
!_I    dvzz       : arg real(ip00      ) ; composante zz du gradient de vitesse
!_I    dtx        : arg real(ip00      ) ; composante x du gradient de temp
!_I    dty        : arg real(ip00      ) ; composante y du gradient de temp
!_I    dtz        : arg real(ip00      ) ; composante z du gradient de temp
!_I    cp         : com real             ; chal spec a pres cste adim
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
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
!_O    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_O                                        contraintes visqueuses
!_O    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_O                                        contraintes visqueuses
!_O    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_O                                        contraintes visqueuses
!_O    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_O                                        contraintes visqueuses
!_O    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_O                                        contraintes visqueuses
!_O    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_O                                        contraintes visqueuses
!_O    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_O    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_O    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    implicit none
    integer          ::    i,  i1,  i2,i2m1,ind1
    integer          :: ind2,   j,  j1,  j2,j2m1
    integer          ::    k,  k1,  k2,k2m1,   l
    integer          ::    m,   n,  n0, nid,nijd
    integer          ::  njd
    double precision ::        ds3, dtx(ip00), dty(ip00), dtz(ip00),dvxx(ip00)
    double precision :: dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00),dvyz(ip00)
    double precision :: dvzx(ip00),dvzy(ip00),dvzz(ip00),  mu(ip12), mut(ip12)
    double precision ::  qcx(ip12), qcy(ip12), qcz(ip12),toxx(ip12),toxy(ip12)
    double precision :: toxz(ip12),toyy(ip12),toyz(ip12),tozz(ip12)
!
!-----------------------------------------------------------------------
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
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    ds3=2./3.
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=ind(i1  ,j,k)
          ind2=ind(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0
             toxx(n)=ds3*(mu(n)+mut(n))*(2*dvxx(m)-dvyy(m)-dvzz(m))
             toyy(n)=ds3*(mu(n)+mut(n))*(2*dvyy(m)-dvxx(m)-dvzz(m))
             tozz(n)=ds3*(mu(n)+mut(n))*(2*dvzz(m)-dvxx(m)-dvyy(m))
             toxy(n)=    (mu(n)+mut(n))*(dvxy(m)+dvyx(m))
             toxz(n)=    (mu(n)+mut(n))*(dvxz(m)+dvzx(m))
             toyz(n)=    (mu(n)+mut(n))*(dvyz(m)+dvzy(m))
             qcx (n)=cp*(mu(n)/pr+mut(n)/prt)*dtx(m)
             qcy (n)=cp*(mu(n)/pr+mut(n)/prt)*dty(m)
             qcz (n)=cp*(mu(n)/pr+mut(n)/prt)*dtz(m)
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
  end subroutine zfluto
end module mod_zfluto
