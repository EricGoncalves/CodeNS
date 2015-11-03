module mod_zpres
  implicit none
contains
  subroutine zpres( &
       l,v, &
       pression,ztemp,cson)
!
!***********************************************************************
!
!_DA  DATE_C :    novembre 2006 -- Eric GONCALVES / LEGI
!
!     ACT
!_A     Calcul de la pression, temperature et vitesse du son
!_A     par la loi d'etat des gaz raides.
!_A
!
!_I    l          : arg int              ; numero de domaine
!_I    v          : arg real(ip11,ip60 ) ; variables de calcul
!
!     OUT
!_O    pression   : arg real(ip11      ) ; pression statique
!_O    ztemp      : arg real(ip11      ) ; temperature
!_O    cson       : arg real(ip11      ) ; vitesse du son
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu 
    implicit none
    integer          ::       i,     i1,   i1m1,   i1p1,     i2
    integer          ::    i2m1,   ind1,   ind2,isortie,      j
    integer          ::      j1,   j1m1,   j1p1,     j2,   j2m1
    integer          ::       k,     k1,   k1m1,   k1p1,     k2
    integer          ::    k2m1,      l,      m,      n,    n0c
    integer          ::     nid,   nijd,    njd
    double precision ::     cson,pression,    rhoe,       v,   ztemp
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension pression(ip11),ztemp(ip11),cson(ip11)
!
    n0c=npc(l)
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
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             rhoe=v(n,5)-0.5*(v(n,2)**2+v(n,3)**2+v(n,4)**2)/v(n,1)
             pression(n)=gam1*(rhoe-pinfl)
!     &      +(5./3.-gam)*v(n,6)
             ztemp(n)=gam*pression(n)/v(n,1)
             cson(n)=sqrt(ztemp(n))
          enddo
       enddo
    enddo
!
    isortie=0
    if(isortie.eq.1) then
       write(6,'("===>zpres")')
       k=1
       i=44
       do j=j1,j2m1
          n=indc(i,j,k)
          m=n-n0c
          write(6,'(i4,i6,4(1pe12.4))') &
               j,n,v(n,1),pression(n),ztemp(n),cson(n)
       enddo
    endif
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine zpres
end module mod_zpres
