module mod_at_lecdist
  implicit none
contains
  subroutine at_lecdist( &
       ldismx, &
       dist,mnpar)
!
!***********************************************************************
!
!     ACT
!_A    lecture des distances pour un domaine
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use maillage
    implicit none
    integer          ::           i,         i1,         i2,       i2m1,        ird
    integer          ::           j,         j1,         j2,       j2m1,          k
    integer          ::          k1,         k2,       k2m1,          l,     ldismx
    integer          :: mnpar(ip12),         n0,        nid,       nijd,        njd
    double precision :: dist(ip12)
!
!-----------------------------------------------------------------------
!
!

!
    ird=98
    open(ird,file='fdist',form='unformatted',status='old',err=50)
!
    ldismx=0
    do l=1,lzx
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
       read(ird,end=10,err=11) &
            (((dist (ind(i,j,k)),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
       read(ird,end=11,err=12) &
            (((mnpar(ind(i,j,k)),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
       ldismx=ldismx+1
    end do
!
10  continue
    close(ird)
    write(imp,'("===>at_lecdist: fin lecture fichier= fdist  nb domaines lus=",i3)')ldismx
    return
!
50  continue
    write(imp,'("!!!at_lecdist: erreur ouverture fichier fdist")')
    stop
11  continue
    write(imp,'("!!!at_lecdist: erreur lecture fichier fdist")')
    stop
12  continue
    write(imp,'("!!!at_lecdist: fin prematuree fichier fdist")')
    stop
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine at_lecdist


end module mod_at_lecdist
