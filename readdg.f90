module mod_readdg
  implicit none
contains
  subroutine readdg( &
       l,kdg,x,y,z)
!
!***********************************************************************
!
!     ACT
!_A    Lecture des coordonnees x, y, z en tout noeud (non fictif)
!_A    d'un domaine structure.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kdg        : arg int              ; unite logique, maillage
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
!
!     OUT
!_O    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_O    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_O    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use maillage
    use mod_mpi
    implicit none
    integer          ::    i,  i1,  i2,   j,  j1
    integer          ::   j2,   k,  k1,  k2, kdg,pos
    integer          ::    l,   n, nid,nijd, njd
    double precision :: x(ip21),y(ip21),z(ip21)
    logical          :: ecri
!
!-----------------------------------------------------------------------
!
    character(len=1) :: coord
!

!
    ecri=.false.
!      ecri=.true.
    pos=FTELL(kdg) 
!
    call START_KEEP_ORDER(pos)
    call my_FSEEK(kdg, pos)
!
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    coord='x'
    read(kdg,err=13) &
         (((x(indn(i,j,k)),i=i1,i2),j=j1,j2),k=k1,k2)
    coord='y'
    read(kdg,err=13) &
         (((y(indn(i,j,k)),i=i1,i2),j=j1,j2),k=k1,k2)
    coord='z'
    read(kdg,err=13) &
         (((z(indn(i,j,k)),i=i1,i2),j=j1,j2),k=k1,k2)
!
    pos=FTELL(kdg)
    call END_KEEP_ORDER(pos)

    return
!
13  continue
    write(imp,'(/,"!!!readdg: probleme lecture maillage ",/,a1,3x,"l=",i3)') coord,l
    write(imp,'(10x,"i2=",i5,3x,"j2=",i5,3x,"k2=",i5)')i2,j2,k2
    stop
!
    if(ecri.and.rank+1==l) then
!       ecriture plaque plane 1 domaine
       open(out  ,file='fout')
       k=1
       do k=1,2
          do i=1,i2,50
             write(out,'("===>readdg: l=",i3)')l
             write(out,'("   i   j   k       n",t25,"x",t36,"y",t47,"z")')
             do j=1,j2,2
                n=indn(i,j,k)
                write(out,'(3i4,i8,5(1pe11.3))') &
                     i,j,k,n,x(n),y(n),z(n)
             enddo
          enddo
       enddo
       close(out)
    endif
!
    return
  contains
    function    indn(i,j,k)
      implicit none
      integer          ::    i,indn,   j,   k
      indn=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indn
  end subroutine readdg
end module mod_readdg
