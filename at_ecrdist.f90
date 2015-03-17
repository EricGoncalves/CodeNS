module mod_at_ecrdist
  implicit none
contains
  subroutine at_ecrdist( &
       l0,         &
       dist,mnpar)
!
!***********************************************************************
!
!     ACT
!_A    ecriture des distances pour un domaine
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use maillage
    implicit none
    integer          ::     i,   i1,   i2, i2m1,  ind
    integer          ::   iwd,    j,   j1,   j2, j2m1
    integer          ::     k,   k1,   k2, k2m1,    l
    integer          ::    l0,mnpar,   n0,  nid, nijd
    integer          ::   njd
    double precision :: dist
!
!-----------------------------------------------------------------------
!
    character(len=8) :: nomfich
    dimension dist(ip12),mnpar(ip12)
!
    ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
    iwd=98
    if(l0.eq.0) then
!       ecriture tous les domaines dans "fdist"
!
       nomfich='fdist   '
       open(iwd,file=nomfich,form='unformatted',err=50)
!
       write(imp,'("===>at_ecrdist: distance tous domaines   fichier=",a8)')nomfich
!
       do l=1,lzx
          n0=npc(l)
          i1=ii1(l)
          i2=ii2(l)
          j1=jj1(l)
          j2=jj2(l)
          k1=kk1(l)
          k2=kk2(l)
          i2m1=i2-1
          j2m1=j2-1
          k2m1=k2-1
          nid = id2(l)-id1(l)+1
          njd = jd2(l)-jd1(l)+1
          nijd = nid*njd
!
          write(iwd)(((dist (ind(i,j,k)),i=i1,i2m1),j=j1,j2m1), &
               k=k1,k2m1)
          write(iwd)(((mnpar(ind(i,j,k)),i=i1,i2m1),j=j1,j2m1), &
               k=k1,k2m1)
       end do
    else
!       ecriture un seul domaine par fichier
!
       l=l0
       if(l.le.9 ) then
          write(nomfich,'("fdist_",i1," ")')l
       else if(l.le.99) then
          write(nomfich,'("fdist_",i2)')l
       else
          write(imp,'("!!!at_ecrdist: plus de 99 domaines. Non prevu")')
          stop
       end if
!
       write(imp,'("===>at_ecrdist: ecriture distance domaine",i2,"   fichier=",a8)')l,nomfich
!
       open(iwd,file=nomfich,form='unformatted',err=50)
!
       n0=npc(l)
       i1=ii1(l)
       i2=ii2(l)
       j1=jj1(l)
       j2=jj2(l)
       k1=kk1(l)
       k2=kk2(l)
       i2m1=i2-1
       j2m1=j2-1
       k2m1=k2-1
       nid = id2(l)-id1(l)+1
       njd = jd2(l)-jd1(l)+1
       nijd = nid*njd
!
       write(iwd)(((dist (ind(i,j,k)),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
       write(iwd)(((mnpar(ind(i,j,k)),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
    end if
!
    close(iwd)
    write(imp,'("===>at_ecrdist: fin ecriture fichier=",a8)')nomfich
!
    if(l.eq.1) then
!       ecriture fichier auxiliaire des donnees necessaires a la relecture
!
       open(iwd,file='fdist-aux',form='formatted',err=60)
       write(imp,'("===>at_ecrdist: ecriture fichier= fdist-aux")')
       write(imp,'(16x,"ip12=",i8)') ip12
       do l=1,lzx
          write(iwd,'(i3,i8,6i5)') &
               l,npc(l),ii1(l),ii2(l),jj1(l),jj2(l),kk1(l),kk2(l)
          write(iwd,'(i3,i8,6i5)') &
               l,npc(l),id1(l),id2(l),jd1(l),jd2(l),kd1(l),kd2(l)
       end do
       close(iwd)
    end if
    return
!
50  continue
    write(imp,'("!!!at_ecrdist: erreur ouverture fichier=",a8)') &
         nomfich
60  continue
    write(imp,'("!!!at_ecrdist: erreur ouverture fichier=fdist-aux")')
    stop
!
  end subroutine at_ecrdist


end module mod_at_ecrdist
