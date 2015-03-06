      subroutine svdual(l,vdual,vdual1,vdual2)
!
!***********************************************************************
!
!_DA  DATE_C : janvier 2002 - Eric Goncalves - SINUMEF
!
!     ACT
!_A    sauvegarde  v  - ordre 3
!
!----------------------------------------------------------------
! v             champ a l'instant n+alpha (RGK)
! vdual         champ a l'instant n
! vdual1        champ a l'instant n-1
! vdual2        champ a l'instant n-2
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use schemanum
      use maillage
!
!-----------------------------------------------------------------------
!
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60)
!
      ind(i,j,k)= n0 + 1 + (i-id1(l)) + (j-jd1(l))*nid + (k-kd1(l))*nijd
!
!            Fichier sauvegarde
!
      open(UNIT=98,FILE='facdual',FORM='unformatted',STATUS='unknown')
!
      n0   = npc(l)
      i1   = ii1(l)
      i2   = ii2(l)
      i2m1 = i2 - 1
      j1   = jj1(l)
      j2   = jj2(l)
      j2m1 = j2 - 1
      k1   = kk1(l)
      k2   = kk2(l)
      k2m1 = k2 - 1
      nid  = id2(l) - id1(l) + 1
      njd  = jd2(l) - jd1(l) + 1
      nijd = nid * njd
!
      write(98) &
         (((vdual(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!
      write(98) &
         (((vdual1(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual1(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual1(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual1(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual1(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
        (((vdual1(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
        (((vdual1(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!
      write(98) &
         (((vdual2(ind(i,j,k),1),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual2(ind(i,j,k),2),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual2(ind(i,j,k),3),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual2(ind(i,j,k),4),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
         (((vdual2(ind(i,j,k),5),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
        (((vdual2(ind(i,j,k),6),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
      write(98) &
        (((vdual2(ind(i,j,k),7),i=i1,i2m1),j=j1,j2m1),k=k1,k2m1)
!
      close(98)
!
      return
      end
