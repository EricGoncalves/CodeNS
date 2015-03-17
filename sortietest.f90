module mod_sortietest
  implicit none
contains
  subroutine sortietest( &
       icycle,ncycl,idcyc, &
       vdual2,dist,vol,mut,mu, &
       x,y,z,l,t,ps,temp)
!
!***********************************************************************
!
!_DA  DATE_C :  mars 2008 -- AUTEUR : Eric Goncalves / LEGI
!
!     ACT
!_A   Sortie pour calculs stationnaires et instationnaires.
!_A   Sigma entree / grandeurs moyennees en temps / pression au plancher
!
!     VAL
!_L    titrt1     : com char             ; titre du calcul
!_L    c          :     char             ; caractere "
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use sortiefichier
    use chainecarac
    use proprieteflu
    use definition
    use schemanum
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2m1
    integer          :: icycle, idcyc,  indc,     j,    j1
    integer          ::   j1m1,    j2,  j2m1,  jmid,     k
    integer          ::     k1,  k1m1,    k2,  k2m1,     l
    integer          ::      m,     n,   n0c, ncycl,   nid
    integer          ::   nijd,   njd,  nmid
    double precision ::   dist,    mu,   mut,    ps,  qinf
    double precision ::      t,  temp,vdual2,   vol,     x
    double precision ::    xcc,     y,   ycc,     z,   zcc
!
!-----------------------------------------------------------------------
!
    character(len=1 ) :: c
!
    dimension vdual2(ip11,ip60),t(ip11,ip60)
    dimension x(ip21),y(ip21),z(ip21)
    dimension ps(ip11),vol(ip11),temp(ip11)
    dimension mu(ip12),mut(ip12),dist(ip12)
!
    indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
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
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
!     double cote
    c=char(34)
!
    qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
!-----initialisation des grandeurs--------------------------------
!
    if(icycle.eq.1) then
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=indc(i,j,k)
                m=n-n0c
                vdual2(n,1)=0.
                vdual2(n,2)=0.
                vdual2(n,3)=0.
                vdual2(n,4)=0.
                vdual2(n,5)=0.
             enddo
          enddo
       enddo
!
!       fichier de sortie du sigma entree
       write(out,'(''TITLE='',a1,a80,a1)')c,titrt1,c
       write(out,'(''VARIABLES = '',a1,4(a,a1,'', '',a1),a,a1)') &
            c,'ite',c, c,'sigmae',c, c,'volvap',c, c,'vit',c
       write(out,'("ZONE F=POINT, I=",i3," J=",i5)')i1,ncycl
!
!       fichier de sortie de la pression
       write(sor2,'(''TITLE='',a1,a80,a1)')c,titrt1,c
       write(sor2,'(''VARIABLES = '',a1,3(a,a1,'', '',a1),a,a1)') &
            c,'x',c, c,'y',c, c,'ps',c
       write(sor2,'("ZONE F=POINT, I=",i3," J=",i3)')i2m1,j2m1
    endif
!
!-----Calculs des moyennes temporelles---------------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          do i=i1,i2m1
             n=indc(i,j,k)
             m=n-n0c
             vdual2(n,1)=vdual2(n,1)+t(n,2)/t(n,1)
             vdual2(n,2)=vdual2(n,2)+t(n,1)
             vdual2(n,3)=vdual2(n,3)+ps(n)
             vdual2(n,4)=vdual2(n,4)+temp(n)
             vdual2(n,5)=vdual2(n,5)+mut(n)/mu(n)
          enddo
       enddo
    enddo
!
!-----Ecriture a la derniere iteration--------------------------------
!
    if(icycle.eq.ncycl) then
!
!       fichier de sortie des grandeurs moyennees en temps
       write(sor1,'(''TITLE='',a1,a80,a1)')c,titrt1,c
       write(sor1,'(''VARIABLES = '',a1,9(a,a1,'', '',a1),a,a1)') &
            c,'x',c, c,'y',c, c,'dist',c, c,'umoy',c, c,'rhomoy',c, &
            c,'Pmoy',c, c,'tempmoy',c, c,'mutsmumoy',c
       write(sor1,'("ZONE F=POINT, I=",i3," J=",i3)')j2m1,i2m1
!
       do k=k1,k2m1
          do i=i1,i2m1
             do j=j1,j2m1
                n=indc(i,j,k)
                m=n-n0c
                xcc=(x(n)     +x(n     +1)+x(n     +nid)+x(n     +nid+1) &
                     +x(n+nijd)+x(n+nijd+1)+x(n+nijd+nid)+x(n+nijd+nid+1))*0.125
                ycc=(y(n)     +y(n     +1)+y(n     +nid)+y(n     +nid+1) &
                     +y(n+nijd)+y(n+nijd+1)+y(n+nijd+nid)+y(n+nijd+nid+1))*0.125
                zcc=(z(n)     +z(n     +1)+z(n     +nid)+z(n     +nid+1) &
                     +z(n+nijd)+z(n+nijd+1)+z(n+nijd+nid)+z(n+nijd+nid+1))*0.125
!
                write(sor1,'(8(1pe15.6))') &
                     xcc,ycc,dist(n),vdual2(n,1)/ncycl,vdual2(n,2)/ncycl, &
                     vdual2(n,3)/ncycl,vdual2(n,4)/ncycl,vdual2(n,5)/ncycl
             enddo
          enddo
       enddo
    endif
!
!------sortie de la pression au plancher
!
    do k=k1,k2m1
       do i=i1,i2m1
!         do j=j1,j2m1
          j=j1
          n=indc(i,j,k)
          m=n-n0c
          xcc=(x(n)     +x(n     +1)+x(n     +nid)+x(n     +nid+1) &
               +x(n+nijd)+x(n+nijd+1)+x(n+nijd+nid)+x(n+nijd+nid+1))*0.125
          ycc=(y(n)     +y(n     +1)+y(n     +nid)+y(n     +nid+1) &
               +y(n+nijd)+y(n+nijd+1)+y(n+nijd+nid)+y(n+nijd+nid+1))*0.125
          zcc=(z(n)     +z(n     +1)+z(n     +nid)+z(n     +nid+1) &
               +z(n+nijd)+z(n+nijd+1)+z(n+nijd+nid)+z(n+nijd+nid+1))*0.125
!
          write(sor2,'(3(1pe15.6))') &
               xcc,ycc,ps(n)
       enddo
    enddo
!
    return
  end subroutine sortietest
end module mod_sortietest
