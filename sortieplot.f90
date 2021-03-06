module mod_sortieplot
  implicit none
contains
  subroutine sortieplot(x,y,z,t,ps,cson)
!
!***********************************************************************
!
!_DA  DATE_C :  novembre 2001 -- AUTEUR : SINUMEF Eric Goncalves
!
!     ACT
!_A   Sortie tecplot pour calculs instationnaires.
!_A   densite, vitesse, pression.
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
    use mod_mpi
    implicit none
    integer          ::    i,  i1,i1m1,  i2,i2m1
    integer          ::    j,  j1,j1m1,  j2,j2m1
    integer          ::    k,  k1,k1m1,  k2,k2m1
    integer          ::    l,   m,   n, n0c, nid
    integer          :: nijd, njd,  ll
    double precision ::   cson(ip11),           e,    ps(ip11),          qq,t(ip11,ip60)
    double precision ::            u,           v,           w,     x(ip21),         xcc
    double precision ::          xme,     y(ip21),         ycc,     z(ip21),         zcc
!
!-----------------------------------------------------------------------
!
    character(len=1 ) :: c
    character(len=50) :: file
!
!

!
   do l=1,lzx
     ll=bl_to_bg(l)
     call START_KEEP_ORDER(ll,bg_to_proc)
                     
    write(file,'(A,I0.4)') 'fsec.',bl_to_bg(l)
    open(sec  ,file=trim(file),status="replace")
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
!     double cote
    c=char(34)
!
!    write(sec,'(''TITLE='',a1,a50,a1)')c,titrt1,c
!    write(sec,'(''VARIABLES = '',a1,8(a,a1,'', '',a1),a,a1)') &
!         c,'x',c, c,'y',c, c,'u',c, c,'v',c, c,'rho',c, c,'Pstat',c, c,'Mach',c
!          c,'dist',c, c,'kk',c, c,'w',c, c,'mu',c, c,'mut',c, c,'mutsmu',c
!    write(sec,'("ZONE F=POINT, I=",i3," J=",i3)')j2m1,i2m1
!
    do k=k1,k2m1
!         do j=j1,j2m1
       do i=i1,i2m1
          do j=j1,j2m1
             n=indc(i,j,k)
             m=n-n0c
!
             xcc=(x(n)     +x(n     +1)+x(n     +nid)+x(n     +nid+1) &
                  +x(n+nijd)+x(n+nijd+1)+x(n+nijd+nid)+x(n+nijd+nid+1))*0.125
             ycc=(y(n)     +y(n     +1)+y(n     +nid)+y(n     +nid+1) &
                  +y(n+nijd)+y(n+nijd+1)+y(n+nijd+nid)+y(n+nijd+nid+1))*0.125
             zcc=(z(n)     +z(n     +1)+z(n     +nid)+z(n     +nid+1) &
                  +z(n+nijd)+z(n+nijd+1)+z(n+nijd+nid)+z(n+nijd+nid+1))*0.125
!
             u=t(n,2)/t(n,1)
             v=t(n,3)/t(n,1)
             w=t(n,4)/t(n,1)
             e=t(n,5)/t(n,1)
             qq=u*u+v*v+w*w
             xme=sqrt(qq)/cson(n)

             write(sec,'(7(1pe15.6))') &
                  xcc,ycc,u,v,t(n,1),ps(n),xme
!             xcc,ycc,u,v,t(n,1),ps(n),dist(n),t(n,6),t(n,7),mu(n),mut(n),mut(n)/mu(n)
          enddo
     write(sec,*) ""
       enddo
    enddo
    close(sec)
    
     call END_KEEP_ORDER(ll,bg_to_proc)
   enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine sortieplot
end module mod_sortieplot
