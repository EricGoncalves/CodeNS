module mod_rscpsv
  implicit none
contains
  subroutine rscpsv( &
       img, &
       u,v,dt, &
       res1xx,res2yy,res3zz,res4,res5,res6,res7, &
       tn8, &
       icyc,ncyc,ncycle, &
       x,y,z,utau)
!
!***********************************************************************
!
!     ACT
!_A    Determination des residus ((u(n+1)-u(n))/dt) en tous points de tous les
!_A    domaines et ecriture sur un fichier d'exploitation ainsi que les
!_A    coordonnees des centres des mailles ou sont connus les residus.
!_A    Calcul des residus moyens en normes L2 et maximaux sur tous les domaines
!_A    et ecriture de ces donnees sur un fichier d'exploitaion.
!
!     INP
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    dt         : arg real(ip11      ) ; pas de temps
!_I    icyc       : arg int              ; cycle courant du calcul
!_I    ncyc       : arg int              ; cycle courant de l'execution
!_I    ncycle     : arg int              ; nbr tot de cycles de l'execution courante
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    out        : com int              ; unite logiq, moyennes des residus
!_I    kdgc       : com int              ; unite logiq, coordonnees des centres
!_I    kres       : com int              ; unite logiq, residus en chaque pt
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    lzx        : com int              ; nbr total de domaines
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!
!     LOC
!_L    res1xx     : arg real(ip00      ) ; residu pour equat 1 ou coordonnee x
!_L                                        au centre de la maille
!_L    res2yy     : arg real(ip00      ) ; residu pour equat 2 ou coordonnee y
!_L                                        au centre de la maille
!_L    res3zz     : arg real(ip00      ) ; residu pour equat 3 ou coordonnee z
!_L                                        au centre de la maille
!_L    res4       : arg real(ip00      ) ; tableau de travail
!_L    res5       : arg real(ip00      ) ; tableau de travail
!_L    tn8        : arg real(ip00      ) ; tableau de travail
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use sortiefichier
    use chainecarac
    use mod_cvccg
    use mod_writdg
    use mod_residu
    use mod_mpi
    implicit none
    integer          ::   icyc,  imax,   img,  imin,  jmax
    integer          ::   jmin,  kmax,  kmin,     l,    lm
    integer          ::      m,    n0,  ncyc,ncycle,    ni
    integer          ::    nid,  nijd,    nj,   njd,    nk
    integer          ::   nmax,  npts,ll
    double precision ::   dt(ip11),res1xx(ip00),res2yy(ip00),res3zz(ip00),res4(ip00)
    double precision ::   res5(ip00),res6(ip00),res7(ip00),tn8(ip00),u(ip11,ip60)
    double precision ::   utau(ip42),v(ip11,ip60),x(ip21),y(ip21),z(ip21)
    integer         ,allocatable :: idumx(:),jdumx(:),kdumx(:)
    double precision,allocatable ::  dumax(:),dumaxg(:), dumy1(:), dumy2(:),dumy2g(:)
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    ALLOCATE(idumx(neqt),jdumx(neqt),kdumx(neqt), &
         dumy1(neqt),dumy2(neqt),dumax(neqt),dumy2g(neqt),dumaxg(neqt))

    dumy1=0.
    dumy2=0.
    dumax=0.
    dumy2g=0.
    dumaxg=0.
!
    if (kimp.ge.1.and.rank==0) then
       form='(/1x,2h--,1x,i7,13h ieme cycle :,3x,' &
            //'20hnb total de cycles =,i6)'
       write(imp,form) icyc,ncycle
    endif
!
    npts=0
    do l=1,lzx
       lm=l+(img-1)*lz
       ll=bl_to_bg(l)
       call START_KEEP_ORDER(ll,bg_to_proc)
       if (kimp.ge.1) then
          form='(/10x,10hzone no : ,i5,5x,12hgrille no : ,i3/)'
          write(imp,form) ll,img
       endif
!
       call residu( &
            img, &
            lm, &
            u,v,dt, &
            res1xx,res2yy,res3zz,res4,res5,res6,res7, &
            icyc, &
            dumy1,dumy2,dumax, &
            idumx,jdumx,kdumx)
       call END_KEEP_ORDER(ll,bg_to_proc)
!
       temp_array(:,1)=dumy1
       temp_array(:,2)=dumy2
       temp_array(:,3)=dumax

       if(img.eq.1) then
          if(ncyc.eq.ncycle) then
             imin=ii1(l)
             imax=ii2(l)-1
             jmin=jj1(l)
             jmax=jj2(l)-1
             kmin=kk1(l)
             kmax=kk2(l)-1
!
!         call writda( &
!                 l,kres,' res   ',utau, &
!                 imin,imax,jmin,jmax,kmin,kmax, &
!                 res1xx,res2yy,res3zz,res4,res5,res6,res7,1,tn8)
!
             nid = id2(l)-id1(l)+1
             njd = jd2(l)-jd1(l)+1
             nijd= nid*njd
             n0  = npc(l)
!
             call cvccg( &
                  l, &
                  x,y,z, &
                  res1xx,res2yy,res3zz)
!             call writdg( &
!                  l,kdgc, &
!                  imin,imax,jmin,jmax,kmin,kmax, &
!                  res1xx,res2yy,res3zz)
          endif
!
          ni = ii2(l)-ii1(l)
          nj = jj2(l)-jj1(l)
          nk = kk2(l)-kk1(l)
!
          nmax=ni*nj*nk
          npts=npts+nmax
!
          if(equat(6:7).eq.'ke') then
             do m=1,7
                dumy2g(m)=dumy2g(m)+nmax*dumy2(m)**2
                dumaxg(m)=max(dumaxg(m),abs(dumax(m)))
             enddo
          else
             do m=1,5
                dumy2g(m)=dumy2g(m)+nmax*dumy2(m)**2
                dumaxg(m)=max(dumaxg(m),abs(dumax(m)))
             enddo
          endif
       endif !gin test img
!
    enddo
    call sum_mpi(dumy2g)
    call max_mpi(dumaxg)
    call sum_mpi(npts)
!
    if(img.eq.1) then
       if(equat(6:7).eq.'ke') then
          do m=1,7
             dumy2g(m)=sqrt(dumy2g(m)/npts)
          enddo
        form='(A,i6,1x,7e11.4)'
        if (rank==0) write(imp,form) "Stationnarite L2  ",icyc,(dumy2g(m),m=1,7)
        if (rank==0) write(imp,form) "Stationnarite Linf",icyc,(dumaxg(m),m=1,7)
       else
          do m=1,5
             dumy2g(m)=sqrt(dumy2g(m)/npts)
          enddo
        form='(i6,1x,10e11.4)'
        if (rank==0) write(out,form) icyc,(dumy2g(m),m=1,5),(dumaxg(m),m=1,5)
       endif
    endif

    DEALLOCATE(idumx,jdumx,kdumx,dumy1,dumy2,dumax,dumy2g,dumaxg)

    return
  end subroutine rscpsv
end module mod_rscpsv
