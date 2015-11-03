module mod_residu
  implicit none
contains
  subroutine residu( &
       img, &
       l, &
       u,v,dt, &
       res1,res2,res3,res4,res5,res6,res7, &
       icyc, &
       dumy1,dumy2,dumax, &
       idumx,jdumx,kdumx)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des residus ((v(n+1)-v(n))/dt) en tout point d'un domaine.
!_A    Calcul des residus moyens en normes L1 et L2 et maximaux avec
!_A    le point ou le maximal est atteint et ecriture de ces donnees
!_A    sur un fichier d'exploitation visuelle et sur un fichier
!_A    d'exploitaion informatique (norme L2 et max).
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kimp       : arg int              ; niveau de sortie sur unite logi imp
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    dt         : arg real(ip11      ) ; pas de temps
!_I    it         : arg int              ; cycle courant du calcul
!_I    out        : com int              ; unite logiq, moyennes des residus
!_I    reelmn     : com real             ; nombre reel petit
!_I    imp        : com int              ; unite logiq, sorties de controle
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
!
!     OUT
!_O    res1       : arg real(ip00      ) ; residu en un point, (instant n+1
!_O                                        - instant n) / pas de temps pour
!_O                                        l'equation 1
!_O    res2       : arg real(ip00      ) ; residu en un point, (instant n+1
!_O                                        - instant n) / pas de temps pour
!_O                                        l'equation 2
!_O    res3       : arg real(ip00      ) ; residu en un point, (instant n+1
!_O                                        - instant n) / pas de temps pour
!_O                                        l'equation 3
!_O    res4       : arg real(ip00      ) ; residu en un point, (instant n+1
!_O                                        - instant n) / pas de temps pour
!_O                                        l'equation 4
!_O    res5       : arg real(ip00      ) ; residu en un point, (instant n+1
!_O                                        - instant n) / pas de temps pour
!_O                                        l'equation 5
!_O    dumy2      : arg real(neqt      ) ; residus quadratiques moyens sur le domaine
!_O    dumax      : arg real(neqt      ) ; residus maximaux sur le domaine
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use sortiefichier
    use constantes
    implicit none
    integer          ::     i,   i1,   i2, i2m1, icyc
    integer          :: idumx,  img,    j,   j1,   j2
    integer          ::  j2m1,jdumx,    k,   k1,   k2
    integer          ::  k2m1,kdumx,    l,    m,    n
    integer          ::    n0,   ni,  nid, nijd,   nj
    integer          ::   njd,   nk, nmax
    double precision ::  coef,   dt,dumax,dumy1,dumy2
    double precision ::  res1, res2, res3, res4, res5
    double precision ::  res6, res7,    s,    u,    v
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
!
    dimension idumx(neqt),jdumx(neqt),kdumx(neqt), &
         dumy1(neqt),dumy2(neqt),dumax(neqt)
    dimension u(ip11,ip60),v(ip11,ip60)
    dimension dt(ip11)
    dimension res1(ip00),res2(ip00),res3(ip00),res4(ip00), &
         res5(ip00),res6(ip00),res7(ip00)
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
    nijd = nid*njd
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    ni = i2-i1
    nj = j2-j1
    nk = k2-k1
!
    nmax=ni*nj*nk
    coef=1./real(nmax)
!
    if(equat(6:7).eq.'ke') then
       do m=1,neqt
          dumy1(m)=0.
          dumy2(m)=0.
          dumax(m)=0.
          idumx(m)=0
          jdumx(m)=0
          kdumx(m)=0
       enddo
    else
       do m=1,5
          dumy1(m)=0.
          dumy2(m)=0.
          dumax(m)=0.
          idumx(m)=0
          jdumx(m)=0
          kdumx(m)=0
       enddo
    end if
!
    do k=k1,k2m1
       do j=j1,j2m1
          do i=i1,i2m1
             n=ind(i,j,k)
             m=n-n0
             res1(m)=(v(n,1)-u(n,1))/dt(n)
             dumy1(1)=dumy1(1)+abs(res1(m))
             dumy2(1)=dumy2(1)+res1(m)*res1(m)
             dumax(1)=max(dumax(1),abs(res1(m)))
             s=0.5*(sign(1.D0,abs(res1(m))-dumax(1))+1.)
             res1(m)=log10(max(abs(res1(m)),reelmn))
             idumx(1)=nint((1.-s)*idumx(1)+s*i)
             jdumx(1)=nint((1.-s)*jdumx(1)+s*j)
             kdumx(1)=nint((1.-s)*kdumx(1)+s*k)
!
             res2(m)=(v(n,2)-u(n,2))/dt(n)
             dumy1(2)=dumy1(2)+abs(res2(m))
             dumy2(2)=dumy2(2)+res2(m)*res2(m)
             dumax(2)=max(dumax(2),abs(res2(m)))
             s=0.5*(sign(1.D0,abs(res2(m))-dumax(2))+1.)
             res2(m)=log10(max(abs(res2(m)),reelmn))
             idumx(2)=nint((1.-s)*idumx(2)+s*i)
             jdumx(2)=nint((1.-s)*jdumx(2)+s*j)
             kdumx(2)=nint((1.-s)*kdumx(2)+s*k)
!
             res3(m)=(v(n,3)-u(n,3))/dt(n)
             dumy1(3)=dumy1(3)+abs(res3(m))
             dumy2(3)=dumy2(3)+res3(m)*res3(m)
             dumax(3)=max(dumax(3),abs(res3(m)))
             s=.5*(sign(1.D0,abs(res3(m))-dumax(3))+1.)
             res3(m)=log10(max(abs(res3(m)),reelmn))
             idumx(3)=nint((1.-s)*idumx(3)+s*i)
             jdumx(3)=nint((1.-s)*jdumx(3)+s*j)
             kdumx(3)=nint((1.-s)*kdumx(3)+s*k)
!
             res4(m)=(v(n,4)-u(n,4))/dt(n)
             dumy1(4)=dumy1(4)+abs(res4(m))
             dumy2(4)=dumy2(4)+res4(m)*res4(m)
             dumax(4)=max(dumax(4),abs(res4(m)))
             s=.5*(sign(1.D0,abs(res4(m))-dumax(4))+1.)
             res4(m)=log10(max(abs(res4(m)),reelmn))
             idumx(4)=nint((1.-s)*idumx(4)+s*i)
             jdumx(4)=nint((1.-s)*jdumx(4)+s*j)
             kdumx(4)=nint((1.-s)*kdumx(4)+s*k)
!
             res5(m)=(v(n,5)-u(n,5))/dt(n)
             dumy1(5)=dumy1(5)+abs(res5(m))
             dumy2(5)=dumy2(5)+res5(m)*res5(m)
             dumax(5)=max(dumax(5),abs(res5(m)))
             s=.5*(sign(1.D0,abs(res5(m))-dumax(5))+1.)
             res5(m)=log10(max(abs(res5(m)),reelmn))
             idumx(5)=nint((1.-s)*idumx(5)+s*i)
             jdumx(5)=nint((1.-s)*jdumx(5)+s*j)
             kdumx(5)=nint((1.-s)*kdumx(5)+s*k)
          enddo
       enddo
    enddo
!
    if (equat(6:7).eq.'ke') then
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=ind(i,j,k)
                m=n-n0
                res6(m)=(v(n,6)-u(n,6))/dt(n)
                dumy1(6)=dumy1(6)+abs(res6(m))
                dumy2(6)=dumy2(6)+res6(m)*res6(m)
                dumax(6)=max(dumax(6),abs(res6(m)))
                s=.5*(sign(1.D0,abs(res6(m))-dumax(6))+1.)
                res6(m)=log10(max(abs(res6(m)),reelmn))
                idumx(6)=nint((1.-s)*idumx(6)+s*i)
                jdumx(6)=nint((1.-s)*jdumx(6)+s*j)
                kdumx(6)=nint((1.-s)*kdumx(6)+s*k)
!
                res7(m)=(v(n,7)-u(n,7))/dt(n)
                dumy1(7)=dumy1(7)+abs(res7(m))
                dumy2(7)=dumy2(7)+res7(m)*res7(m)
                dumax(7)=max(dumax(7),abs(res7(m)))
                s=.5*(sign(1.D0,abs(res7(m))-dumax(7))+1.)
                res7(m)=log10(max(abs(res7(m)),reelmn))
                idumx(7)=nint((1.-s)*idumx(7)+s*i)
                jdumx(7)=nint((1.-s)*jdumx(7)+s*j)
                kdumx(7)=nint((1.-s)*kdumx(7)+s*k)
             enddo
          enddo
       enddo
    end if
!
    form='(10x,2hdu,i1,19h/dt : moyenne l1 = ,1pe10.3,' &
         //'16h - moyenne l2 = ,1pe10.3,13h - maximum = ,' &
         //'1pe10.3,20h atteint au point i=,i3,3h j=,' &
         //'i3,3h k=,i3)'
!
    if(equat(6:7).eq.'ke') then
       do m=1,neqt
          dumy1(m)=coef*dumy1(m)
          dumy2(m)=sqrt(coef*dumy2(m))
          if(kimp.ge.1) then
             write(imp,form) m,dumy1(m),dumy2(m),dumax(m), &
                  idumx(m),jdumx(m),kdumx(m)
          endif
       enddo
       if(kimp.ge.1) then
          form='(/)'
          write(imp,form)
       endif
    else
       do m=1,5
          dumy1(m)=coef*dumy1(m)
          dumy2(m)=sqrt(coef*dumy2(m))
          if(kimp.ge.1) then
             write(imp,form) m,dumy1(m),dumy2(m),dumax(m), &
                  idumx(m),jdumx(m),kdumx(m)
          endif
       enddo
       if(kimp.ge.1) then
          form='(/)'
          write(imp,form)
       endif

    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
  end subroutine residu
end module mod_residu
