      subroutine tcmd_dfpmdsd( &
                 mot,imot,nmot, &
                 ldom,ldomd, &
                 lgr,lgrd)
!
!***********************************************************************
!
!     ACT
!_A    Traduction des mots lus en donnees neccessaires a l'action dfpmdsd.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use schemanum
      use maillage
      use kcle
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension lgr(nobj)
!
      do icmt=1,32
       comment(icmt:icmt)=' '
      enddo
!
      do l=1,lt
       if(kki2(l).eq.2) kki2(l)=3
       if(kki4(l).eq.2) kki4(l)=3
      enddo
!
      if(nmot.eq.2)then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      nm=2
!
      nm=nm+1
      if((imot(nm).eq.4).and.(mot(nm).eq.'ldom')) then
        nm=nm+1
        if(nmot.lt.nm) then
          comment=cm
          call synterr(mot,imot,nmot,comment)
        else
          call vallent(mot,imot,nm,ldom,ldomd,lzx,klzx)
          nm=nm+1
          if((imot(nm).eq.3).and.(mot(nm).eq.'lgr')) then
            nm=nm+1
            if(nmot.lt.nm) then
              comment=cm
              call synterr(mot,imot,nmot,comment)
            else
              call vallent(mot,imot,nm,lgr,lgrd,lgx,klgx)
!
               do while(nm.lt.nmot)
              nm=nm+1
              if((imot(nm).eq.3).and.(mot(nm).eq.'ki2')) then
                nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
                  call valreel(mot,imot,nm,rree,krree)
                  do nl=1,ldomd
                   l=ldom(nl)
                   do ng=1,lgrd
                    img=lgr(ng)
                    lm=l+(img-1)*lz
                    ki2(lm)=rree
                    kki2(lm)=krree
                   enddo
                  enddo
                endif
              elseif((imot(nm).eq.3).and.(mot(nm).eq.'ki4')) then
                nm=nm+1
                if(nmot.lt.nm) then
                  comment=cr
                  call synterr(mot,imot,nmot,comment)
                else
                  call valreel(mot,imot,nm,rree,krree)
                  do nl=1,ldomd
                   l=ldom(nl)
                   do ng=1,lgrd
                    img=lgr(ng)
                    lm=l+(img-1)*lz
                    ki4(lm)=rree
                    kki4(lm)=krree
                   enddo
                  enddo
                endif
              else if(imot(nm).eq.0) then
                comment=cs
                call synterr(mot,imot,nm,comment)
              else
                comment=cb
                call synterr(mot,imot,nm,comment)
              end if
             enddo
            endif
          else
            comment=cb
            call synterr(mot,imot,nm,comment)
          endif
        endif
      else
        comment=cb
        call synterr(mot,imot,nm,comment)
      endif
!
      return
      end
