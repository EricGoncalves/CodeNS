      subroutine tcmd_dfst( &
                 mot,imot,nmot, &
                 nst)
!
!***********************************************************************
!
!     traduction des mots lus en donnees neccessaires a
!     l'action dfst
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use chainecarac
      use kcle
      use definition
implicit none
integer :: imot
integer :: nmot
integer :: nst
integer :: icmt
integer :: kval
integer :: lst
integer :: nm
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  comment
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      do icmt=1,32
       comment(icmt:icmt)=' '
      enddo
      kval=0
!
      do nst=1,nsta
       do lst=1,lsta
        if (kvarst(nst,lst).eq.2) kvarst(nst,lst)=3
       enddo
      enddo
!
      if(nmot.eq.2) then
        comment=cb
        call synterr(mot,imot,2,comment)
      endif
!
      if(nmot.gt.2) then
        nm=2
!
        nm=nm+1
          call valenti(mot,imot,nm,nst,kval)
!
        nm=nm+1
        if((imot(nm).eq.5).and.(mot(nm).eq.'ktype')) then
          nm=nm+1
          if(nmot.lt.nm) then
            comment=ci
            call synterr(mot,imot,nmot,comment)
          endif
        else
          comment=cb
          call synterr(mot,imot,nm,comment)
        endif
!
        do while(nm.lt.nmot)
        nm=nm+1
          if((imot(nm).eq.3).and.(mot(nm).eq.'roi')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,1),kvarst(nst,1))
!
          else if((imot(nm).eq.2).and.(mot(nm).eq.'ai')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,2),kvarst(nst,2))
!
          else if((imot(nm).eq.2).and.(mot(nm).eq.'ti')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,3),kvarst(nst,3))
!
          else if((imot(nm).eq.4).and.(mot(nm).eq.'long')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,4),kvarst(nst,4))
!
          else if((imot(nm).eq.4).and.(mot(nm).eq.'mach')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,5),kvarst(nst,5))
!
          else if((imot(nm).eq.5).and.(mot(nm).eq.'alpha')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,6),kvarst(nst,6))
!
          else if((imot(nm).eq.4).and.(mot(nm).eq.'beta')) then
              nm=nm+1
            call valreel(mot,imot,nm,varst(nst,7),kvarst(nst,7))
!
          else if(imot(nm).eq.0) then
            comment=cs
            call synterr(mot,imot,nm,comment)
!
          else
            comment=cb
            call synterr(mot,imot,nm,comment)
          end if
!
         enddo
        endif
!
      return
      end
