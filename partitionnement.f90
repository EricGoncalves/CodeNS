module mod_partitionnement
  use tools
  implicit none
contains
  subroutine partitionnement(x,y,z,mot,imot,nmot,ncbd,mnc,ncin,bceqt,exs1,exs2)
    use mod_valenti
    use boundary
    use chainecarac
    use para_var
    use schemanum
    use kcle
    use modeleturb
    use sortiefichier
    use mod_crbds
    use mod_c_crbds
    use mod_inbdc
    use mod_c_inbdc
    use mod_inbdb
    use mod_c_inbdb
    use mod_crdms
    use mod_c_crdms
    use mod_mpi
    !
    !***********************************************************************
    !
    !     act
    !_a    realisation du partitionnement
    !      ecrit par Alexandre Poux
    !
    !***********************************************************************
    !-----parameters figes--------------------------------------------------
    !
    implicit none
    integer             :: nblocks,nxyza,i,j,k,xyz,xs,ys,zs,nmin,nmax,nmin1,nmax1
    integer             :: imot(nmx),nmot,fr,imax,imin,jmax,jmin,kmax,kmin
    integer             :: l2,mfbe,nid,njd,nijd,sblock,l3
    integer             :: fr2,i2,j2,k2,l4,fr3,fri2
    double precision    :: exs1,exs2,vbc(ista*lsta),sub_bc1(2,6)
    double precision,allocatable :: x(:),y(:),z(:)
    integer,allocatable :: nblock2(:),nblockd(:,:),ni(:,:,:,:),nj(:,:,:,:),nk(:,:,:,:)
    integer,allocatable :: new2old_b(:),num_cf2(:,:,:)
    integer,allocatable :: ni1(:,:,:),nj1(:,:,:),nk1(:,:,:),ncbd(:)

    integer             :: save_lt,l,verbosity,l1,ll2,ll3
    integer             :: save_lzx,save_mdimtbx,save_mdimubx
    integer             :: save_mtb,save_mtt,save_ndimctbx
    integer             :: save_ndimubx
    integer             :: val(3),fr4,fri,save_mtbx
    double precision,allocatable :: save_x(:),save_y(:),save_z(:),save_bc(:,:),save_bceqt(:,:),save_vbc(:)
    integer,allocatable :: save_ii1(:),save_jj1(:),save_kk1(:)
    integer,allocatable :: save_ii2(:),save_jj2(:),save_kk2(:)
    integer,allocatable :: save_id1(:),save_jd1(:),save_kd1(:)
    integer,allocatable :: save_id2(:),save_jd2(:),save_kd2(:)
    integer,allocatable :: save_npn(:),save_ndlb(:)
    integer,allocatable :: save_iminb(:),save_jminb(:),save_kminb(:)
    integer,allocatable :: save_imaxb(:),save_jmaxb(:),save_kmaxb(:)
    integer,allocatable :: save_nnn(:),save_nnc(:),save_nnfb(:)
    integer,allocatable :: save_npc(:),save_npfb(:),save_nfei(:)
    integer,allocatable :: save_mpb(:),save_ncbd(:),save_mmb(:)
    integer,allocatable :: save_ndcc(:),save_nbdc(:),save_nfbc(:)
    integer,allocatable :: save_mdnc(:),save_mper(:),save_mpc(:)
    integer,allocatable :: save_ncin(:),save_mnc(:)
    integer,allocatable :: save_bcg_to_proc(:),save_bcg_to_bcl(:),save_bcg_to_bci(:),save_bcl_to_bcg(:)
    integer,allocatable :: save_bg_to_proc(:),save_bg_to_bl(:),save_bg_to_bi(:),save_bl_to_bg(:)
    integer,allocatable :: ii2g(:),jj2g(:),kk2g(:),nblockdg(:,:),save_bcg_to_bg(:),sub_bc(:,:)
    integer             :: save_ip41,save_num_bcg,save_num_bci,save_num_bcl,save_num_bg
    integer             :: save_num_bl,save_num_bi,proc,orig,dest,xe,ye,ze,fr1,orig1,orig2,nsub
    character(len=2),allocatable :: save_indfl(:)
    character(len=4),allocatable :: save_cl(:)

    character(len=50)::fich
    character(len=32) ::  mot(nmx)
    character(len=2) :: indmf

    integer         ,allocatable ::   mnc(:),ncin(:)
    double precision,allocatable ::  bceqt(:,:)
    logical :: test

    !############################################################################################
    !############################## GET PARAMETERS ##############################################
    !############################################################################################

    verbosity=2 ! from 0 to 3
    nblocks=max(nprocs,num_bg)
if(.true.)then
!if(nblocks/=num_bg) then ! there is some partitionning to do

    !############################################################################################
    !############################# SAVE OLD MESH FOR CHECKING PURPOSE ###########################
    !############################################################################################

    if(verbosity>=3) then
    ! write old grid
      do l=1,lt
         write(fich,'(A,I0.2,A)') "origmesh_",bl_to_bg(l),".dat"
         open(42,file=fich,status="replace")

         do k=kk1(l),kk2(l)
            do j=jj1(l),jj2(l)
               do i=ii1(l),ii2(l)

                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),bl_to_bg(l)

               enddo
               write(42,*) ""
            enddo
         enddo
         close(42)
      enddo

    ! write old boundaries
      do fr=1,mtb

      write(fich,'(A,I0.2,A)') "origbnd_",bcl_to_bcg(fr),".dat"
      open(42,file=fich,status="replace")

         do k=kminb(fr),kmaxb(fr)
            do j=jminb(fr),jmaxb(fr)
               do i=iminb(fr),imaxb(fr)
                  l=ndlb(fr)
                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),bcl_to_bcg(fr)

               enddo
               if(indfl(fr)(1:1)/="i") write(42,*) ""
            enddo
            if(indfl(fr)(1:1)=="i") write(42,*) ""
         enddo
      enddo
      close(42)
    endif


    !############################################################################################
    !################# SAVE ALL MESH AND BOUNDARIES VARIABLES ###################################
    !#################            AND CLEAR EVERYTHING        ###################################
    !#################       WE WILL RECREATE EVERYTHING      ###################################
    !############################################################################################

    ! allocate statments are non necessary with newer fortran norm, but ifort don't support it yet
    allocate(save_x(size(x)))
    allocate(save_y(size(y)))
    allocate(save_z(size(z)))
    allocate(save_ndlb(size(ndlb)))
    allocate(save_nfei(size(nfei)))
    allocate(save_indfl(size(indfl)))
    allocate(save_iminb(size(iminb)))
    allocate(save_imaxb(size(imaxb)))
    allocate(save_jminb(size(jminb)))
    allocate(save_jmaxb(size(jmaxb)))
    allocate(save_kminb(size(kminb)))
    allocate(save_kmaxb(size(kmaxb)))
    allocate(save_mpb(size(mpb)))
    allocate(save_mmb(size(mmb)))
    allocate(save_ncbd(size(ncbd)))
    allocate(save_ii1(size(ii1)))
    allocate(save_jj1(size(jj1)))
    allocate(save_kk1(size(kk1)))
    allocate(save_ii2(size(ii2)))
    allocate(save_jj2(size(jj2)))
    allocate(save_kk2(size(kk2)))
    allocate(save_id1(size(id1)))
    allocate(save_jd1(size(jd1)))
    allocate(save_kd1(size(kd1)))
    allocate(save_id2(size(id2)))
    allocate(save_jd2(size(jd2)))
    allocate(save_kd2(size(kd2)))
    allocate(save_nnn(size(nnn)))
    allocate(save_nnc(size(nnc)))
    allocate(save_nnfb(size(nnfb)))
    allocate(save_npn(size(npn)))
    allocate(save_npc(size(npc)))
    allocate(save_npfb(size(npfb)))
    allocate(save_cl(size(cl)))
    allocate(save_ndcc(size(ndcc)))
    allocate(save_nbdc(size(nbdc)))
    allocate(save_nfbc(size(nfbc)))
    allocate(save_bc(size(bc,1),size(bc,2)))
    allocate(save_mdnc(size(mdnc)))
    allocate(save_mper(size(mper)))
    allocate(save_mpc(size(mpc)))
    allocate(save_ncin(size(ncin)))
    allocate(save_bceqt(size(bceqt,1),size(bceqt,2)))
    allocate(save_mnc(size(mnc)))
    allocate(save_bcg_to_proc(size(bcg_to_proc)))
    allocate(save_bcg_to_bcl(size(bcg_to_bcl)))
    allocate(save_bcg_to_bci(size(bcg_to_bci)))
    allocate(save_bcl_to_bcg(size(bcl_to_bcg)))
    allocate(save_bg_to_proc(size(bg_to_proc)))
    allocate(save_bg_to_bl(size(bg_to_bl)))
    allocate(save_bg_to_bi(size(bg_to_bi)))
    allocate(save_bl_to_bg(size(bl_to_bg)))
    allocate(save_bcg_to_bg(size(bcg_to_bg)))
    allocate(save_vbc(size(vbc)))

    ! Save old split
    save_lzx=lzx
    save_lt=lt

    ! Save old grid
    save_x=x
    save_y=y
    save_z=z
    save_ii1 = ii1
    save_jj1 = jj1
    save_kk1 = kk1
    save_ii2 = ii2
    save_jj2 = jj2
    save_kk2 = kk2
    save_id1 = id1
    save_jd1 = jd1
    save_kd1 = kd1
    save_id2 = id2
    save_jd2 = jd2
    save_kd2 = kd2

    save_nnn  = nnn
    save_nnc  = nnc
    save_nnfb = nnfb

    save_npn = npn
    save_npc = npc
    save_npfb=npfb
    save_ndimubx = ndimubx
    save_ndimctbx=ndimctbx

    ! Save old boundary
    save_ndlb=ndlb
    save_nfei=nfei
    save_indfl=indfl
    save_iminb=iminb
    save_imaxb=imaxb
    save_jminb=jminb
    save_jmaxb=jmaxb
    save_kminb=kminb
    save_kmaxb=kmaxb
    save_mpb=mpb
    save_mmb=mmb
    save_ncbd=ncbd
    save_cl=cl
    save_ndcc=ndcc
    save_nbdc=nbdc
    save_nfbc=nfbc
    save_bc=bc
    save_mdnc=mdnc
    save_mper=mper
    save_mpc=mpc
    save_ncin=ncin
    save_bceqt=bceqt
    save_mnc=mnc

    save_mtbx=mtbx
    save_mtb=mtb
    save_mtt=mtt
    save_mdimubx=mdimubx
    save_mdimtbx=mdimtbx

    save_ip41=ip41
    save_vbc=vbc
    save_num_bcg=num_bcg
    save_num_bci=num_bci
    save_num_bcl=num_bcl
    save_num_bg=num_bg
    save_num_bi=num_bi
    save_num_bl=num_bl

    save_bcg_to_proc=bcg_to_proc
    save_bcg_to_bcl=bcg_to_bcl
    save_bcg_to_bci=bcg_to_bci
    save_bcl_to_bcg=bcl_to_bcg
    save_bg_to_proc=bg_to_proc
    save_bg_to_bl=bg_to_bl
    save_bg_to_bi=bg_to_bi
    save_bl_to_bg=bl_to_bg
    save_bcg_to_bg=bcg_to_bg

    !   reinitialisation

    call reallocate(x,0)
    call reallocate(y,0)
    call reallocate(z,0)
    call reallocate(ndlb,0)
    call reallocate(nfei,0)
    call reallocate(indfl,0)
    call reallocate(iminb,0)
    call reallocate(imaxb,0)
    call reallocate(jminb,0)
    call reallocate(jmaxb,0)
    call reallocate(kminb,0)
    call reallocate(kmaxb,0)
    call reallocate(mpb,0)
    call reallocate(mmb,0)
    call reallocate(ncbd,0)
    call reallocate(ii1,0)
    call reallocate(jj1,0)
    call reallocate(kk1,0)
    call reallocate(ii2,0)
    call reallocate(jj2,0)
    call reallocate(kk2,0)
    call reallocate(id1,0)
    call reallocate(jd1,0)
    call reallocate(kd1,0)
    call reallocate(id2,0)
    call reallocate(jd2,0)
    call reallocate(kd2,0)
    call reallocate(nnn,0)
    call reallocate(nnc,0)
    call reallocate(nnfb,0)
    call reallocate(npn,0)
    call reallocate(npc,0)
    call reallocate(npfb,0)
    call reallocate(cl,0)
    call reallocate(ndcc,0)
    call reallocate(nbdc,0)
    call reallocate(nfbc,0)
    call reallocate(bc,0,0)
    call reallocate(mdnc,0)
    call reallocate(mper,0)
    call reallocate(mpc,0)
    call reallocate(ncin,0)
    call reallocate(bceqt,0,0)
    call reallocate(mnc,0)

    call reallocate(bcg_to_proc,0)
    call reallocate(bcg_to_bcl,0)
    call reallocate(bcg_to_bci,0)
    call reallocate(bcl_to_bcg,0)
    call reallocate(bg_to_proc,0)
    call reallocate(bg_to_bl,0)
    call reallocate(bg_to_bi,0)
    call reallocate(bl_to_bg,0)
    call reallocate(bcg_to_bg,0)

    ndimubx = 0
    ndimctbx=0
    ndimntbx=0
    lzx =0
    klzx=0
    lt=0
    mtbx=0
    mtcx=0
    kmtbx=0
    mtb=0
    mtt=0
    mdimubx=0
    mdimtbx=0
    mfbe=0
    ip41=0
    vbc=0.
    num_bcg=0
    num_bci=0
    num_bcl=0
    num_bg=0
    num_bi=0
    num_bl=0

    !############################################################################################
    !####################### COMPUTE SPLITTING ##################################################
    !############################################################################################


    ! nouveaux tableaux
    allocate(nblock2(save_num_bg),nblockd(3,save_lt),ni(1,1,1,nblocks),nj(1,1,1,nblocks),nk(1,1,1,nblocks))
    allocate(ii2g(save_num_bg),jj2g(save_num_bg),kk2g(save_num_bg),nblockdg(3,save_num_bg))
    allocate(sub_bc(0,0),ni1(1,1,1),nj1(1,1,1),nk1(1,1,1))

    nblock2=1    ! initial number of splitting for each existing block
    nblockd=1    ! initial number of splitting for each existing block
    nblockdg=0

    do l=1,save_num_bg
      l1=save_bg_to_bl(l)
      if (l1/=0) then 
        val(1)=save_ii2(l1)
        val(2)=save_jj2(l1)
        val(3)=save_kk2(l1)
      endif
      call bcast(val,save_bg_to_proc(l))
      ii2g(l)=val(1)
      jj2g(l)=val(2)
      kk2g(l)=val(3)
    enddo
    nxyza=sum(ii2g*jj2g*kk2g)  ! total number of points

    ! routine calculant combiens de fois splitter chaque block
    ! sortie : nblock2
    ! entrée : tout le reste
    call num_split(nblock2,save_num_bg,nxyza,nblocks,ii2g,jj2g,kk2g)

    do l=1,save_lt
      l1=save_bl_to_bg(l)

       ! calcule le split pour un block, c'est à dire 
       ! le nombre de découpe par direction            (nblockd)
       ! le nombre de points dans chaque nouveau block (nicv2,njcv2,nkcv2)
       ! une estimation des communications             (num_cf2)
       call triv_split(nblock2(l1),l,nxyza,save_ii2 ,save_jj2 ,save_kk2, &
            nblockd(:,l),num_cf2,ni1,nj1,nk1)

       ! ajoute le split avec les splits des autres blocks
       call reallocate_s(ni,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nblocks)
       ni(:size(ni1,1),:size(ni1,2),:size(ni1,3),l1)=ni1

       call reallocate_s(nj,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nblocks)
       nj(:size(nj1,1),:size(nj1,2),:size(nj1,3),l1)=nj1

       call reallocate_s(nk,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nblocks)
       nk(:size(nk1,1),:size(nk1,2),:size(nk1,3),l1)=nk1
    enddo

    do l=1,save_num_bg
       l1=save_bg_to_bl(l)
       proc=save_bg_to_proc(l)
       if(l1/=0) nblockdg(:,l)=nblockd(:,l1)
       call bcast(nblockdg(:,l),proc)
       call reallocate_s(ni,maxval(nblockdg(1,:)),maxval(nblockdg(2,:)),maxval(nblockdg(3,:)),nblocks)
       call reallocate_s(nj,maxval(nblockdg(1,:)),maxval(nblockdg(2,:)),maxval(nblockdg(3,:)),nblocks)
       call reallocate_s(nk,maxval(nblockdg(1,:)),maxval(nblockdg(2,:)),maxval(nblockdg(3,:)),nblocks)
       call bcast(ni(:,:,:,l),proc)
       call bcast(nj(:,:,:,l),proc)
       call bcast(nk(:,:,:,l),proc)
    enddo

    sblock=sum(nblockdg(1,:)*nblockdg(2,:)*nblockdg(3,:))
    if(sblock/=nblocks) then
       stop 'partitionnement impossible'
    else
       if(verbosity>=1.and.rank==0) then
         print*,'découpage réussis : '
         nmax=0
         nmin=nxyza
         do l=1,save_num_bg
            nmax1=maxval(ni(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l)* &
                          nj(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l)* &
                          nk(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l))
            nmin1=minval(ni(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l)* &
                          nj(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l)* &
                          nk(:nblockdg(1,l),:nblockdg(2,l),:nblockdg(3,l),l))
            print*, l,nblockdg(:,l),( nmax1- nmin1 )*100./ nmin1,"%"
            nmax=max(nmax,nmax1)
            nmin=min(nmin,nmin1)
         enddo
         print*, "Total : ",sblock,( nmax- nmin )*100./ nmin,"%"
       endif
       if(verbosity>=2.and.rank==0) then
          do l=1,save_num_bg
            do k=1,nblockdg(3,l)
              do j=1,nblockdg(2,l)
                  print*, l,k,j,ni(:nblockdg(1,l),j,k,l)*nj(:nblockdg(1,l),j,k,l)*nk(:nblockdg(1,l),j,k,l)
              enddo
            enddo
          enddo
       endif
    end if


    !############################################################################################
    !######################### RECREATE GRID ####################################################
    !############################################################################################
    allocate(new2old_b(nblocks))

    if(verbosity>=2) then
      mot="" ; nmot=7 ; imot=0
      mot(1)="create" ; imot(1)=6
      mot(2)="dom"    ; imot(2)=3
      mot(3)="st"     ; imot(3)=2
    endif
!    print*,'recreate grid '
    do l=1,save_num_bg  
      l1=save_bg_to_bl(l)
       do k=1,nblockdg(3,l)
          do j=1,nblockdg(2,l)
             do i=1,nblockdg(1,l)!     create new grid and initialize it

                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockdg(1,l)+(k-1)*nblockdg(2,l)*nblockdg(1,l)

                if(verbosity>=2) then
                  call str(mot,imot,nmx,4 ,l2)
                  call str(mot,imot,nmx,5 ,ni(i,j,k,l))
                  call str(mot,imot,nmx,6 ,nj(i,j,k,l))
                  call str(mot,imot,nmx,7 ,nk(i,j,k,l))
                  call c_crdms( mot,imot,nmot,save_bg_to_bi(l))
                else
                  call crdms( l2,ni(i,j,k,l),nj(i,j,k,l),nk(i,j,k,l),save_bg_to_bi(l))
                endif

                orig=save_bg_to_proc(l)
                dest=bg_to_proc(l2)
                l3=bg_to_bl(l2)


                if(rank==orig.or.rank==dest) then
                  if(orig==dest) then ! sending to myself

                    ! offset, don't forget to count the interface twice
                    call new2old_p(0,0,0,xs,ys,zs,i,j,k,l)

                    call reallocate_s(x,ndimntbx)
                    call reallocate_s(y,ndimntbx)
                    call reallocate_s(z,ndimntbx)

                    call copy_grid(l1,l3,xs,ys,zs,x,y,z,save_x,save_y,save_z,              &
                         save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn)

                  elseif(rank==orig) then

                    ! offset, don't forget to count the interface twice
                    call new2old_p(1,1,1,xs,ys,zs,i,j,k,l)
                    call new2old_p(ni(i,j,k,l),nj(i,j,k,l),nk(i,j,k,l),xe,ye,ze,i,j,k,l)

                    nid = save_id2(l1)-save_id1(l1)+1
                    njd = save_jd2(l1)-save_jd1(l1)+1
                    nijd = nid*njd

                    call send_grid(save_x,save_y,save_z,xs,ys,zs,xe,ye,ze,       &
                        save_id1(l1),save_jd1(l1),save_kd1(l1),nid,nijd,save_npn(l1),orig,dest)

                  elseif(rank==dest) then

                    call reallocate_s(x,ndimntbx)
                    call reallocate_s(y,ndimntbx)
                    call reallocate_s(z,ndimntbx)

                    nid = id2(l3)-id1(l3)+1
                    njd = jd2(l3)-jd1(l3)+1
                    nijd = nid*njd

                    call recv_grid(x,y,z,ii1(l3),jj1(l3),kk1(l3),ii2(l3),jj2(l3),kk2(l3),       &
                            id1(l3),jd1(l3),kd1(l3),nid,nijd,npn(l3),orig,dest)
                  endif
                endif

                new2old_b(l2)= l

             enddo
          enddo
       enddo
    enddo

    !############################################################################################
    !######################### RECREATE OLD BOUNDARIES ##########################################
    !############################################################################################
    if(verbosity>=2) then
      mot="" ; nmot=13  ; imot=0
      mot(1)="create"   ; imot(1)=6
      mot(2)="boundary" ; imot(2)=8
      mot(3)="st"       ; imot(3)=2
    endif

!    print*,'recreate old boundaries '
    do fr1=1,save_num_bcg
       fr=save_bcg_to_bcl(fr1)
       l=save_bcg_to_bg(fr1)
       l1=save_bg_to_bl(l)
       orig1=save_bg_to_proc(l)
       do k=1,nblockdg(3,l)
          do j=1,nblockdg(2,l)
             do i=1,nblockdg(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockdg(1,l)+(k-1)*nblockdg(2,l)*nblockdg(1,l)
                dest=bg_to_proc(l2)

                nsub=0 ! count sub_boundaries
                call reallocate(sub_bc,1,6)

                if (rank==orig1) then

                    call new2old_p(1,1,1,xs,ys,zs,i,j,k,l)
                    call new2old_p(ni(i,j,k,l),nj(i,j,k,l),nk(i,j,k,l),xe,ye,ze,i,j,k,l)

                    if ( save_iminb(fr)<=xe .and. &
                         save_imaxb(fr)>=xs .and. &
                         save_jminb(fr)<=ye .and. &
                         save_jmaxb(fr)>=ys .and. & ! there is a part of the boundary in this block
                         save_kminb(fr)<=ze .and. &
                         save_kmaxb(fr)>=zs ) then

                         nsub=1

                       ! part of the boundary which concern this block

                        sub_bc(1,1)=min(xe,max(xs,save_iminb(fr))) ! coordinate of the new boundary
                        sub_bc(1,2)=min(xe,max(xs,save_imaxb(fr))) ! in the old block ref
                        sub_bc(1,3)=min(ye,max(ys,save_jminb(fr)))
                        sub_bc(1,4)=min(ye,max(ys,save_jmaxb(fr)))
                        sub_bc(1,5)=min(ze,max(zs,save_kminb(fr)))
                        sub_bc(1,6)=min(ze,max(zs,save_kmaxb(fr)))

                        indmf=save_indfl(fr)
                   endif
                endif
                call bcast(nsub,orig1)

                if (nsub==0) cycle

                 if(tab_raccord(fr1)/=0) then ! if boundary is shared with another block
                    fr2=tab_raccord(fr1)      ! a split on the other side induce split here
                    fr3=save_bcg_to_bcl(fr2)
                    l3=save_bcg_to_bg(fr2)
                    l4=save_bg_to_bl(l3)
                    orig2=save_bg_to_proc(l3)

                    if(rank==orig1) then
                      sub_bc(1,1)=sub_bc(1,1)-save_iminb(fr)
                      sub_bc(1,2)=sub_bc(1,2)-save_iminb(fr) ! coordinate of the new boundary
                      sub_bc(1,3)=sub_bc(1,3)-save_jminb(fr) ! in the old boundary ref
                      sub_bc(1,4)=sub_bc(1,4)-save_jminb(fr)
                      sub_bc(1,5)=sub_bc(1,5)-save_kminb(fr)
                      sub_bc(1,6)=sub_bc(1,6)-save_kminb(fr)
                    endif

                    call MPI_TRANS(sub_bc,sub_bc,orig1,orig2)

                    if(rank==orig2) then
                      nsub=0

                      imin=sub_bc(1,1)+save_iminb(fr3)
                      imax=sub_bc(1,2)+save_iminb(fr3) ! coordinate of the new boundary
                      jmin=sub_bc(1,3)+save_jminb(fr3) ! in the old block ref
                      jmax=sub_bc(1,4)+save_jminb(fr3)
                      kmin=sub_bc(1,5)+save_kminb(fr3)
                      kmax=sub_bc(1,6)+save_kminb(fr3)

                      do k2=1,nblockdg(3,l3)
                         do j2=1,nblockdg(2,l3)
                            do i2=1,nblockdg(1,l3)

                                call new2old_p(imin,jmin,kmin,xs,ys,zs,i2,j2,k2,l3)
                                call new2old_p(imax,jmax,kmax,xe,ye,ze,i2,j2,k2,l3)

                                if ( save_iminb(fr3)<=xe .and. &
                                     save_imaxb(fr3)>=xs .and. &
                                     save_jminb(fr3)<=ye .and. &
                                     save_jmaxb(fr3)>=ys .and. & ! there is a part of the boundary in this block
                                     save_kminb(fr3)<=ze .and. &
                                     save_kmaxb(fr3)>=zs ) then

                                   ! part of the boundary which concern this block

                                    nsub=nsub+1
                                    call reallocate_s(sub_bc,nsub,6)
                                    sub_bc(nsub,1)=min(xe,max(xs,save_iminb(fr3)))-save_iminb(fr3) ! coordinate of the new boundary
                                    sub_bc(nsub,2)=min(xe,max(xs,save_imaxb(fr3)))-save_iminb(fr3) ! in the old boundary ref
                                    sub_bc(nsub,3)=min(ye,max(ys,save_jminb(fr3)))-save_jminb(fr3)
                                    sub_bc(nsub,4)=min(ye,max(ys,save_jmaxb(fr3)))-save_jminb(fr3)
                                    sub_bc(nsub,5)=min(ze,max(zs,save_kminb(fr3)))-save_kminb(fr3)
                                    sub_bc(nsub,6)=min(ze,max(zs,save_kmaxb(fr3)))-save_kminb(fr3)
                               endif
                            enddo
                         enddo
                      enddo
                    endif
                    call MPI_TRANS(nsub,nsub,orig2,orig1)
                    if(rank==orig1) call reallocate_s(sub_bc,nsub,6)
                    call MPI_TRANS(sub_bc,sub_bc,orig2,orig1)
                    if(rank==orig1) then
                      sub_bc(:,1)=sub_bc(:,1)+save_iminb(fr)
                      sub_bc(:,2)=sub_bc(:,2)+save_iminb(fr) ! coordinate of the new boundary
                      sub_bc(:,3)=sub_bc(:,3)+save_jminb(fr) ! in the old block ref
                      sub_bc(:,4)=sub_bc(:,4)+save_jminb(fr)
                      sub_bc(:,5)=sub_bc(:,5)+save_kminb(fr)
                      sub_bc(:,6)=sub_bc(:,6)+save_kminb(fr)
                    endif
                  endif
                  call bcast(nsub,orig1)
                  if(rank/=orig1) call reallocate(sub_bc,nsub,6)
                  call bcast(sub_bc,orig1)
                  call bcast(indmf,orig1)
                  do i2=1,nsub
                      mfbe=mfbe+1
                      call old2new_p(sub_bc(i2,1),sub_bc(i2,3),sub_bc(i2,5),sub_bc(i2,1),sub_bc(i2,3),sub_bc(i2,5),i,j,k,l)
                      call old2new_p(sub_bc(i2,2),sub_bc(i2,4),sub_bc(i2,6),sub_bc(i2,2),sub_bc(i2,4),sub_bc(i2,6),i,j,k,l)
                      if(verbosity>=2) then
                        call str(mot,imot,nmx,4 ,mfbe)
                        call str(mot,imot,nmx,5 ,1)
                        call str(mot,imot,nmx,6 ,l2)
                        call str(mot,imot,nmx,7 ,sub_bc(i2,1))
                        call str(mot,imot,nmx,8 ,sub_bc(i2,2))
                        call str(mot,imot,nmx,9 ,sub_bc(i2,3))
                        call str(mot,imot,nmx,10,sub_bc(i2,4))
                        call str(mot,imot,nmx,11,sub_bc(i2,5))
                        call str(mot,imot,nmx,12,sub_bc(i2,6))
                        mot(13)=indmf ; imot(13)=2
                        call c_crbds( mot,imot,nmot, ncbd,save_bcg_to_bci(fr1))
                      else
                        call crbds( &
                             mfbe,1,l2, &
                             sub_bc(i2,1),sub_bc(i2,2),sub_bc(i2,3),sub_bc(i2,4),sub_bc(i2,5),sub_bc(i2,6), &
                             indmf, &
                             ncbd,save_bcg_to_bci(fr1))
                      endif
                  enddo
            enddo
         enddo
      enddo
    enddo

    !############################################################################################
    !########################### CREATE NEW BOUNDARIES ##########################################
    !############################################################################################

!    print*,'create new boundaries '
    call reallocate(sub_bc,2,6)
    do l=1,save_num_bg  
      l1=save_bg_to_bl(l)
       do k=1,nblockdg(3,l)
          do j=1,nblockdg(2,l)
             do i=1,nblockdg(1,l) !           New coincident boundaries
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockdg(1,l)+(k-1)*nblockdg(2,l)*nblockdg(1,l)
                ll2=bg_to_bl(l2)
                orig1=bg_to_proc(l2)
                sub_bc=0

                if (i>1) then
                  l3=sum(nblock2(1:l-1))+i-1+(j-1)*nblockdg(1,l)+(k-1)*nblockdg(2,l)*nblockdg(1,l)
                  ll3=bg_to_bl(l3)
                  orig2=bg_to_proc(l3)

                  if (rank==orig1) then
                      sub_bc(1,1)=ii1(ll2)
                      sub_bc(1,2)=ii1(ll2)
                      sub_bc(1,3)=jj1(ll2)
                      sub_bc(1,4)=jj2(ll2)
                      sub_bc(1,5)=kk1(ll2)
                      sub_bc(1,6)=kk2(ll2)
                  elseif (rank==orig2) then
                      sub_bc(2,1)=ii2(ll3)
                      sub_bc(2,2)=ii2(ll3)
                      sub_bc(2,3)=jj1(ll3)
                      sub_bc(2,4)=jj2(ll3)
                      sub_bc(2,5)=kk1(ll3)
                      sub_bc(2,6)=kk2(ll3)
                  endif

                  call bcast(sub_bc(1,:),orig1)
                  call bcast(sub_bc(2,:),orig2)

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,sub_bc(1,1))
                     call str(mot,imot,nmx,8 ,sub_bc(1,2))
                     call str(mot,imot,nmx,9 ,sub_bc(1,3))
                     call str(mot,imot,nmx,10,sub_bc(1,4))
                     call str(mot,imot,nmx,11,sub_bc(1,5))
                     call str(mot,imot,nmx,12,sub_bc(1,6))
                     mot(13)="i1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                           sub_bc(1,1),sub_bc(1,2),sub_bc(1,3),sub_bc(1,4),sub_bc(1,5),sub_bc(1,6), &
                          'i1', &
                          ncbd,0)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,sub_bc(2,1))
                     call str(mot,imot,nmx,8 ,sub_bc(2,2))
                     call str(mot,imot,nmx,9 ,sub_bc(2,3))
                     call str(mot,imot,nmx,10,sub_bc(2,4))
                     call str(mot,imot,nmx,11,sub_bc(2,5))
                     call str(mot,imot,nmx,12,sub_bc(2,6))
                     mot(13)="i2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          sub_bc(2,1),sub_bc(2,2),sub_bc(2,3),sub_bc(2,4),sub_bc(2,5),sub_bc(2,6), &
                          'i2', &
                          ncbd,0)
                    endif
                endif
                if (j>1) then
                  l3=sum(nblock2(1:l-1))+i+(j-2)*nblockdg(1,l)+(k-1)*nblockdg(2,l)*nblockdg(1,l)
                  ll3=bg_to_bl(l3)
                  orig2=bg_to_proc(l3)

                  if (rank==orig1) then
                      sub_bc(1,1)=ii1(ll2)
                      sub_bc(1,2)=ii2(ll2)
                      sub_bc(1,3)=jj1(ll2)
                      sub_bc(1,4)=jj1(ll2)
                      sub_bc(1,5)=kk1(ll2)
                      sub_bc(1,6)=kk2(ll2)
                  elseif (rank==orig2) then
                      sub_bc(2,1)=ii1(ll3)
                      sub_bc(2,2)=ii2(ll3)
                      sub_bc(2,3)=jj2(ll3)
                      sub_bc(2,4)=jj2(ll3)
                      sub_bc(2,5)=kk1(ll3)
                      sub_bc(2,6)=kk2(ll3)
                  endif

                  call bcast(sub_bc(1,:),orig1)
                  call bcast(sub_bc(2,:),orig2)

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,sub_bc(1,1))
                     call str(mot,imot,nmx,8 ,sub_bc(1,2))
                     call str(mot,imot,nmx,9 ,sub_bc(1,3))
                     call str(mot,imot,nmx,10,sub_bc(1,4))
                     call str(mot,imot,nmx,11,sub_bc(1,5))
                     call str(mot,imot,nmx,12,sub_bc(1,6))
                     mot(13)="j1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                           sub_bc(1,1),sub_bc(1,2),sub_bc(1,3),sub_bc(1,4),sub_bc(1,5),sub_bc(1,6), &
                          'j1', &
                          ncbd,0)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,sub_bc(2,1))
                     call str(mot,imot,nmx,8 ,sub_bc(2,2))
                     call str(mot,imot,nmx,9 ,sub_bc(2,3))
                     call str(mot,imot,nmx,10,sub_bc(2,4))
                     call str(mot,imot,nmx,11,sub_bc(2,5))
                     call str(mot,imot,nmx,12,sub_bc(2,6))
                     mot(13)="j2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          sub_bc(2,1),sub_bc(2,2),sub_bc(2,3),sub_bc(2,4),sub_bc(2,5),sub_bc(2,6), &
                          'j2', &
                          ncbd,0)
                    endif
                endif
                if (k>1) then
                  l3=sum(nblock2(1:l-1))+i+(j-1)*nblockdg(1,l)+(k-2)*nblockdg(2,l)*nblockdg(1,l)
                  ll3=bg_to_bl(l3)
                  orig2=bg_to_proc(l3)

                  if (rank==orig1) then
                      sub_bc(1,1)=ii1(ll2)
                      sub_bc(1,2)=ii2(ll2)
                      sub_bc(1,3)=jj1(ll2)
                      sub_bc(1,4)=jj2(ll2)
                      sub_bc(1,5)=kk1(ll2)
                      sub_bc(1,6)=kk1(ll2)
                  elseif (rank==orig2) then
                      sub_bc(2,1)=ii1(ll3)
                      sub_bc(2,2)=ii2(ll3)
                      sub_bc(2,3)=jj1(ll3)
                      sub_bc(2,4)=jj2(ll3)
                      sub_bc(2,5)=kk2(ll3)
                      sub_bc(2,6)=kk2(ll3)
                  endif

                  call bcast(sub_bc(1,:),orig1)
                  call bcast(sub_bc(2,:),orig2)

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,sub_bc(1,1))
                     call str(mot,imot,nmx,8 ,sub_bc(1,2))
                     call str(mot,imot,nmx,9 ,sub_bc(1,3))
                     call str(mot,imot,nmx,10,sub_bc(1,4))
                     call str(mot,imot,nmx,11,sub_bc(1,5))
                     call str(mot,imot,nmx,12,sub_bc(1,6))
                     mot(13)="k1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                           sub_bc(1,1),sub_bc(1,2),sub_bc(1,3),sub_bc(1,4),sub_bc(1,5),sub_bc(1,6), &
                          'k1', &
                          ncbd,0)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,sub_bc(2,1))
                     call str(mot,imot,nmx,8 ,sub_bc(2,2))
                     call str(mot,imot,nmx,9 ,sub_bc(2,3))
                     call str(mot,imot,nmx,10,sub_bc(2,4))
                     call str(mot,imot,nmx,11,sub_bc(2,5))
                     call str(mot,imot,nmx,12,sub_bc(2,6))
                     mot(13)="k2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd,0)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          sub_bc(2,1),sub_bc(2,2),sub_bc(2,3),sub_bc(2,4),sub_bc(2,5),sub_bc(2,6), &
                          'k2', &
                          ncbd,0)
                    endif
                endif
             enddo
          enddo
       enddo
    enddo



    !############################################################################################
    !############################# SAVE NEW MESH FOR CHECKING PURPOSE ###########################
    !############################################################################################

   if(verbosity>=3) then
    ! write new grid
      do l=1,lt
         write(fich,'(A,I0.2,A)') "testmesh_",bl_to_bg(l),".dat"
         open(42,file=fich,status="replace")

         do k=kk1(l),kk2(l)
            do j=jj1(l),jj2(l)
               do i=ii1(l),ii2(l)

                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),bl_to_bg(l)

               enddo
               write(42,*) ""
            enddo
         enddo
         close(42)
      enddo

    ! write new boundaries
      do fr=1,mtb

      write(fich,'(A,I0.2,A)') "testbnd_",bcl_to_bcg(fr),".dat"
      open(42,file=fich,status="replace")

         do k=kminb(fr),kmaxb(fr)
            do j=jminb(fr),jmaxb(fr)
               do i=iminb(fr),imaxb(fr)
                  l=ndlb(fr)
                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),bcl_to_bcg(fr)

               enddo
               if(indfl(fr)(1:1)/="i") write(42,*) ""
            enddo
            if(indfl(fr)(1:1)=="i") write(42,*) ""
         enddo
      enddo
      close(42)
    endif

 deallocate(save_x,save_y,save_z,save_ndlb,save_nfei,save_indfl,save_mpb,save_mmb,save_ncbd,save_ii1) ! TESTED AND OK
 deallocate(save_jj1,save_kk1,save_ii2,save_jj2,save_kk2,save_id1,save_jd1,save_kd1,save_id2,save_jd2)
 deallocate(save_kd2,save_nnn,save_nnc,save_nnfb,save_npn,save_npc,save_npfb)

    !############################################################################################
    !################### INITIALIZE COINCIDENT BOUNDARIES #######################################
    !############################################################################################
    ip21=ndimntbx
    ip40=mdimubx            ! Nb point frontiere
    ip41=mdimtbx            ! Nb point frontiere
    ip42=mdimtbx            ! Nb point frontiere
    ip43=mdimtbx            ! Nb point frontiere
    ip44=0!mdimubx            !TODO
    test=.false.
    call reallocate(cl,mtb)
    call reallocate(ndcc,mtb)
    call reallocate(nbdc,mtb)
    call reallocate(nfbc,mtb)
    call reallocate(bc,mtb,ista*lsta)
    call reallocate(mdnc,mtt)
    call reallocate(mper,mtt)
    call reallocate(mpc,mtt)
    call reallocate(ncin,ip41)
    call reallocate(bceqt,ip41,neqt)
    call reallocate(mnc,ip43)
    call reallocate(sub_bc,2,6)

!    print*,'initialization '
    do fr1=1,num_bcg
       fr=bcg_to_bcl(fr1)
       l=bcg_to_bg(fr1)
       l1=bg_to_bl(l)
       fri=bcg_to_bci(fr1)
       test=.false.
       orig1=bcg_to_proc(fr1)
       if (fri==0) then  ! new boundary , it's easy
         test=.true.
         if (rank==orig1) then
           if (indfl(fr)(2:2)=="1") fr2=fr1+1
           if (indfl(fr)(2:2)=="2") fr2=fr1-1
         endif
      elseif(tab_raccord(fri)/=0) then ! old raccord boundary
         test=.true.
         if (rank==orig1) &
            call get_coords_box(sub_bc1(1,1),sub_bc1(1,2),sub_bc1(1,3),sub_bc1(1,4),sub_bc1(1,5),sub_bc1(1,6),                 &
                  iminb(fr),imaxb(fr),jminb(fr),jmaxb(fr),kminb(fr),kmaxb(fr), &
                  id1(l1),id2(l1),jd1(l1),jd2(l1),kd1(l1),kd2(l1),npn(l1),       &
                  x,y,z)
         call bcast(sub_bc1,orig1)

         fr4=0
         find_otherblock: do fr2=1,mtb
            fr3=bcl_to_bcg(fr2)
            l2=bcg_to_bg(fr3)
            l3=bg_to_bl(l2)
            fri2=bcg_to_bci(fr3)
            if(fri2/=0)then
            if(tab_raccord(fri2)==fri) then  ! potential new boundary number

              call get_coords_box(sub_bc1(2,1),sub_bc1(2,2),sub_bc1(2,3),sub_bc1(2,4),sub_bc1(2,5),sub_bc1(2,6),                 &
                    iminb(fr2),imaxb(fr2),jminb(fr2),jmaxb(fr2),kminb(fr2),kmaxb(fr2), &
                    id1(l3),id2(l3),jd1(l3),jd2(l3),kd1(l3),kd2(l3),npn(l3),       &
                    x,y,z)

              if (abs(sub_bc1(1,1)-sub_bc1(2,1))<=1d-10 .and. &
                  abs(sub_bc1(1,2)-sub_bc1(2,2))<=1d-10 .and. &
                  abs(sub_bc1(1,3)-sub_bc1(2,3))<=1d-10 .and. & ! It's me ! 
                  abs(sub_bc1(1,4)-sub_bc1(2,4))<=1d-10 .and. &
                  abs(sub_bc1(1,5)-sub_bc1(2,5))<=1d-10 .and. &
                  abs(sub_bc1(1,6)-sub_bc1(2,6))<=1d-10) then

                  fr4=fr3
                  exit find_otherblock
              endif
            endif
            endif
         enddo find_otherblock
         call sum_mpi(fr4) ! there should be only one non zero value in this sum
         fr2=fr4
      endif
      if (test) then   ! raccord boundary

         l2=bcg_to_bg(fr2)
         if (rank==bg_to_proc(l2)) then
           sub_bc(1,1)=iminb(bcg_to_bcl(fr2))
           sub_bc(1,2)=jminb(bcg_to_bcl(fr2))
           sub_bc(1,3)=kminb(bcg_to_bcl(fr2))
         endif
         call bcast(sub_bc(1,1:3),bcg_to_proc(fr2))


           if(verbosity>=2) then
             mot="" ; nmot=6   ; imot=0
             mot(1)="init"     ; imot(1)=4
             mot(2)="boundary" ; imot(2)=8
             mot(3)="basic"    ; imot(3)=5
             call str(mot,imot,nmx,4 ,fr1)
             mot(5)="rc"       ; imot(5)=2
             call str(mot,imot,nmx,6 ,1)
             call c_inbdb( mot,imot,nmot,ncbd,ncin,bceqt,partition=.true.)
           else
            call inbdb( &
                 ncbd,ncin, &
                 fr1,"rc  ",1, &
                 0,0,0,0,vbc,bceqt)
            endif

         if (rank==orig1) indmf=indfl(fr)
         call bcast(indmf,orig1)

          select case(indmf(1:1))
          case("i")
               mot(16)="fa"      ; imot(16)=2
               mot(17)="+j"      ; imot(17)=2
               mot(18)="+k"      ; imot(18)=2
          case("j")
             mot(16)="+i"      ; imot(16)=2
             mot(17)="fa"      ; imot(17)=2
             mot(18)="+k"      ; imot(18)=2
          case("k")
             mot(16)="+i"      ; imot(16)=2
             mot(17)="+j"      ; imot(17)=2
             mot(18)="fa"      ; imot(18)=2
          end select

           if(verbosity>=2) then
             nmot=18  
             mot(1)="init"     ; imot(1)=4
             mot(2)="boundary" ; imot(2)=8
             mot(3)="coin"     ; imot(3)=4
             call str(mot,imot,nmx,4 ,fr1)
             mot(5)="frc"      ; imot(5)=3
             call str(mot,imot,nmx,6 ,fr2)
             mot(7)="kibdc"    ; imot(7)=5
             call str(mot,imot,nmx,8 ,1)
             mot(9)="krr"      ; imot(9)=3
             call str(mot,imot,nmx,10 ,0)
             mot(11)="ptc"     ; imot(11)=3
             call str(mot,imot,nmx,12 ,sub_bc(1,1))
             call str(mot,imot,nmx,13 ,sub_bc(1,2))
             call str(mot,imot,nmx,14 ,sub_bc(1,3))
             mot(15)="dir"     ; imot(15)=3
             call c_inbdc(  mot,imot,nmot, exs1,exs2, x,y,z, ncbd,ncin,mnc)
           else
            call inbdc( &
                 exs1,exs2, &
                 x,y,z, &
                 ncbd,ncin,mnc, &
                 0,fr1,fr2,1,0., &
                 sub_bc(1,1),sub_bc(1,2),sub_bc(1,3),mot(16)(1:2),mot(17)(1:2),mot(18)(1:2))
           endif
      endif
!    print*,l,' : filling done'
   enddo
endif
    write(stderr,*) "done"
 return


contains

subroutine old2new_p(old_i,old_j,old_k,new_i,new_j,new_k,i,j,k,l)
implicit none
integer,intent(in) :: old_i,old_j,old_k,i,j,k,l
integer,intent(out)  :: new_i,new_j,new_k
new_i=old_i-sum(ni(:i-1,j,k,l))+i-1
new_j=old_j-sum(nj(i,:j-1,k,l))+j-1
new_k=old_k-sum(nk(i,j,:k-1,l))+k-1
end subroutine old2new_p

subroutine new2old_p(new_i,new_j,new_k,old_i,old_j,old_k,i,j,k,l)
implicit none
integer,intent(out) :: old_i,old_j,old_k
integer,intent(in)  :: new_i,new_j,new_k,i,j,k,l
old_i=new_i+sum(ni(:i-1,j,k,l))-i+1
old_j=new_j+sum(nj(i,:j-1,k,l))-j+1
old_k=new_k+sum(nk(i,j,:k-1,l))-k+1
end subroutine new2old_p


end subroutine partitionnement

subroutine num_split(nblock2,lt,nxyza,nblocks,ii2,jj2,kk2)
 implicit none
 integer,intent(in)  :: lt,nxyza,nblocks
 integer,intent(in)  :: ii2(lt),jj2(lt),kk2(lt)
 integer,intent(out) :: nblock2(lt)
 integer             :: rsize,sblock(lt),i,j,k,nblock(lt)
 double precision    :: unbalance,unbalance1

 ! TODO
 ! switch to the alternative version which permit to have ideal blocks size
 ! need a criteria to avoid too small block, and need to manage the residual block
 ! todo

 !   compute number of spliting of each blocks with the best equilibrium
! do i=lt,nblocks-1                       ! split until lt>=nblocks
!    unbalance=10000 ! a lot
!    do j=1,lt
!      nblock=nblock2                    ! try every split
!      nblock(j)=nblock(j)+1
!      sblock=ceiling(ii2*jj2*kk2*1./nblock) ! compute the current size of blocks (approx.)
!      unbalance1=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
!      if (unbalance1< unbalance) then ! if better remember it
!        unbalance = unbalance1
!        k=j
!      endif
!    end do
!    nblock2(k)=nblock2(k)+1
! end do


 !   compute number of spliting of each blocks with the best equilibrium
 do i=lt,nblocks-1                       ! split until lt>=nblocks
    sblock=ceiling(ii2*jj2*kk2*1./nblock2) ! compute the current size of blocks
    j=maxloc(sblock,1)                        ! split the first bigest block
    do k=j+1,lt
       if(sblock(k)==sblock(j) &          ! if more than one bigest block
            .and.nblock2(k)>nblock2(j)) &      ! split the most splitted
            j=k
    end do
    nblock2(j)=nblock2(j)+1
 end do


!    compute number of spliting of each blocks with the ideal equilibrium
!     rsize=nint(nxyza*1./nblocks)              ! ideal size of a block
!     nblock2=ceiling(ii2*jj2*kk2*1./rsize) ! number of split needed
!!     sblock=mod(ii2*jj2*kk2,rsize)         ! size of the smallest block
!!     do i=1,lt
!!       if (sblock(i) <= 3*3*3) &          ! allow for small imbalance in order to avoid too small blocks
!!            nblock2(i)=nblock2(i)-1
!!     end do
! do i=sum(nblock2),nblocks-1,-1
!    sblock=ceiling(ii2*jj2*kk2*1./nblock2) ! compute the current size of blocks (approx.)
!    unbalance=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
!    k=0
!    do j=1,lt
!      nblock=nblock2                    ! try every split
!      nblock(j)=nblock(j)-1
!      sblock=ceiling(ii2*jj2*kk2*1./nblock) ! compute the current size of blocks (approx.)
!      unbalance1=(maxval(sblock)-minval(sblock))*1./maxval(sblock) ! and unbalance
!      if (unbalance1< unbalance) then ! if better remember it
!        unbalance = unbalance1
!        k=j
!      endif
!    end do
!    print*,sum(nblock2),unbalance
!    if (k==0) exit
!    nblock2(k)=nblock2(k)-1
! end do


end subroutine num_split

subroutine triv_split(nblock2,nbl,nxyza,ii2,jj2,kk2, &
    nblockd,num_cf2,new_ii2,new_jj2,new_kk2)
 implicit none
 integer,allocatable,intent(in)  :: ii2(:),jj2(:),kk2(:)
 integer,intent(in)              :: nxyza,nbl,nblock2

 integer,allocatable,intent(inout) :: new_ii2(:,:,:),new_jj2(:,:,:),new_kk2(:,:,:), num_cf2(:,:,:)
 integer,intent(out)             :: nblockd(3)

 integer             :: i,j,k,i1,j1,k1
 integer,allocatable :: tmp_ii2(:,:,:),tmp_jj2(:,:,:),tmp_kk2(:,:,:),num_cft(:,:,:)

 !   trivial spliting : divide my block in nblock2 subblock
 !                      test all possiblities constisting in dividing
 !                      i times in the x direction, j times in the y direction and k times in the z direction
 call reallocate(num_cf2,1,1,1)
 num_cf2=nxyza*nblock2 ! useless initial big value
 do k=1,nblock2
    do j=1,nblock2
       do i=1,nblock2
          if(i*j*k==nblock2) then !           if we get the right number of blocks
             allocate(tmp_jj2(i,j,k),tmp_ii2(i,j,k),tmp_kk2(i,j,k), num_cft(i,j,k))
             !       compute sizes of sub-blocks
             do k1=1,k
                do j1=1,j
                   do i1=1,i
                      tmp_ii2(i1,j1,k1)=nint(i1*(ii2(nbl)+i-1)*1./i) - nint((i1-1.)*(ii2(nbl)+i-1)*1./i)
                      tmp_jj2(i1,j1,k1)=nint(j1*(jj2(nbl)+j-1)*1./j) - nint((j1-1.)*(jj2(nbl)+j-1)*1./j)    ! count interface twice
                      tmp_kk2(i1,j1,k1)=nint(k1*(kk2(nbl)+k-1)*1./k) - nint((k1-1.)*(kk2(nbl)+k-1)*1./k)
                   end do
                end do
             end do
             if (min(minval(tmp_ii2),minval(tmp_jj2),minval(tmp_kk2))>1) then ! if the splitting is acceptable   !  todo : criteria may be different elsewhere
                !             compute sizes of new communication, must be over evaluated (including boundary condition)
                num_cft=2*(tmp_jj2*tmp_ii2 + tmp_ii2*tmp_kk2 + tmp_jj2*tmp_kk2)
                !             choose the best splitting (less comm)
                if(sum(num_cft)<sum(num_cf2)) then  !  TODO : is sum better than maxval ?
                   nblockd=(/i,j,k/)
                   call reallocate(  new_jj2,i,j,k) ;   new_jj2=tmp_jj2
                   call reallocate(  new_ii2,i,j,k) ;   new_ii2=tmp_ii2
                   call reallocate(  new_kk2,i,j,k) ;   new_kk2=tmp_kk2
                   call reallocate(num_cf2,i,j,k) ; num_cf2=num_cft
                end if
             end if
             deallocate(tmp_jj2,tmp_ii2,tmp_kk2,num_cft)
          end if
       end do
    end do
 end do
end subroutine triv_split

subroutine copy_grid(l,l2,xs,ys,zs,x,y,z,save_x,save_y,save_z,       &
    save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn)
 use maillage,only : kk1,kk2,jj1,jj2,ii1,ii2,kd1,kd2,jd1,jd2,id1,id2,npn
 implicit none
 integer,intent(in)             :: l,l2,xs,ys,zs,save_npn(:)
 integer,intent(in)             :: save_id1(:),save_id2(:),save_jd1(:),save_jd2(:),save_kd1(:),save_kd2(:)
 double precision,intent(in)    :: save_x(:),save_y(:),save_z(:)
 double precision,intent(inout) :: x(:),y(:),z(:)

 integer :: xi,yi,zi,nid,njd,nijd,xyz
 integer :: save_xi,save_yi,save_zi,save_nid,save_njd,save_nijd,save_xyz

  nid = id2(l2)-id1(l2)+1
  njd = jd2(l2)-jd1(l2)+1
  nijd = nid*njd

  save_nid = save_id2(l)-save_id1(l)+1
  save_njd = save_jd2(l)-save_jd1(l)+1
  save_nijd = save_nid*save_njd

 do zi=kk1(l2),kk2(l2)
    do yi=jj1(l2),jj2(l2)
       do xi=ii1(l2),ii2(l2)

          save_xi=xi+xs
          save_yi=yi+ys
          save_zi=zi+zs

          xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
          save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd

          ! fill grid
          x(xyz)=save_x(save_xyz)
          y(xyz)=save_y(save_xyz)
          z(xyz)=save_z(save_xyz)
       enddo
    enddo
 enddo
end subroutine copy_grid


subroutine send_grid(x,y,z,xs,ys,zs,xe,ye,ze,       &
    id1,jd1,kd1,nid,nijd,npn,orig,dest)
  use mod_mpi
 implicit none
 integer,intent(in)             :: xs,ys,zs,xe,ye,ze,id1,jd1,kd1,nid,nijd,npn,orig,dest
 double precision,intent(in)    :: x(:),y(:),z(:)
 double precision,allocatable   :: buff(:,:,:,:)
 integer :: xi,yi,zi,xyz

 allocate(buff(3,xe-xs+1,ye-ys+1,xe-xs+1))
 do zi=zs,ze
    do yi=ys,ye
       do xi=xs,xe

          xyz =npn+1+(xi-id1)+(yi-jd1)*nid+(zi-kd1)*nijd

          ! fill buffer
          buff(1,xi-xs+1,yi-ys+1,zi-zs+1)=x(xyz)
          buff(2,xi-xs+1,yi-ys+1,zi-zs+1)=y(xyz)
          buff(3,xi-xs+1,yi-ys+1,zi-zs+1)=z(xyz)
       enddo
    enddo
 enddo

  call MPI_TRANS(buff,buff,orig,DEST)

 deallocate(buff)
end subroutine send_grid

subroutine recv_grid(x,y,z,xs,ys,zs,xe,ye,ze,       &
    id1,jd1,kd1,nid,nijd,npn,orig,dest)
  use mod_mpi
 implicit none
 integer,intent(in)             :: xs,ys,zs,xe,ye,ze,id1,jd1,kd1,nid,nijd,npn,orig,dest
 double precision,intent(inout) :: x(:),y(:),z(:)
 double precision,allocatable   :: buff(:,:,:,:)
 integer :: xi,yi,zi,xyz

 allocate(buff(3,xe-xs+1,ye-ys+1,xe-xs+1))

  call MPI_TRANS(buff,buff,orig,dest)

 do zi=zs,ze
    do yi=ys,ye
       do xi=xs,xe

          xyz      =     npn+1+(     xi-     id1)+(     yi-     jd1)*     nid+(     zi-     kd1)*     nijd

          ! fill grid
          x(xyz)=buff(1,xi-xs+1,yi-ys+1,zi-zs+1)
          y(xyz)=buff(2,xi-xs+1,yi-ys+1,zi-zs+1)
          z(xyz)=buff(3,xi-xs+1,yi-ys+1,zi-zs+1)
       enddo
    enddo
 enddo
 deallocate(buff)
end subroutine recv_grid

subroutine get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,  &
      a,b,c,d,e,f,id1,id2,jd1,jd2,kd1,kd2,npn,x,y,z)
implicit none
  double precision,intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in)           :: a,b,c,d,e,f,id1,id2,jd1,jd2,kd1,kd2,npn
  double precision,intent(in)  :: x(:),y(:),z(:)

  integer nid,njd,nijd

  nid = id2-id1+1
  njd = jd2-jd1+1
  nijd = nid*njd

  xmin=x( npn+1+(a - id1)+(c - jd1)*nid+(e - kd1)*nijd )
  xmax=x( npn+1+(b - id1)+(c - jd1)*nid+(e - kd1)*nijd )
  ymin=y( npn+1+(a - id1)+(c - jd1)*nid+(e - kd1)*nijd )
  ymax=y( npn+1+(a - id1)+(d - jd1)*nid+(e - kd1)*nijd )
  zmin=z( npn+1+(a - id1)+(c - jd1)*nid+(e - kd1)*nijd )
  zmax=z( npn+1+(a - id1)+(c - jd1)*nid+(f - kd1)*nijd )

end subroutine get_coords_box

subroutine iniraccord(mot,imot,nmot)
 use chainecarac,only : ci
 use para_fige,only : nmx
 use boundary,only : tab_raccord
 use mod_valenti
 implicit none
 character(len=32) ::  mot(nmx),comment
 integer             :: imot(nmx),nmot,fr1,fr2,nm,kval,icmt

 do icmt=1,32
    comment(icmt:icmt)=' '
 enddo
 kval=0
 !
 nm=3
 if(nmot.lt.nm) then ! read number first index of coincident boundary
    comment=ci
    call synterr(mot,imot,nmot,comment)
 else
    call valenti(mot,imot,nm,fr1,kval)
 endif
 nm=nm+1
 if(nmot.lt.nm) then ! read number second index of coincident boundary
    comment=ci
    call synterr(mot,imot,nmot,comment)
 else
    call valenti(mot,imot,nm,fr2,kval)
 endif

 ! fr1 and fr2 are coincident boundary
 tab_raccord(fr1)=fr2
 tab_raccord(fr2)=fr1
end subroutine iniraccord

subroutine str(mot,imot,nmx,lmot,val)
 implicit none
 integer          ,intent(in)    :: nmx,lmot,val
 integer          ,intent(inout) :: imot(nmx)
 character(len=32),intent(inout) ::  mot(nmx)

 write(mot(lmot) ,*) val
 mot(lmot)=adjustl(mot(lmot))
 imot(lmot) =len_trim(adjustl(mot(lmot)))
end subroutine str

subroutine get_param(mot,nmot,imot,nblocks)
 use chainecarac,only : ci
 use para_fige,only : nmx
 use mod_valenti
 implicit none
 integer,intent(in)  :: nmot,imot(nmx)
 integer,intent(out) :: nblocks
 character(len=32),intent(in) ::  mot(nmx)
 integer :: icmt,kval,nm
 character(len=32) ::  comment

 ! get number of block we want from flec TODO : replace valenti by mpi
 do icmt=1,32
    comment(icmt:icmt)=' '
 enddo
 kval=0
 !
 nm=2
 if(nmot.lt.nm) then ! read number of block at the end  TODO : replace valenti by mpi
    comment=ci
    call synterr(mot,imot,nmot,comment)
 else
    call valenti(mot,imot,nm,nblocks,kval)
 endif
end subroutine get_param

end module mod_partitionnement
