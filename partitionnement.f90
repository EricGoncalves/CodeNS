module mod_partitionnement
  implicit none
  interface reallocate
     module procedure reallocate_1r,reallocate_1i, &
          reallocate_2r,reallocate_2i, &
          reallocate_3r,reallocate_3i, &
          reallocate_4r,reallocate_4i, &
          reallocate_1c
  end interface reallocate
  interface reallocate_s
     module procedure reallocate_s_1r,reallocate_s_4i,reallocate_s_1i
  end interface reallocate_s
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
    integer             :: icmt,nprocs,nxyza,i,j,k,xyz,nm,xs,ys,zs
    integer             :: imot(nmx),nmot,fr,imax,imin,jmax,jmin,kmax,kmin,kval
    integer             :: l2,mfbe,nid,njd,nijd,xi,yi,zi,sblock,l3,old_mtb
    integer             :: imax2,imin2,jmax2,jmin2,kmax2,kmin2,fr2,i2,j2,k2,l4,fr3
    double precision    :: exs1,exs2,vbc(ista*lsta),xmin,xmax,ymin,ymax,zmin,zmax
    double precision    :: save_xmin,save_xmax,save_ymin,save_ymax,save_zmin,save_zmax
    double precision,allocatable :: x(:),y(:),z(:)
    integer,allocatable :: nblock2(:),nblockd(:,:),ni(:,:,:,:),nj(:,:,:,:),nk(:,:,:,:),tmp(:,:,:,:)
    integer,allocatable :: old2new_p(:),new2old_p(:),new2old_b(:),num_cf2(:,:,:)
    integer,allocatable :: ni1(:,:,:),nj1(:,:,:),nk1(:,:,:),ncbd(:)

    integer             :: save_lt,save_ndimntbx,l,verbosity
    integer             :: save_nid,save_njd,save_nijd,save_klzx
    integer             :: save_kmtbx,save_lzx,save_mdimtbx,save_mdimubx
    integer             :: save_mtb,save_mtt,save_ndimctbx
    integer             :: save_ndimubx,save_xi,save_xyz,save_yi,save_zi
    double precision,allocatable :: save_x(:),save_y(:),save_z(:),save_bc(:,:),save_bceqt(:,:)
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
    character(len=2),allocatable :: save_indfl(:)
    character(len=4),allocatable :: save_cl(:)

    character(len=50)::fich
    character(len=32) ::  mot(nmx)
    character(len=32) ::  comment

    integer         ,allocatable ::   mnc(:),ncin(:)
    double precision,allocatable ::  bceqt(:,:)
    logical :: test


    !############################################################################################
    !############################## GET PARAMETERS ##############################################
    !############################################################################################

    verbosity=2 ! from 0 to 3
    call get_param(mot,nmot,imot,nprocs) ! nprocs is the minimum number of blocks we want

    !############################################################################################
    !############################# SAVE OLD MESH FOR CHECKING PURPOSE ###########################
    !############################################################################################

    if(verbosity>=3) then
      ! write old grid
      do l=1,lt
         write(fich,'(A,I0.2,A)') "origmesh_",l,".dat"
         open(42,file=fich,status="replace")

         do k=kk1(l),kk2(l)
            do j=jj1(l),jj2(l)
               do i=ii1(l),ii2(l)

                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),l

               enddo
               write(42,*) ""
            enddo
         enddo
         close(42)
      enddo
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

    ! Save old split
    save_lzx=lzx
    save_klzx=klzx
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
    save_ndimntbx=ndimntbx

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
    save_kmtbx=kmtbx
    save_mtb=mtb
    save_mtt=mtt
    save_mdimubx=mdimubx
    save_mdimtbx=mdimtbx


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
    vbc=0

    !############################################################################################
    !####################### COMPUTE SPLITTING ##################################################
    !############################################################################################

    nxyza=sum(save_ii2*save_jj2*save_kk2)  ! total number of points

    ! nouveaux tableaux
    allocate(nblock2(save_lt),nblockd(3,save_lt),ni(1,1,1,nprocs),nj(1,1,1,nprocs),nk(1,1,1,nprocs))
    nblock2=1    ! initial number of splitting for each existing block
    nblockd=1    ! initial number of splitting for each existing block


    ! routine calculant combiens de fois splitter chaque block
    ! sortie : nblock2
    ! entrée : tout le reste
    call num_split(nblock2,save_lt,nxyza,nprocs,save_ii2,save_jj2,save_kk2)

    do l=1,save_lt

       ! calcule le split pour un block, c'est à dire 
       ! le nombre de découpe par direction            (nblockd)
       ! le nombre de points dans chaque nouveau block (nicv2,njcv2,nkcv2)
       ! une estimation des communications             (num_cf2)
       call triv_split(nblock2(l),l,nxyza,save_ii2 ,save_jj2 ,save_kk2, &
            nblockd(:,l),num_cf2,ni1,nj1,nk1)

       ! ajoute le split avec les splits des autres blocks
       call reallocate_s(ni,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nprocs)
       ni(:size(ni1,1),:size(ni1,2),:size(ni1,3),l)=ni1

       call reallocate_s(nj,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nprocs)
       nj(:size(nj1,1),:size(nj1,2),:size(nj1,3),l)=nj1

       call reallocate_s(nk,maxval(nblockd(1,:)),maxval(nblockd(2,:)),maxval(nblockd(3,:)),nprocs)
       nk(:size(nk1,1),:size(nk1,2),:size(nk1,3),l)=nk1
    enddo

    sblock=sum(nblockd(1,:)*nblockd(2,:)*nblockd(3,:))
    if(sblock/=nprocs) then
       stop 'partitionnement impossible'
    else
       if(verbosity>=1) then
         print*,'découpage réussis : '
         do l=1,save_lt
            print*, l,nblockd(:,l)
         enddo
       endif
       if(verbosity>=2) then
          do l=1,save_lt
            do k=1,nblockd(3,l)
              do j=1,nblockd(2,l)
                  print*, l,k,j,ni(:nblockd(1,l),j,k,l)*nj(:nblockd(1,l),j,k,l)*nk(:nblockd(1,l),j,k,l)
              enddo
            enddo
          enddo
       endif
    end if


    !############################################################################################
    !######################### RECREATE GRID ####################################################
    !############################################################################################
    allocate(old2new_p(save_ndimntbx),new2old_p(0)) 
    allocate(new2old_b(nprocs))

    if(verbosity>=2) then
      mot="" ; nmot=7 ; imot=0
      mot(1)="create" ; imot(1)=6
      mot(2)="dom"    ; imot(2)=3
      mot(3)="st"     ; imot(3)=2
    endif
!    print*,'recreate grid '
    do l=1,save_lt  
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)!     create new grid and initialize it

                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)

                ! +1 for duplicate interface
                if(i>1) ni(i,j,k,l)=ni(i,j,k,l)+1 
                if(j>1) nj(i,j,k,l)=nj(i,j,k,l)+1 ! TODO put this in triv_split
                if(k>1) nk(i,j,k,l)=nk(i,j,k,l)+1 

                if(verbosity>=2) then
                  call str(mot,imot,nmx,4 ,l2)
                  call str(mot,imot,nmx,5 ,ni(i,j,k,l))
                  call str(mot,imot,nmx,6 ,nj(i,j,k,l))
                  call str(mot,imot,nmx,7 ,nk(i,j,k,l))
                  call c_crdms( mot,imot,nmot)
                else
                  call crdms( l2,ni(i,j,k,l),nj(i,j,k,l),nk(i,j,k,l))
                endif

                call reallocate_s(new2old_p,ndimntbx)
                call reallocate_s(x,ndimntbx)
                call reallocate_s(y,ndimntbx)
                call reallocate_s(z,ndimntbx)

                ! offset, don't forget to count the interface twice
                xs=sum(ni(:i-1,j,k,l))-i+1
                ys=sum(nj(i,:j-1,k,l))-j+1
                zs=sum(nk(i,j,:k-1,l))-k+1

                call copy_grid(l,l2,xs,ys,zs,x,y,z,save_x,save_y,save_z,              &
                     save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn,  &
                     old2new_p,new2old_p)

                new2old_b(l2)= l

             enddo
          enddo
       enddo
    enddo

    !############################################################################################
    !######################### RECREATE OLD BOUNDARIES ##########################################
    !############################################################################################
    allocate(new2old_f(0))

    if(verbosity>=2) then
      mot="" ; nmot=13  ; imot=0
      mot(1)="create"   ; imot(1)=6
      mot(2)="boundary" ; imot(2)=8
      mot(3)="st"       ; imot(3)=2
    endif

!    print*,'recreate old boundaries ' ! TODO there might be something simpler with old2new_p and new2old_p
    do fr=1,save_mtb
       l=save_ndlb(fr)
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)

                call get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,             &
                      ii1(l2),ii2(l2),jj1(l2),jj2(l2),kk1(l2),kk2(l2),         &
                      id1(l2),id2(l2),jd1(l2),jd2(l2),kd1(l2),kd2(l2),npn(l2), &
                      x,y,z)

                call get_coords_box(save_xmin,save_xmax,save_ymin,save_ymax,save_zmin,save_zmax,                 &
                      save_iminb(fr),save_imaxb(fr),save_jminb(fr),save_jmaxb(fr),save_kminb(fr),save_kmaxb(fr), &
                      save_id1(l),save_id2(l),save_jd1(l),save_jd2(l),save_kd1(l),save_kd2(l),save_npn(l),       &
                      save_x,save_y,save_z)

                if ( save_xmin-1d-10<=xmax .and. &
                     save_xmax+1d-10>=xmin .and. &
                     save_ymin-1d-10<=ymax .and. & ! there is a part of the boundary in this block
                     save_ymax+1d-10>=ymin .and. &
                     save_zmin-1d-10<=zmax .and. &
                     save_zmax+1d-10>=zmin) then

                   ! search the part of the boundary which concern this block

                    nid = id2(l2)-id1(l2)+1
                    njd = jd2(l2)-jd1(l2)+1
                    nijd = nid*njd

                    save_nid = save_id2(l)-save_id1(l)+1
                    save_njd = save_jd2(l)-save_jd1(l)+1
                    save_nijd = save_nid*save_njd


                   xi=ii1(l2)
                   yi=jj1(l2)
                   zi=kk1(l2)
                   save_xi=save_iminb(fr)
                   save_yi=save_jminb(fr)
                   save_zi=save_kminb(fr)

                   ! walk on the block boundaries to find limits 
                   do imin=ii1(l2)+1,ii2(l2)
                      xi=imin ; save_xi=save_iminb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(x(xyz)>save_x(save_xyz)+1d-10) exit
                   enddo
                   do imax=imin,ii2(l2)
                      xi=imax ; save_xi=save_imaxb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(x(xyz)>save_x(save_xyz)+1d-10) exit
                   enddo
                   do jmin=jj1(l2)+1,jj2(l2)
                      yi=jmin ; save_yi=save_jminb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(y(xyz)>save_y(save_xyz)+1d-10) exit
                   enddo
                   do jmax=jmin,jj2(l2)
                      yi=jmax ; save_yi=save_jmaxb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(y(xyz)>save_y(save_xyz)+1d-10) exit
                   enddo
                   do kmin=kk1(l2)+1,kk2(l2)
                      zi=kmin ; save_zi=save_kminb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(z(xyz)>save_z(save_xyz)+1d-10) exit
                   enddo
                   do kmax=kmin,kk2(l2)
                      zi=kmax ; save_zi=save_kmaxb(fr)
                      xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                      save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                      if(z(xyz)>save_z(save_xyz)+1d-10) exit
                   enddo

                   ! we got out one point after the aim
                   imin=imin-1 ; jmin=jmin-1 ; kmin=kmin-1
                   imax=imax-1 ; jmax=jmax-1 ; kmax=kmax-1

                   if(tab_raccord(fr)/=0) then ! if boundary shared with another block
                      fr2=tab_raccord(fr)      ! a split on the other side induce split here
                      l3=save_ndlb(fr2) 

                      call get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,             &
                            imin,imax,jmin,jmax,kmin,kmax,         &
                            id1(l2),id2(l2),jd1(l2),jd2(l2),kd1(l2),kd2(l2),npn(l2), &
                            x,y,z)

                      do k2=1,nblockd(3,l3)
                         do j2=1,nblockd(2,l3)
                            do i2=1,nblockd(1,l3)
                               l4=sum(nblock2(1:l3-1))+i2+(j2-1)*nblockd(1,l3)+(k2-1)*nblockd(2,l3)*nblockd(1,l3)

                              call get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,             &
                                    ii1(l4),ii2(l4),jj1(l4),jj2(l4),kk1(l4),kk2(l4),         &
                                    id1(l4),id2(l4),jd1(l4),jd2(l4),kd1(l4),kd2(l4),npn(l4), &
                                    x,y,z)

                                  if ( save_xmin-1d-10<=xmax .and. &
                                       save_xmax+1d-10>=xmin .and. &
                                       save_ymin-1d-10<=ymax .and. & ! there is a part of the boundary in this block
                                       save_ymax+1d-10>=ymin .and. &
                                       save_zmin-1d-10<=zmax .and. &
                                       save_zmax+1d-10>=zmin) then

                                  ! search indexs

                               save_nid = id2(l4)-id1(l4)+1
                               save_njd = jd2(l4)-jd1(l4)+1
                               save_nijd = save_nid*save_njd

                                  xi=imin
                                  yi=jmin
                                  zi=kmin
                                  save_xi=ii1(l4)
                                  save_yi=jj1(l4)
                                  save_zi=kk1(l4)

                                   ! walk on the block boundaries to find limits
                                  do imin2=imin+1,imax
                                     xi=imin2 ; save_xi=ii1(l4)
                                     xyz =     npn(l2)+1+(     xi-id1(l2))+(     yi-jd1(l2))*     nid+(     zi-kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(x(xyz)>x(save_xyz)+1d-10) exit
                                  enddo
                                  do imax2=imin2,imax
                                     xi=imax2 ; save_xi=ii2(l4)
                                     xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(x(xyz)>x(save_xyz)+1d-10) exit
                                  enddo
                                  do jmin2=jmin+1,jmax
                                     yi=jmin2 ; save_yi=jj1(l4)
                                     xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(y(xyz)>y(save_xyz)+1d-10) exit
                                  enddo
                                  do jmax2=jmin2,jmax
                                     yi=jmax2 ; save_yi=jj2(l4)
                                     xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(y(xyz)>y(save_xyz)+1d-10) exit
                                  enddo
                                  do kmin2=kmin+1,kmax
                                     zi=kmin2 ; save_zi=kk1(l4)
                                     xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(z(xyz)>z(save_xyz)+1d-10) exit
                                  enddo
                                  do kmax2=kmin2,kmax
                                     zi=kmax2 ; save_zi=kk2(l4)
                                     xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                     save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                     if(z(xyz)>z(save_xyz)+1d-10) exit
                                  enddo

                                 ! we got out one point after the aim
                                  imin2=imin2-1 ; jmin2=jmin2-1 ; kmin2=kmin2-1
                                  imax2=imax2-1 ; jmax2=jmax2-1 ; kmax2=kmax2-1

                                  mfbe=mfbe+1
                                  if(verbosity>=2) then
                                    call str(mot,imot,nmx,4 ,mfbe)
                                    call str(mot,imot,nmx,5 ,1)
                                    call str(mot,imot,nmx,6 ,l2)
                                    call str(mot,imot,nmx,7 ,imin2)
                                    call str(mot,imot,nmx,8 ,imax2)
                                    call str(mot,imot,nmx,9 ,jmin2)
                                    call str(mot,imot,nmx,10,jmax2)
                                    call str(mot,imot,nmx,11,kmin2)
                                    call str(mot,imot,nmx,12,kmax2)
                                    mot(13)=save_indfl(fr) ; imot(13)=2
                                    call c_crbds( mot,imot,nmot, ncbd)
                                  else
                                    call crbds( &
                                         mfbe,1,l2, &
                                         imin2,imax2,jmin2,jmax2,kmin2,kmax2, &
                                         save_indfl(fr), &
                                         ncbd)
                                  endif
                                  call reallocate_s(new2old_f,mtb)
                                  new2old_f(mfbe)= fr
                               endif
                            enddo
                         enddo
                      enddo
                   else
                      mfbe=mfbe+1
                      if(verbosity>=2) then
                        call str(mot,imot,nmx,4 ,mfbe)
                        call str(mot,imot,nmx,5 ,1)
                        call str(mot,imot,nmx,6 ,l2)
                        call str(mot,imot,nmx,7 ,imin)
                        call str(mot,imot,nmx,8 ,imax)
                        call str(mot,imot,nmx,9 ,jmin)
                        call str(mot,imot,nmx,10,jmax)
                        call str(mot,imot,nmx,11,kmin)
                        call str(mot,imot,nmx,12,kmax)
                        mot(13)=save_indfl(fr) ; imot(13)=2
                        call c_crbds( mot,imot,nmot, ncbd)
                      else
                        call crbds( &
                             mfbe,1,l2, &
                             imin,imax,jmin,jmax,kmin,kmax, &
                             save_indfl(fr), &
                             ncbd)
                      endif
                      call reallocate_s(new2old_f,mtb)
                      new2old_f(mfbe)= fr
                   endif
               endif
            enddo
         enddo
      enddo
    enddo

    !############################################################################################
    !########################### CREATE NEW BOUNDARIES ##########################################
    !############################################################################################

!    print*,'create new boundaries '
    do l=1,save_lt  !           New coincident boundaries
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)

                if (i>1) then
                   l3=sum(nblock2(1:l-1))+i-1+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,ii1(l2))
                     call str(mot,imot,nmx,8 ,ii1(l2))
                     call str(mot,imot,nmx,9 ,jj1(l2))
                     call str(mot,imot,nmx,10,jj2(l2))
                     call str(mot,imot,nmx,11,kk1(l2))
                     call str(mot,imot,nmx,12,kk2(l2))
                     mot(13)="i1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                          ii1(l2),ii1(l2),jj1(l2),jj2(l2),kk1(l2),kk2(l2), &
                          'i1', &
                          ncbd)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,ii2(l3))
                     call str(mot,imot,nmx,8 ,ii2(l3))
                     call str(mot,imot,nmx,9 ,jj1(l3))
                     call str(mot,imot,nmx,10,jj2(l3))
                     call str(mot,imot,nmx,11,kk1(l3))
                     call str(mot,imot,nmx,12,kk2(l3))
                     mot(13)="i2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          ii2(l3),ii2(l3),jj1(l3),jj2(l3),kk1(l3),kk2(l3), &
                          'i2', &
                          ncbd)
                    endif

                   call reallocate_s(new2old_f,mtb)
                   new2old_f(mfbe-1)= 0
                   new2old_f(mfbe  )= 0
                endif
                if (j>1) then
                   l3=sum(nblock2(1:l-1))+i+(j-2)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,ii1(l2))
                     call str(mot,imot,nmx,8 ,ii2(l2))
                     call str(mot,imot,nmx,9 ,jj1(l2))
                     call str(mot,imot,nmx,10,jj1(l2))
                     call str(mot,imot,nmx,11,kk1(l2))
                     call str(mot,imot,nmx,12,kk2(l2))
                     mot(13)="j1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                          ii1(l2),ii2(l2),jj1(l2),jj1(l2),kk1(l2),kk2(l2), &
                          'j1', &
                          ncbd)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,ii1(l3))
                     call str(mot,imot,nmx,8 ,ii2(l3))
                     call str(mot,imot,nmx,9 ,jj2(l3))
                     call str(mot,imot,nmx,10,jj2(l3))
                     call str(mot,imot,nmx,11,kk1(l3))
                     call str(mot,imot,nmx,12,kk2(l3))
                     mot(13)="j2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          ii1(l3),ii2(l3),jj2(l3),jj2(l3),kk1(l3),kk2(l3), &
                          'j2', &
                          ncbd)
                    endif

                   call reallocate_s(new2old_f,mtb)
                   new2old_f(mfbe-1)= 0
                   new2old_f(mfbe  )= 0
                endif
                if (k>1) then
                   l3=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-2)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l2)
                     call str(mot,imot,nmx,7 ,ii1(l2))
                     call str(mot,imot,nmx,8 ,ii2(l2))
                     call str(mot,imot,nmx,9 ,jj1(l2))
                     call str(mot,imot,nmx,10,jj2(l2))
                     call str(mot,imot,nmx,11,kk1(l2))
                     call str(mot,imot,nmx,12,kk1(l2))
                     mot(13)="k1" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l2, &
                          ii1(l2),ii2(l2),jj1(l2),jj2(l2),kk1(l2),kk1(l2), &
                          'k1', &
                          ncbd)
                    endif

                   mfbe=mfbe+1
                   if(verbosity>=2) then
                     call str(mot,imot,nmx,4 ,mfbe)
                     call str(mot,imot,nmx,5 ,1)
                     call str(mot,imot,nmx,6 ,l3)
                     call str(mot,imot,nmx,7 ,ii1(l3))
                     call str(mot,imot,nmx,8 ,ii2(l3))
                     call str(mot,imot,nmx,9 ,jj1(l3))
                     call str(mot,imot,nmx,10,jj2(l3))
                     call str(mot,imot,nmx,11,kk2(l3))
                     call str(mot,imot,nmx,12,kk2(l3))
                     mot(13)="k2" ; imot(13)=2
                     call c_crbds( mot,imot,nmot, ncbd)
                   else
                     call crbds( &
                          mfbe,1,l3, &
                          ii1(l3),ii2(l3),jj1(l2),jj2(l2),kk2(l2),kk2(l2), &
                          'k2', &
                          ncbd)
                    endif

                   call reallocate_s(new2old_f,mtb)
                   new2old_f(mfbe-1)= 0
                   new2old_f(mfbe  )= 0
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
         write(fich,'(A,I0.2,A)') "testmesh_",l,".dat"
         open(42,file=fich,status="replace")

         do k=kk1(l),kk2(l)
            do j=jj1(l),jj2(l)
               do i=ii1(l),ii2(l)

                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),l

               enddo
               write(42,*) ""
            enddo
         enddo
         close(42)
      enddo
      write(fich,'(A,I0.2,A)') "testbnd.dat"
      open(42,file=fich,status="replace")
    ! write new boundaries
      do fr=1,mfbe
         do k=kminb(fr),kmaxb(fr)
            do j=jminb(fr),jmaxb(fr)
               do i=iminb(fr),imaxb(fr)
                  l=ndlb(fr)
                  nid = id2(l)-id1(l)+1
                  njd = jd2(l)-jd1(l)+1
                  nijd = nid*njd

                  xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

                  write(42,'(3e11.3,i8)') x(xyz),y(xyz),z(xyz),fr

               enddo
               if(indfl(fr)(1:1)/="i") write(42,*) ""
            enddo
            if(indfl(fr)(1:1)=="i") write(42,*) ""
         enddo
      enddo
      close(42)
    endif

    !############################################################################################
    !################### INITIALIZE COINCIDENT BOUNDARIES #######################################
    !############################################################################################

    ip21=ndimntbx
    ip40=mdimubx            ! Nb point frontiere
    ip41=mdimtbx            ! Nb point frontiere
    ip42=mdimtbx            ! Nb point frontiere
    ip43=mdimtbx            ! Nb point frontiere
    ip44=0!mdimubx            !TODO
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

!    print*,'initialization '
    do l=1,save_lt
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                do fr=1,mfbe
                   if(ndlb(fr)==l2) then
                      test=.false.
                      if (new2old_f(fr)==0) then  ! new boundary 
                         test=.true.
                         if (indfl(fr)(2:2)=="1") fr2=fr+1
                         if (indfl(fr)(2:2)=="2") fr2=fr-1
                      else                        ! old raccord boundary
                         if(tab_raccord(new2old_f(fr))/=0) then
                            test=.true.
                            fr2=tab_raccord(new2old_f(fr)) ! old other boundary number
                            l3=save_ndlb(fr2)              ! old other block number

                            call get_coords_box(xmin,xmax,ymin,ymax,zmin,zmax,                 &
                                  iminb(fr),imaxb(fr),jminb(fr),jmaxb(fr),kminb(fr),kmaxb(fr), &
                                  id1(l2),id2(l2),jd1(l2),jd2(l2),kd1(l2),kd2(l2),npn(l2),       &
                                  x,y,z)

                            find_otherblock:&
                                 do k2=1,nblockd(3,l3)
                            do j2=1,nblockd(2,l3)
                               do i2=1,nblockd(1,l3)
                                  l4=sum(nblock2(1:l3-1))+i2+(j2-1)*nblockd(1,l3)+(k2-1)*nblockd(2,l3)*nblockd(1,l3) ! potential new block number
                                  do fr3=1,mfbe
                                     if(ndlb(fr3)==l4) then  ! potential new boundary number

                                    call get_coords_box(save_xmin,save_xmax,save_ymin,save_ymax,save_zmin,save_zmax,                 &
                                          iminb(fr3),imaxb(fr3),jminb(fr3),jmaxb(fr3),kminb(fr3),kmaxb(fr3), &
                                          id1(l4),id2(l4),jd1(l4),jd2(l4),kd1(l4),kd2(l4),npn(l4),       &
                                          x,y,z)

                                        if (abs(save_xmin-xmin)<=1d-10 .and. &
                                             abs(save_xmax-xmax)<=1d-10 .and. &
                                             abs(save_ymin-ymin)<=1d-10 .and. & ! there is a part of the boundary in this block
                                             abs(save_ymax-ymax)<=1d-10 .and. &
                                             abs(save_zmin-zmin)<=1d-10 .and. &
                                             abs(save_zmax-zmax)<=1d-10) then

                                           fr2=fr3
                                           exit find_otherblock
                                        endif
                                     endif
                                  enddo
                               enddo
                            enddo
                         enddo find_otherblock
                      endif
                   endif
                   if (test) then   ! raccord boundary
                       if(verbosity>=2) then
                         mot="" ; nmot=6   ; imot=0
                         mot(1)="init"     ; imot(1)=4
                         mot(2)="boundary" ; imot(2)=8
                         mot(3)="basic"    ; imot(3)=5
                         call str(mot,imot,nmx,4 ,fr)
                         mot(5)="rc"       ; imot(5)=2
                         call str(mot,imot,nmx,6 ,1)
                         call c_inbdb( mot,imot,nmot,ncbd,ncin,bceqt,partition=.true.)
                       else
                        call inbdb( &
                             ncbd,ncin, &
                             fr,"rc  ",1, &
                             0,0,0,0,vbc,bceqt)
                        endif

                      select case(indfl(fr)(1:1))
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
                         call str(mot,imot,nmx,4 ,fr)
                         mot(5)="frc"      ; imot(5)=3
                         call str(mot,imot,nmx,6 ,fr2)
                         mot(7)="kibdc"    ; imot(7)=5
                         call str(mot,imot,nmx,8 ,1)
                         mot(9)="krr"      ; imot(9)=3
                         call str(mot,imot,nmx,10 ,0)
                         mot(11)="ptc"     ; imot(11)=3
                         call str(mot,imot,nmx,12 ,iminb(fr2))
                         call str(mot,imot,nmx,13 ,jminb(fr2))
                         call str(mot,imot,nmx,14 ,kminb(fr2))
                         mot(15)="dir"     ; imot(15)=3
                         call c_inbdc(  mot,imot,nmot, exs1,exs2, x,y,z, ncbd,ncin,mnc)
                       else
                        call inbdc( &
                             exs1,exs2, &
                             x,y,z, &
                             ncbd,ncin,mnc, &
                             0,fr,fr2,1,0., &
                             iminb(fr2),jminb(fr2),kminb(fr2),mot(16)(1:2),mot(17)(1:2),mot(18)(1:2))
                       endif
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
!    print*,l,' : filling done'
 enddo

 deallocate(save_x,save_y,save_z,save_ndlb,save_nfei,save_indfl,save_mpb,save_mmb,save_ncbd,save_ii1)
 deallocate(save_jj1,save_kk1,save_ii2,save_jj2,save_kk2,save_id1,save_jd1,save_kd1,save_id2,save_jd2)
 deallocate(save_kd2,save_nnn,save_nnc,save_nnfb,save_npn,save_npc,save_npfb)

 return
end subroutine partitionnement

subroutine num_split(nblock2,lt,nxyza,nprocs,ii2,jj2,kk2)
 implicit none
 integer,intent(in)  :: lt,nxyza,nprocs
 integer,intent(in)  :: ii2(lt),jj2(lt),kk2(lt)
 integer,intent(out) :: nblock2(lt)
 integer             :: rsize,sblock(lt),i,j,k

 ! TODO
 ! switch to the alternative version which permit to have ideal blocks size
 ! need a criteria to avoid too small block, and need to manage the residual block
 ! todo

 !   compute number of spliting of each blocks with the best equilibrium
 do i=lt,nprocs-1                       ! split until lt>=nprocs
    sblock=ceiling(ii2*jj2*kk2*1./nblock2) ! compute the current size of blocks
    j=maxloc(sblock,1)                        ! split the first bigest block
    do k=j+1,lt
       if(sblock(k)==sblock(j) &          ! if more than one bigest block
            .and.nblock2(k)>nblock2(j)) &      ! split the most splitted
            j=k
    end do
    nblock2(j)=nblock2(j)+1
 end do


 !   compute number of spliting of each blocks with the ideal equilibrium
 !    rsize=nint(nxyza*1./nprocs)              ! ideal size of a block
 !    nblock2=ceiling(ii2*jj2*kk2*1./rsize) ! number of split needed
 !    sblock=mod(ii2*jj2*kk2,rsize)         ! size of the smallest block
 !    do i=1,lt
 !      if (sblock(i) <= something) &          ! allow for small imbalance in order to avoid too small blocks
 !           nblock2(i)=nblock2(i)-1
 !    end do

end subroutine num_split

subroutine triv_split(nblock2,nbl,nxyza,ii2,jj2,kk2, &
    nblockd,num_cf2,new_ii2,new_jj2,new_kk2)
 implicit none
 integer,allocatable,intent(in)  :: ii2(:),jj2(:),kk2(:)
 integer,intent(in)              :: nxyza,nbl,nblock2

 integer,allocatable,intent(out) :: new_ii2(:,:,:),new_jj2(:,:,:),new_kk2(:,:,:), num_cf2(:,:,:)
 integer,intent(out)             :: nblockd(3)

 integer             :: i,j,k,i1,j1,k1
 integer,allocatable :: tmp_ii2(:,:,:),tmp_jj2(:,:,:),tmp_kk2(:,:,:),num_cft(:,:,:)

 !   trivial spliting : divide my block in nblock2 subblock
 !                      test all possiblities constisting in dividing
 !                      i times in the x direction, j times in the y direction and k times in the z direction
 allocate(num_cf2(1,1,1)) ; num_cf2=nxyza*nblock2 ! useless initial big value
 do k=1,nblock2
    do j=1,nblock2
       do i=1,nblock2
          if(i*j*k==nblock2) then !           if we get the right number of blocks
             allocate(tmp_jj2(i,j,k),tmp_ii2(i,j,k),tmp_kk2(i,j,k), num_cft(i,j,k))
             !       compute sizes of sub-blocks
             do k1=1,k
                do j1=1,j
                   do i1=1,i
                      tmp_ii2(i1,j1,k1)=nint(i1*ii2(nbl)*1./i) - nint((i1-1.)*ii2(nbl)*1./i)
                      tmp_jj2(i1,j1,k1)=nint(j1*jj2(nbl)*1./j) - nint((j1-1.)*jj2(nbl)*1./j) ! TODO count interface twice
                      tmp_kk2(i1,j1,k1)=nint(k1*kk2(nbl)*1./k) - nint((k1-1.)*kk2(nbl)*1./k)
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
    save_id1,save_id2,save_jd1,save_jd2,save_kd1,save_kd2,save_npn,  &
    old2new_p,new2old_p)
 use maillage,only : kk1,kk2,jj1,jj2,ii1,ii2,kd1,kd2,jd1,jd2,id1,id2,npn
 implicit none
 integer,intent(in)             :: l,l2,xs,ys,zs,save_npn(:)
 integer,intent(in)             :: save_id1(:),save_id2(:),save_jd1(:),save_jd2(:),save_kd1(:),save_kd2(:)
 integer,intent(inout)          :: old2new_p(:),new2old_p(:)
 double precision,intent(in)    :: save_x(:),save_y(:),save_z(:)
 double precision,intent(inout) :: x(:),y(:),z(:)

 integer :: xi,yi,zi,nid,njd,nijd,xyz
 integer :: save_xi,save_yi,save_zi,save_nid,save_njd,save_nijd,save_xyz

 do zi=kk1(l2),kk2(l2)
    do yi=jj1(l2),jj2(l2)
       do xi=ii1(l2),ii2(l2)

          save_xi=xi+xs
          save_yi=yi+ys
          save_zi=zi+zs

          nid = id2(l2)-id1(l2)+1
          njd = jd2(l2)-jd1(l2)+1
          nijd = nid*njd

          save_nid = save_id2(l)-save_id1(l)+1
          save_njd = save_jd2(l)-save_jd1(l)+1
          save_nijd = save_nid*save_njd

          xyz      =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
          save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd

          ! fill grid
          x(xyz)=save_x(save_xyz)
          y(xyz)=save_y(save_xyz)
          z(xyz)=save_z(save_xyz)

          ! fill temporary arrays
          old2new_p(save_xyz)=     xyz ! incomplete on interfaces
          new2old_p(     xyz)=save_xyz
       enddo
    enddo
 enddo
end subroutine copy_grid

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

subroutine reallocate_1r(in,size)
 implicit none
 double precision,allocatable,intent(inout) :: in(:)
 integer,intent(in)             :: size

 if(allocated(in)) deallocate(in)
 allocate(in(size))

end subroutine reallocate_1r

subroutine reallocate_2r(in,size1,size2)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:)
 integer,intent(in)             :: size1,size2

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2))

end subroutine reallocate_2r

subroutine reallocate_3r(in,size1,size2,size3)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:,:)
 integer,intent(in)             :: size1,size2,size3

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2,size3))

end subroutine reallocate_3r

subroutine reallocate_4r(in,size1,size2,size3,size4)
 implicit none
 double precision,allocatable,intent(inout) :: in(:,:,:,:)
 integer,intent(in)                :: size1,size2,size3,size4

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2,size3,size4))

end subroutine reallocate_4r

subroutine reallocate_1i(in,size)
 implicit none
 integer,allocatable,intent(inout) :: in(:)
 integer,intent(in)                :: size

 if(allocated(in)) deallocate(in)
 allocate(in(size))

end subroutine reallocate_1i

subroutine reallocate_2i(in,size1,size2)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:)
 integer,intent(in)                :: size1,size2

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2))

end subroutine reallocate_2i

subroutine reallocate_3i(in,size1,size2,size3)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:,:)
 integer,intent(in)                :: size1,size2,size3

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2,size3))

end subroutine reallocate_3i

subroutine reallocate_4i(in,size1,size2,size3,size4)
 implicit none
 integer,allocatable,intent(inout) :: in(:,:,:,:)
 integer,intent(in)                :: size1,size2,size3,size4

 if(allocated(in)) deallocate(in)
 allocate(in(size1,size2,size3,size4))

end subroutine reallocate_4i

subroutine reallocate_1c(in,size)
 implicit none
 character(*),allocatable,intent(inout) :: in(:)
 integer,intent(in)                :: size

 if(allocated(in)) deallocate(in)
 allocate(in(size))

end subroutine reallocate_1c

subroutine reallocate_s_1i(in,newsize)
 implicit none
 integer,allocatable::in(:),tmp(:)
 integer :: newsize
 allocate(tmp(size(in)))
 tmp=in
 call reallocate(in,newsize)
 in=0
 in(1:size(tmp))=tmp
 deallocate(tmp)
end subroutine reallocate_s_1i

subroutine reallocate_s_4i(in,size1,size2,size3,size4)
 implicit none
 integer,allocatable::in(:,:,:,:),tmp(:,:,:,:)
 integer :: size1,size2,size3,size4
 allocate(tmp(size(in,1),size(in,2),size(in,3),size(in,4)))
 tmp=in
 call reallocate(in,size1,size2,size3,size4)
 in=0
 in(:size(tmp,1),:size(tmp,2),:size(tmp,3),:size(tmp,4))=tmp
 deallocate(tmp)
end subroutine reallocate_s_4i

subroutine reallocate_s_1r(in,newsize)
 implicit none
 double precision,allocatable::in(:),tmp(:)
 integer :: newsize
 allocate(tmp(size(in)))
 tmp=in
 call reallocate(in,newsize)
 in=0.
 in(1:size(tmp))=tmp
 deallocate(tmp)
end subroutine reallocate_s_1r


subroutine str(mot,imot,nmx,lmot,val)
 implicit none
 integer          ,intent(in)    :: nmx,lmot,val
 integer          ,intent(inout) :: imot(nmx)
 character(len=32),intent(inout) ::  mot(nmx)

 write(mot(lmot) ,*) val
 mot(lmot)=adjustl(mot(lmot))
 imot(lmot) =len_trim(adjustl(mot(lmot)))
end subroutine str

subroutine get_param(mot,nmot,imot,nprocs)
 use chainecarac,only : ci
 use para_fige,only : nmx
 use mod_valenti
 implicit none
 integer,intent(in)  :: nmot,imot(nmx)
 integer,intent(out) :: nprocs
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
    call valenti(mot,imot,nm,nprocs,kval)
 endif
end subroutine get_param

end module mod_partitionnement
