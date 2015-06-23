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
    use mod_crdms
    use mod_inbdc
    use mod_inbdb
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
    integer             :: icmt,nprocs,nxyza,i,j,k,xyz,nm
    integer             :: imot(nmx),nmot,fr,imax,imin,jmax,jmin,kmax,kmin,kval
    integer             :: l2,mfbe,nid,njd,nijd,xi,yi,zi,sblock,l3,old_mtb
    integer             :: imax2,imin2,jmax2,jmin2,kmax2,kmin2,fr2,i2,j2,k2,l4
    double precision    :: exs1,exs2,vbc(ista*lsta),xmin,xmax,ymin,ymax,zmin,zmax
    double precision    :: save_xmin,save_xmax,save_ymin,save_ymax,save_zmin,save_zmax
    double precision,allocatable :: x(:),y(:),z(:)
    integer,allocatable :: nblock2(:),nblockd(:,:),ni(:,:,:,:),nj(:,:,:,:),nk(:,:,:,:),tmp(:,:,:,:)
    integer,allocatable :: old2new_p(:),new2old_p(:),new2old_b(:),num_cf2(:,:,:)
    integer,allocatable :: ni1(:,:,:),nj1(:,:,:),nk1(:,:,:),ncbd(:)

    integer             :: save_lt,save_ndimntbx,l
    integer             :: save_nid,save_njd,save_nijd,save_klzx
    integer             :: save_kmtbx,save_lzx,save_mdimtbx,save_mdimubx
    integer             :: save_mtb,save_mtbx,save_mtt,save_ndimctbx
    integer             :: save_ndimubx,save_xi,save_xyz,save_yi,save_zi
    double precision,allocatable :: save_x(:),save_y(:),save_z(:)
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
    character(len=2),allocatable :: save_indfl(:)

    character(len=50)::fich
    character(len=32) ::  mot(nmx)
    character(len=32) ::  comment

    integer         ,allocatable ::   mnc(:),ncin(:)
    double precision,allocatable ::  bceqt(:,:)
    logical :: test


    !############################################################################################
    !############################## GET PARAMETERS ##############################################
    !############################################################################################

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

    !############################################################################################
    !############################# SAVE OLD MESH FOR CHECKING PURPOSE ###########################
    !############################################################################################

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


    !############################################################################################
    !################# SAVE ALL MESH AND BOUNDARIES VARIABLES ###################################
    !############################ AND CLEAR EVERYTHING ##########################################
    !####################### WE WILL RECREATE EVERYTHING ########################################
    !############################################################################################

    ! Save old split
    save_lzx=lzx
    save_klzx=klzx

    ! lz=lzx
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
       print*,'découpage réussis : '
       do l=1,save_lt
          print*, l,nblockd(:,l)
       enddo
       !      do l=1,save_lt
       !        do k=1,nblockd(3,l)
       !          do j=1,nblockd(2,l)
       !              print*, l,k,j,ni(:nblockd(1,l),j,k,l)*nj(:nblockd(1,l),j,k,l)*nk(:nblockd(1,l),j,k,l)
       !          enddo
       !        enddo
       !      enddo
    end if


    !############################################################################################
    !######################### RECREATE GRID ####################################################
    !############################################################################################


    allocate(old2new_p(save_ndimntbx),new2old_p(0)) 
    allocate(new2old_b(nprocs))

    do l=1,save_lt  

       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)!     create new grid and initialize it

                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)

                ! duplicate interface

                if(i>1) ni(i,j,k,l)=ni(i,j,k,l)+1 
                if(j>1) nj(i,j,k,l)=nj(i,j,k,l)+1 
                if(k>1) nk(i,j,k,l)=nk(i,j,k,l)+1 

                call crdms( l2,ni(i,j,k,l),nj(i,j,k,l),nk(i,j,k,l))

                call reallocate_s(new2old_p,ndimntbx)
                call reallocate_s(x,ndimntbx)
                call reallocate_s(y,ndimntbx)
                call reallocate_s(z,ndimntbx)

                do zi=kk1(l2),kk2(l2)
                   do yi=jj1(l2),jj2(l2)
                      do xi=ii1(l2),ii2(l2)

                          ! don't forget to count the interface twice
                         save_xi=xi+sum(ni(:i-1,j,k,l))-i+1
                         save_yi=yi+sum(nj(i,:j-1,k,l))-j+1
                         save_zi=zi+sum(nk(i,j,:k-1,l))-k+1

                         nid = id2(l2)-id1(l2)+1
                         njd = jd2(l2)-jd1(l2)+1
                         nijd = nid*njd

                         save_nid = save_id2(l)-save_id1(l)+1
                         save_njd = save_jd2(l)-save_jd1(l)+1
                         save_nijd = save_nid*save_njd

                         xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                         save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd

                         !           fill grid
                         x(xyz)=save_x(save_xyz)
                         y(xyz)=save_y(save_xyz)
                         z(xyz)=save_z(save_xyz)

                         !           fill temporary arrays
                         old2new_p(save_xyz)=    xyz
                         new2old_p(    xyz)=save_xyz
                         new2old_b(l2)= l

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !############################################################################################
    !######################### RECREATE OLD BOUNDARIES ##########################################
    !############################################################################################

    allocate(new2old_f(0))
    do l=1,save_lt  
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                do fr=1,save_mtb
                   if(save_ndlb(fr)==l) then
                      l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                      nid = id2(l2)-id1(l2)+1
                      njd = jd2(l2)-jd1(l2)+1
                      nijd = nid*njd

                      save_nid = save_id2(l)-save_id1(l)+1
                      save_njd = save_jd2(l)-save_jd1(l)+1
                      save_nijd = save_nid*save_njd

                      xmin=x( npn(l2)+1+(ii1(l2)-     id1(l2))+(jj1(l2)-     jd1(l2))*     nid+(kk1(l2)-     kd1(l2))*     nijd )
                      xmax=x( npn(l2)+1+(ii2(l2)-     id1(l2))+(jj1(l2)-     jd1(l2))*     nid+(kk1(l2)-     kd1(l2))*     nijd )
                      ymin=y( npn(l2)+1+(ii1(l2)-     id1(l2))+(jj1(l2)-     jd1(l2))*     nid+(kk1(l2)-     kd1(l2))*     nijd )
                      ymax=y( npn(l2)+1+(ii1(l2)-     id1(l2))+(jj2(l2)-     jd1(l2))*     nid+(kk1(l2)-     kd1(l2))*     nijd )
                      zmin=z( npn(l2)+1+(ii1(l2)-     id1(l2))+(jj1(l2)-     jd1(l2))*     nid+(kk1(l2)-     kd1(l2))*     nijd )
                      zmax=z( npn(l2)+1+(ii1(l2)-     id1(l2))+(jj1(l2)-     jd1(l2))*     nid+(kk2(l2)-     kd1(l2))*     nijd )

                      save_xmin=save_x( save_npn(l)+1+(save_iminb(fr)-save_id1(l))+(save_jminb(fr)-save_jd1(l))*save_nid+(save_kminb(fr)-save_kd1(l))*save_nijd )
                      save_xmax=save_x( save_npn(l)+1+(save_imaxb(fr)-save_id1(l))+(save_jminb(fr)-save_jd1(l))*save_nid+(save_kminb(fr)-save_kd1(l))*save_nijd )
                      save_ymin=save_y( save_npn(l)+1+(save_iminb(fr)-save_id1(l))+(save_jminb(fr)-save_jd1(l))*save_nid+(save_kminb(fr)-save_kd1(l))*save_nijd )
                      save_ymax=save_y( save_npn(l)+1+(save_iminb(fr)-save_id1(l))+(save_jmaxb(fr)-save_jd1(l))*save_nid+(save_kminb(fr)-save_kd1(l))*save_nijd )
                      save_zmin=save_z( save_npn(l)+1+(save_iminb(fr)-save_id1(l))+(save_jminb(fr)-save_jd1(l))*save_nid+(save_kminb(fr)-save_kd1(l))*save_nijd )
                      save_zmax=save_z( save_npn(l)+1+(save_iminb(fr)-save_id1(l))+(save_jminb(fr)-save_jd1(l))*save_nid+(save_kmaxb(fr)-save_kd1(l))*save_nijd )

                      if (save_xmin<=xmax .and. &
                          save_xmax>=xmin .and. &
                          save_ymin<=ymax .and. & ! there is a part of the boundary in this block
                          save_ymax>=ymin .and. &
                          save_zmin<=zmax .and. &
                          save_zmax>=zmin) then

                          ! search indexs

                          xi=ii1(l2)
                          yi=jj1(l2)
                          zi=kk1(l2)
                          save_xi=save_iminb(fr)
                          save_yi=save_jminb(fr)
                          save_zi=save_kminb(fr)


                          do imin=ii1(l2),ii2(l2)-1
                             xi=imin ; save_xi=save_iminb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(x(xyz)>=save_x(save_xyz)) exit
                          enddo
                          do imax=imin,ii2(l2)-1
                             xi=imax ; save_xi=save_imaxb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(x(xyz)>=save_x(save_xyz)) exit
                          enddo
                          do jmin=jj1(l2),jj2(l2)-1
                             yi=jmin ; save_yi=save_jminb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(y(xyz)>=save_y(save_xyz)) exit
                          enddo
                          do jmax=jmin,jj2(l2)-1
                             yi=jmax ; save_yi=save_jmaxb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(y(xyz)>=save_y(save_xyz)) exit
                          enddo
                          do kmin=kk1(l2),kk2(l2)-1
                             zi=kmin ; save_zi=save_kminb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(z(xyz)>=save_z(save_xyz)) exit
                          enddo
                          do kmax=kmin,kk2(l2)-1
                             zi=kmax ; save_zi=save_kmaxb(fr)
                             xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                             save_xyz =save_npn(l )+1+(save_xi-save_id1(l ))+(save_yi-save_jd1(l ))*save_nid+(save_zi-save_kd1(l ))*save_nijd
                             if(z(xyz)>=save_z(save_xyz)) exit
                          enddo

                         if(tab_raccord(fr)/=0) then ! if boundary shared with another block
                            fr2=tab_raccord(fr)      ! a split on the other side induce split here
                            l3=save_ndlb(fr2) 
                            do k2=1,nblockd(3,l3)
                               do j2=1,nblockd(2,l3)
                                  do i2=1,nblockd(1,l3)
                                     l4=sum(nblock2(1:l3-1))+i2+(j2-1)*nblockd(1,l3)+(k2-1)*nblockd(2,l3)*nblockd(1,l3)

                                     save_nid = id2(l4)-id1(l4)+1
                                     save_njd = jd2(l4)-jd1(l4)+1
                                     save_nijd = save_nid*save_njd

                                      save_xmin=x( npn(l4)+1+(ii1(l4)-     id1(l4))+(jj1(l4)-     jd1(l4))*save_nid+(kk1(l4)-     kd1(l4))*save_nijd )
                                      save_xmax=x( npn(l4)+1+(ii2(l4)-     id1(l4))+(jj1(l4)-     jd1(l4))*save_nid+(kk1(l4)-     kd1(l4))*save_nijd )
                                      save_ymin=y( npn(l4)+1+(ii1(l4)-     id1(l4))+(jj1(l4)-     jd1(l4))*save_nid+(kk1(l4)-     kd1(l4))*save_nijd )
                                      save_ymax=y( npn(l4)+1+(ii1(l4)-     id1(l4))+(jj2(l4)-     jd1(l4))*save_nid+(kk1(l4)-     kd1(l4))*save_nijd )
                                      save_zmin=z( npn(l4)+1+(ii1(l4)-     id1(l4))+(jj1(l4)-     jd1(l4))*save_nid+(kk1(l4)-     kd1(l4))*save_nijd )
                                      save_zmax=z( npn(l4)+1+(ii1(l4)-     id1(l4))+(jj1(l4)-     jd1(l4))*save_nid+(kk2(l4)-     kd1(l4))*save_nijd )

                                      if (save_xmin<=xmax .and. &
                                          save_xmax>=xmin .and. &
                                          save_ymin<=ymax .and. & ! there is a part of the boundary in this block
                                          save_ymax>=ymin .and. &
                                          save_zmin<=zmax .and. &
                                          save_zmax>=zmin) then

                                    ! search indexs

                                     xi=imin
                                     yi=jmin
                                     zi=kmin
                                     save_xi=ii1(l4)
                                     save_yi=jj1(l4)
                                     save_zi=kk1(l4)

                                     do imin2=imin,imax-1
                                        xi=imin2 ; save_xi=ii1(l4)
                                        xyz =     npn(l2)+1+(     xi-id1(l2))+(     yi-jd1(l2))*     nid+(     zi-kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(x(xyz)>=x(save_xyz)) exit
                                     enddo
                                     do imax2=imin2,imax-1
                                        xi=imax2 ; save_xi=ii2(l4)
                                        xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(x(xyz)>=x(save_xyz)) exit
                                     enddo
                                     do jmin2=jmin,jmax-1
                                        yi=jmin2 ; save_yi=jj1(l4)
                                        xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(y(xyz)>=y(save_xyz)) exit
                                     enddo
                                     do jmax2=jmin2,jmax-1
                                        yi=jmax2 ; save_yi=jj2(l4)
                                        xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(y(xyz)>=y(save_xyz)) exit
                                     enddo
                                     do kmin2=kmin,kmax-1
                                        zi=kmin2 ; save_zi=kk1(l4)
                                        xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(z(xyz)>=z(save_xyz)) exit
                                     enddo
                                     do kmax2=kmin2,kmax-1
                                        zi=kmax2 ; save_zi=kk2(l4)
                                        xyz =     npn(l2)+1+(     xi-     id1(l2))+(     yi-     jd1(l2))*     nid+(     zi-     kd1(l2))*     nijd
                                        save_xyz =npn(l4)+1+(save_xi-id1(l4))+(save_yi-jd1(l4))*save_nid+(save_zi-kd1(l4))*save_nijd
                                        if(z(xyz)>=z(save_xyz)) exit
                                     enddo

                                        mfbe=mfbe+1
                                        call crbds( &
                                             mfbe,1,l2, &
                                             imin2,imax2,jmin2,jmax2,kmin2,kmax2, &
                                             save_indfl(fr), &
                                             ncbd)
                                        call reallocate_s(new2old_f,mtb)
                                        new2old_f(mfbe)= fr

                                     endif
                                  enddo
                               enddo
                            enddo
                         else
                            mfbe=mfbe+1
                            call crbds( &
                                 mfbe,1,l2, &
                                 imin,imax,jmin,jmax,kmin,kmax, &
                                 save_indfl(fr), &
                                 ncbd)
                            call reallocate_s(new2old_f,mtb)
                            new2old_f(mfbe)= fr
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    old_mtb=mfbe    ! the first mfbe boundaries are old ones

    !############################################################################################
    !########################### CREATE NEW BOUNDARIES ##########################################
    !############################################################################################

    do l=1,save_lt  !           New coincident boundaries
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)

                if (i>1) then
                   l3=sum(nblock2(1:l-1))+i-1+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+2

                   call crbds( &
                        mfbe-1,1,l2, &
                        ii1(l2),ii1(l2),jj1(l2),jj2(l2),kk1(l2),kk2(l2), &
                        'i1', &
                        ncbd)

                   call crbds( &
                        mfbe,1,l3, &
                        ii2(l3),ii2(l3),jj1(l3),jj2(l3),kk1(l3),kk2(l3), &
                        'i2', &
                        ncbd)

                   call reallocate_s(new2old_f,mtb)
                   new2old_f(mfbe-1)= 0
                   new2old_f(mfbe  )= 0
                endif
                if (j>1) then
                   l3=sum(nblock2(1:l-1))+i+(j-2)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+2

                   call crbds( &
                        mfbe-1,1,l2, &
                        ii1(l2),ii2(l2),jj1(l2),jj1(l2),kk1(l2),kk2(l2), &
                        'j1', &
                        ncbd)

                   call crbds( &
                        mfbe,1,l3, &
                        ii1(l3),ii2(l3),jj2(l3),jj2(l3),kk1(l3),kk2(l3), &
                        'j2', &
                        ncbd)

                   call reallocate_s(new2old_f,mtb)
                   new2old_f(mfbe-1)= 0
                   new2old_f(mfbe  )= 0
                endif
                if (k>1) then
                   l3=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-2)*nblockd(2,l)*nblockd(1,l)
                   mfbe=mfbe+2

                   call crbds( &
                        mfbe-1,1,l2, &
                        ii1(l2),ii2(l2),jj1(l2),jj2(l2),kk1(l2),kk1(l2), &
                        'k1', &
                        ncbd)

                   call crbds( &
                        mfbe,1,l3, &
                        ii1(l3),ii2(l3),jj1(l2),jj2(l2),kk2(l2),kk2(l2), &
                        'k2', &
                        ncbd)

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

    !############################################################################################
    !################### INITIALIZE COINCIDENT BOUNDARIES #######################################
    !############################################################################################


    ip21=ndimntbx
    ip40=mdimubx            ! Nb point frontiere
    ip41=mdimtbx            ! Nb point frontiere
    ip42=mdimtbx            ! Nb point frontiere
    ip43=2*mdimtbx            ! Nb point frontiere
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

    do l=1,save_lt
       do k=1,nblockd(3,l)
          do j=1,nblockd(2,l)
             do i=1,nblockd(1,l)
                l2=sum(nblock2(1:l-1))+i+(j-1)*nblockd(1,l)+(k-1)*nblockd(2,l)*nblockd(1,l)
                do fr=1,mfbe
                   test=fr>old_mtb ! new boundary
                   if (new2old_f(fr)/=0) test=test.or.tab_raccord(new2old_f(fr))/=0   ! old raccord boundary
                   if (test) then   ! raccord boundary
                      select case(indfl(fr))
                      case("i1")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr+1,1,0., &
                              ii1(l2),jj1(l2),kk1(l2),'fa','+j','+k')

                      case("i2")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr-1,1,0., &
                              ii2(l3),jj1(l2),kk1(l2),'fa','+j','+k')

                      case("j1")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr+1,1,0., &
                              ii1(l2),jj1(l2),kk1(l2),'+i','fa','+k')

                      case("j2")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr-1,1,0., &
                              ii1(l2),jj2(l2),kk1(l2),'+i','fa','+k')

                      case("k1")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr+1,1,0., &
                              ii1(l2),jj1(l2),kk1(l2),'+i','+j','fa')

                      case("k2")
                         call inbdb( &
                              ncbd,ncin, &
                              fr,"rc  ",1, &
                              0,0,0,0,vbc,bceqt)

                         call inbdc( &
                              exs1,exs2, &
                              x,y,z, &
                              ncbd,ncin,mnc, &
                              0,fr,fr-1,1,0., &
                              ii1(l2),jj1(l2),kk2(l2),'+i','+j','fa')
                      end select
                   endif
                enddo
             enddo
          enddo
       enddo
       print*,l,' : filling done'
    enddo

    return
  contains

  end subroutine partitionnement

  subroutine num_split(nblock2,lt,nxyza,nprocs,ii2,jj2,kk2)
    implicit none
    integer,intent(in)  :: lt,nxyza,nprocs
    integer,intent(in)  :: ii2(lt),jj2(lt),kk2(lt)
    integer,intent(out) :: nblock2(lt)
    integer             :: rsize,sblock(lt),i,j,k

    ! todo
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
                         tmp_jj2(i1,j1,k1)=nint(j1*jj2(nbl)*1./j) - nint((j1-1.)*jj2(nbl)*1./j)
                         tmp_kk2(i1,j1,k1)=nint(k1*kk2(nbl)*1./k) - nint((k1-1.)*kk2(nbl)*1./k)
                      end do
                   end do
                end do
                if (min(minval(tmp_ii2),minval(tmp_jj2),minval(tmp_kk2))>1) then ! if the splitting is acceptable   !  todo : criteria may be different elsewhere
                   !             compute sizes of new communication, must be over evaluated (including boundary condition)
                   num_cft=2*(tmp_jj2*tmp_ii2 + tmp_ii2*tmp_kk2 + tmp_jj2*tmp_kk2)
                   !             choose the best splitting (less comm)
                   if(sum(num_cft)<sum(num_cf2)) then  !  todo : is sum better than maxval ?
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

  subroutine iniraccord(mot,imot,nmot)
    use chainecarac
    use boundary
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
  end subroutine reallocate_s_1r

end module mod_partitionnement



!subroutine partition(nblock,nblockt,nprocs,nxyza,num_cf_all,nicv,njcv,nkcv,num_cf,num_dr,x,y,z, &
!                     num_typ,num_num,num_pr,num_cfi,nijkbl,nxmax,nzmax,li,lk,rank)
!  implicit none
!! todo  simplify !!!!!
!! todo using subroutines
!    real   ,allocatable,intent(inout) :: x(:),y(:),z(:)
!    integer,allocatable,intent(inout) :: nicv(:),njcv(:),nkcv(:),nijkbl(:)
!    integer,allocatable,intent(inout) :: num_cf(:),num_dr(:),num_cfi(:,:),li(:),lk(:)
!    integer,allocatable,intent(inout) :: num_typ(:),num_pr(:),num_num(:)
!    integer,intent(inout)             :: nblock,nblockt,nprocs,num_cf_all,nxyza,nxmax,nzmax,rank

!    real                :: rsize
!    integer             :: nblockd(3),num_bf_all2,num_cf_all2,num_dr_all2
!    integer             :: i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3,l,nblock2(nblockt),nblock4(0:nprocs-1)
!    integer             :: inb,idr,inr,inp,ine,inf,nbl1,nbl2,nbl3,inp1
!    integer             :: nblock3(nprocs),nbl,sblock(nblockt)
!    integer,allocatable :: nicv2(:,:,:),njcv2(:,:,:),nkcv2(:,:,:),num_typ2(:),num_num2(:),   loc2glob2(:,:,:)
!    integer,allocatable :: nibl2(:,:,:),njbl2(:,:,:),nkbl2(:,:,:),nijkbl2(:,:,:),loc2glob1(:,:,:),glob2loc(:,:,:)
!    integer,allocatable ::  num_dr2(:), num_bf(:), num_pr2(:), num_bf2(:,:,:), num_cf2(:,:,:),num_ind2(:,:)
!    integer,allocatable :: num_sdr2(:),num_sbf(:),num_spr2(:),num_sbf2(:)    ,num_scf2(:)    ,num_cfi2(:,:)
!    double precision,allocatable    :: x1(:,:,:)  ,y1(:,:,:)  ,z1(:,:,:)
!    double precision,allocatable    :: x2(:,:,:,:),y2(:,:,:,:),z2(:,:,:,:)
!    logical                 :: test

!!   what are we doing ?

!!   todo    merge geom blocks -> macro-blocks
!!           splitting macro-block by testing all possibility
!!   todo      more splitting possibilities
!!             choose a splitting
!!             compute size of blocks
!!             compute size of new com
!!             if better then keep it
!!           for each new block
!!             compute index
!!             extract the corresponding block of mesh
!!             generate a correspondance table
!!             for each region of boundary
!!               for each point of a region
!!                 if inside the considered new block
!!                   memorize it
!!               if necessary, count the region
!!             for each point on old interface
!!               if inside the considered new block
!!                 memorize it
!!             for each new interfaces
!!               memorize each point 
!!           distributing the blocks
!!           gathering info on the work done by other
!!           distribute mesh
!!           distribute boundary and interface
!!           finishing
!! todo      reorder blocks with metis

!! todo      tests
!! todo        one block - multiples procs : ok
!! todo        multiples block - one procs
!! todo        nblock > nprocs
!! todo        nblock < nprocs

!!print*,'start partitionning'
!    nbl=0
!    if (rank+1<=nblockt) nbl=blockl2blockg(1)
!    nblock4=0

!    allocate(glob2loc(maxval(nibl*njbl*nkbl),nblockt,3)) ; glob2loc=0
!!    allocate(loc2glob1(maxval(nibl*njbl*nkbl),nprocs,3)) ; loc2glob1=0
!!    allocate(loc2glob2(maxval(nibl*njbl*nkbl),nprocs,3))
!!    allocate(loc2glob3(nprocs,maxval(nibl*njbl*nkbl),3))
!    allocate(  njcv2(1,1,1)) ; njcv2=nxyza*nprocs ! useless initial big value
!    allocate(  nicv2(1,1,1)) ; nicv2=nxyza*nprocs ! useless initial big value
!    allocate(  nkcv2(1,1,1)) ; nkcv2=nxyza*nprocs ! useless initial big value
!    nblockd=0

!! todo
!! here we could split user's blocks into elementary pieces
!! and merge them into macro blocks
!! in order to be more independant from user's partition
!! todo

!    nxyza=sum(nicv*njcv*nkcv)  ! total number of points
!    nblock2=1 ! number of splitting

!    call num_split(nblock2,nblockt,nxyza,nprocs,nicv,njcv,nkcv)

!    if (nbl/=0) then 
!!     spliting my block
!      call triv_split(nblock2(nbl),nbl,nxyza,nicv,njcv,nkcv, &
!                      nblockd,num_cf2,nicv2,njcv2,nkcv2)

!! todo
!! here we could use an iterative process to remove impossible split
!! todo
!!                                      __ _
!!                                     |  | |
!!   todo test more possibilities like |--| | (seems usefull when nblock2 not /3)
!!                                     |__|_| (maybe use a recursive algorithme )
!!
!      nblock4(rank)=nblock2(nbl)
!    end if

!!     splitting is done
!    call sum_mpi((/nblockd(1)*nblockd(2)*nblockd(3)/),sblock(1:1))
!    if(sblock(1)/=nprocs) then
!      if(rank==0) print*, "splitting fail",sblock(1),'=/=',nprocs
!      if(rank==0) print*, "do you have enough points per blocks ?",nxyza,nprocs,(nxyza/nprocs)**(1./3.)
!      call end
!    end if

!    if (nbl/=0) then  !     now we have to fill all arrays
!      allocate(  nibl2(nblockd(1),nblockd(2),nblockd(3)), &
!                 njbl2(nblockd(1),nblockd(2),nblockd(3)), &
!                 nkbl2(nblockd(1),nblockd(2),nblockd(3)), &
!               nijkbl2(nblockd(1),nblockd(2),nblockd(3)), &
!               num_bf2(nblockd(1),nblockd(2),nblockd(3)))
!      allocate(num_sbf2(nblock2(nbl)+1),num_scf2(nblock2(nbl)+1),&
!                num_dr2(nblock2(nbl)+1),num_sdr2(nblock2(nbl)+1))

!      njbl2=njcv2+2 ; nibl2=nicv2+2 ; nkbl2=nkcv2+2
!      num_bf2=2*(nicv2*njcv2 + nicv2*nkcv2 + njcv2*nkcv2)

!      num_bf_all2=sum(num_bf2)
!      num_cf_all2=sum(num_cf2)
!      num_dr_all2=sum(num_dr)*(nprocs+1) ! overestimated
!      allocate(num_cfi2(num_cf_all2,8))
!      allocate(num_ind2(num_bf_all2,2))
!      allocate(num_spr2(num_dr_all2))
!      allocate( num_pr2(num_dr_all2))
!      allocate(num_typ2(num_dr_all2))
!      allocate(num_num2(num_dr_all2))

!      allocate(x1(njbl(nbl),nibl(nbl),nkbl(nbl))) ; x1=reshape(x,shape(x1))
!      allocate(y1(njbl(nbl),nibl(nbl),nkbl(nbl))) ; y1=reshape(y,shape(y1))
!      allocate(z1(njbl(nbl),nibl(nbl),nkbl(nbl))) ; z1=reshape(z,shape(z1))
!      allocate(x2(maxval(njbl2),maxval(nibl2),maxval(nkbl2),nblock2(nbl)))
!      allocate(y2(maxval(njbl2),maxval(nibl2),maxval(nkbl2),nblock2(nbl)))
!      allocate(z2(maxval(njbl2),maxval(nibl2),maxval(nkbl2),nblock2(nbl)))

!       num_dr2=0 ;  num_pr2=0 ;  num_cf2=0
!      num_sdr2=0 ; num_spr2=0 ; num_scf2=0
!      num_cfi2=0 ; inr=1
!      ine=0 ; inf=0
!      do k=1,nblockd(3)
!        do i=1,nblockd(2)
!          do j=1,nblockd(1)
!            inb=j+(i-1)*nblockd(1)+(k-1)*nblockd(2)*nblockd(1)

!!           i is also the number of boundaries in the left direction
!!           the last mesh point in each direction is useless (thus n?cv2(j,i,k)+1 and not n?cv2(j,i,k)+2)
!            j1=sum(njcv2(1:j-1,i,k)-1)+j ; j2=j1+njcv2(j,i,k) ; j3=njcv2(j,i,k)+1
!            i1=sum(nicv2(j,1:i-1,k)-1)+i ; i2=i1+nicv2(j,i,k) ; i3=nicv2(j,i,k)+1
!            k1=sum(nkcv2(j,i,1:k-1)-1)+k ; k2=k1+nkcv2(j,i,k) ; k3=nkcv2(j,i,k)+1

!!           fill x,y,z
!            x2(1:j3,1:i3,1:k3,inb)=x1(j1:j2,i1:i2,k1:k2)
!            y2(1:j3,1:i3,1:k3,inb)=y1(j1:j2,i1:i2,k1:k2)
!            z2(1:j3,1:i3,1:k3,inb)=z1(j1:j2,i1:i2,k1:k2)

!!           fill temporary arrays
!            do k3=k1+1,k2
!              do i3=i1+1,i2
!                do j3=j1+1,j2
!                  inp =j3 + &
!                      (i3-1)*njbl(nbl) + &                 ! old inp of the current point
!                      (k3-1)*njbl(nbl)*nibl(nbl)
!                  inp1=j3-j1+1 + &
!                      (i3-i1)*njbl2(j,i,k) + &             ! new inp of the current point
!                      (k3-k1)*njbl2(j,i,k)*nibl2(j,i,k)
!                  glob2loc (inp ,nbl,1:2)=(/inp1,inb/)
!                  if(glob2loc (inp ,nbl,3)==0) glob2loc (inp ,nbl,3)=1
!!                  loc2glob1(inp1,inb,1:3)=(/inp ,nbl,1/)
!                end do
!              end do
!            end do

!!           fill boundary 
!            do idr=num_sdr(nbl)+1,num_sdr(nbl)+num_dr(nbl)
!              test=.false.
!              inr=num_sdr2(inb)+num_dr2(inb)+1
!              do l=num_spr(idr)+1,num_spr(idr)+num_pr(idr) 
!                inp=num_ind(l,1)                           ! point on the boundary
!                k3=floor((mod(inp,nijk)-1.)/nij)+1
!                i3=floor((mod(mod(inp,nijk),nij)-1.)/nj)+1
!                j3=mod(mod(mod(inp,nijk),nij)-1,nj)+1
!                  if(glob2loc (inp ,nbl,2)==inb) then     ! if it is a boundary inside the current block
!                      inp1=j3-j1+1 + &
!                          (i3-i1)*njbl2(j,i,k) + &             ! new inp of the current point
!                          (k3-k1)*njbl2(j,i,k)*nibl2(j,i,k)
!                      test=.true.
!                       ine=ine+1                           ! count it globally
!                       num_pr2(inr)  =num_pr2(inr)+1       ! count it for the current region
!                      num_ind2(ine,2)=num_ind(l,2)         ! the face have not changed
!                      num_ind2(ine,1)=inp1
!                      glob2loc (inp ,nbl,1:3)=(/inp1,inb,1+num_typ(idr)/)
!!                      loc2glob1(inp1,inb,1:3)=(/inp ,nbl,1+num_typ(idr)/)
!                end if
!              end do 
!              if(test) then                        ! this block is concerned by this region
!                 num_dr2(inb)=num_dr2(inb)+1       ! add the region to the block
!                num_typ2(inr)=num_typ(idr)         ! the typ have not changed
!                num_num2(inr)=num_num(idr)         ! the num have not changed
!                num_sdr2(inb+1:)=num_sdr2(inb)+num_dr2(inb)
!                num_spr2(inr+1:)=num_spr2(inr)+num_pr2(inr)
!              end if
!            end do 
!   
!!           fill fictive interface 
!!           old ones
!            do l=num_scf(nbl)+1,num_scf(nbl)+num_cf(nbl)
!              inp=num_cfi(l,1)                           ! point on the boundary
!              k3=floor((mod(inp,nijk)-1.)/nij)+1
!              i3=floor((mod(mod(inp,nijk),nij)-1.)/nj)+1
!              j3=mod(mod(mod(inp,nijk),nij)-1,nj)+1
!                  if(glob2loc (inp ,nbl,2)==inb) then     ! if it is an interface inside the current block
!                    inp1=glob2loc (inp ,nbl,1)            ! new inp of the current point
!                    inf=inf+1                             ! count it globally
!                    num_cf2(j,i,k)=num_cf2(j,i,k)+1       ! count it for the current block
!                    if(inb<nblock2(nbl)) num_scf2(inb+1:)=inf
!                    num_cfi2(inf,2)=num_cfi(l,2)          ! the face have not changed
!                    glob2loc (inp ,nbl,1:3)=(/inp1,inb,-1-num_cfi(l,2)/)
!!                    loc2glob1(inp1,inb,1:3)=(/inp ,nbl,-1-num_cfi(l,2)/)
!                    num_cfi2(inf,1)=j3-j1+1 + &
!                                   (i3-i1)*njbl2(j,i,k) + &       ! new inp of the current point
!                                   (k3-k1)*njbl2(j,i,k)*nibl2(j,i,k)
!                    select case(num_cfi2(inf,2))
!                      case (1); i3=i3-1; case (2); i3=i3+1
!                      case (3); j3=j3-1; case (4); j3=j3+1
!                      case (5); k3=k3-1; case (6); k3=k3+1
!                    end select
!                    num_cfi2(inf,4:5)=-num_cfi(l,4:5)     ! negative value means old inp and inb, to be corrected later
!                    num_cfi2(inf,3)=j3-j1+1 + &
!                                   (i3-i1)*njbl2(j,i,k) + &       ! new inp of the point on the boundary
!                                   (k3-k1)*njbl2(j,i,k)*nibl2(j,i,k)
!              end if
!            end do

!!           fill fictive interface 
!!           new ones
!            if (i>1) &
!              call fict_interf(1,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)
!            if (i<nblockd(2)) &
!              call fict_interf(2,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)
!            if (j>1) &
!              call fict_interf(3,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)
!            if (j<nblockd(1)) &
!              call fict_interf(4,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)
!            if (k>1) &
!              call fict_interf(5,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)
!            if (k<nblockd(3)) &
!              call fict_interf(6,nibl2,njbl2,nkbl2,nbl,&
!                               i,j,k,inb,nblock2,nblockd, &
!                               inf,num_cf2,num_scf2,num_cfi2)

!          end do
!        end do
!      end do
!!print*,'filling done'

!    else ! initialize some variables to avoid unecessary tests
!      allocate(num_cf2(0,0,0))
!      allocate(num_dr2(1),num_pr2(1),num_typ2(1),num_num2(1),num_cfi2(1,8),num_ind2(1,2))
!      njcv2=0 ; num_cf2=0
!      nicv2=0 ; num_dr2=0
!      nkcv2=0 ; num_cfi2=0
!      ine=0
!    end if

!! todo manage situation where a proc rcv multibles blocks

!!print*,'allocating'
!    if (allocated(nibl2)) deallocate(nibl2)
!    if (allocated(njbl2)) deallocate(njbl2)
!    if (allocated(nkbl2)) deallocate(nkbl2)
!    if (allocated(x1)) deallocate(x1  )
!    if (allocated(y1)) deallocate(y1  )
!    if (allocated(z1)) deallocate(z1)

!!   reallocate final arrays
!    call reallocate(nicv,nprocs)    ; call reallocate(nibl,nprocs)    ; call reallocate(num_cf,nprocs)
!    call reallocate(njcv,nprocs)    ; call reallocate(njbl,nprocs)    ; call reallocate(num_bf,nprocs)
!    call reallocate(nkcv,nprocs)    ; call reallocate(nkbl,nprocs)    ; call reallocate(num_dr,nprocs)
!    call reallocate(num_scf,nprocs) ; call reallocate(num_sbf,nprocs) ; call reallocate(num_sdr,nprocs)
!    call reallocate (nijkbl,nprocs) ; call reallocate (nbl_st,nprocs) 

!!print*,'gather splitting'
!!   gather splitting from all cpus
!    call gather(reshape(nicv2,(/nblock4(rank)/)),nicv,nblock4(rank))
!    call gather(reshape(njcv2,(/nblock4(rank)/)),njcv,nblock4(rank))
!    call gather(reshape(nkcv2,(/nblock4(rank)/)),nkcv,nblock4(rank))

!    call gather_p(nblock4(rank),nblock3)

!    if (allocated(nicv2)) deallocate(nicv2)
!    if (allocated(njcv2)) deallocate(njcv2)
!    if (allocated(nkcv2)) deallocate(nkcv2)

!    njbl=njcv+2
!    nibl=nicv+2
!    nkbl=nkcv+2
!    nijkbl=nibl*njbl*nkbl
!    call reallocate(x,nijkbl(rank+1)) ; call reallocate(y,nijkbl(rank+1)) ; call reallocate(z,nijkbl(rank+1))

!!print*,'spread mesh block'
!!   send mesh block to final cpu
!    do nbl1=1,nprocs
!      do nbl2=1,nprocs
!         if(nbl1-sum(nblock3(1:nbl2-1))<=nblock3(nbl2).and. &
!            nbl1-sum(nblock3(1:nbl2-1))>=1) exit
!      end do
!      nbl2=nbl2-1
!      nbl3=nbl1-sum(nblock3(1:nbl2))
!      if(rank==nbl2) then  ! nbl2 have to send block nbl1, locally named nbl3, to nbl1-1
!        call mpi_trans(       x, &
!                      reshape(x2(1:njbl(nbl1),1:nibl(nbl1),1:nkbl(nbl1),nbl3),(/nijkbl(nbl1)/)), &
!                      nbl2,nbl1-1)
!        call mpi_trans(       y, &
!                      reshape(y2(1:njbl(nbl1),1:nibl(nbl1),1:nkbl(nbl1),nbl3),(/nijkbl(nbl1)/)), &
!                      nbl2,nbl1-1)
!        call mpi_trans(       z, &
!                      reshape(z2(1:njbl(nbl1),1:nibl(nbl1),1:nkbl(nbl1),nbl3),(/nijkbl(nbl1)/)), &
!                      nbl2,nbl1-1)
!      else if(rank==nbl1-1) then    ! nbl1-1 have to rcv block nbl1, locally named nbl3, from nbl2
!        call mpi_trans( x, (/0._rk/), nbl2,nbl1-1)
!        call mpi_trans( y, (/0._rk/), nbl2,nbl1-1)
!        call mpi_trans( z, (/0._rk/), nbl2,nbl1-1)
!      end if
!    end do

!!print*,'gather boundaries and interfaces'
!!   gather all boundaries and interfaces
!    call gather(reshape(num_cf2,(/nblock4(rank)/)),num_cf,nblock4(rank))
!    call gather(num_dr2,num_dr,nblock4(rank))
!    num_bf=2*(nicv*njcv + nicv*nkcv + njcv*nkcv)

!    num_bf_all=max(sum(num_bf),1)
!    num_cf_all=max(sum(num_cf),1)
!    num_dr_all=max(sum(num_dr),1)
!    call reallocate(num_cfi,num_cf_all,8) ; num_cfi=0
!    call reallocate(num_ind,num_bf_all,2)
!    call reallocate( num_pr,num_dr_all)
!    call reallocate(num_spr,num_dr_all)
!    call reallocate(num_typ,num_dr_all)
!    call reallocate(num_num,num_dr_all)

!    call gather(num_pr2      ,num_pr      ,sum(num_dr2))
!    call gather(num_typ2     ,num_typ     ,sum(num_dr2))
!    call gather(num_num2     ,num_num     ,sum(num_dr2))

!    do i=1,2
!      call gather(num_ind2(:,i),num_ind(:,i),ine)
!    end do
!    do i=1,8
!      call gather(num_cfi2(:,i),num_cfi(:,i),sum(num_cf2))
!    end do

!!   finalize num_cfi2(l,4:5) for old fictive interfaces using the temporary array
!    do i=1,nblockt
!      call bcast(glob2loc(:,i,:),i-1)
!      glob2loc(:,i,2)=glob2loc(:,i,2)+sum(nblock3(1:i-1))
!    end do

!    do i=1,num_cf_all
!      if(num_cfi(i,4)<0) &
!         num_cfi(i,4:5)=glob2loc(-num_cfi(i,4),-num_cfi(i,5),1:2)
!    end do

!!   compute some index
!    num_sdr=0
!    num_scf=0
!    num_spr=0
!    do i=1,nprocs
!       if (i>1) then
!         num_sdr(i)=num_sdr(i-1)+num_dr(i-1)
!         num_scf(i)=num_scf(i-1)+num_cf(i-1)
!       end if
!      do idr=num_sdr(i)+1,num_sdr(i)+num_dr(i)
!         if(idr>1) num_spr(idr)=num_spr(idr-1)+num_pr(idr-1)
!      end do
!    end do

!    if (nbl/=0) write(*,'(a,i2,a,i2,a,3i3)') 'splitting done',rank,'/',nblock2(nbl),' : ',nblockd

!!print*,'finishing'
!    nblockt=nprocs
!    nblock=1
!    nbl=rank+1
!    call reallocate(blockl2blockg,nblock)
!    call reallocate(block2proc,nprocs)
!    blockl2blockg=rank+1
!    do i=1,nprocs
!      block2proc(i)=i-1
!    end do
!    nbl_st=0
!    nxyza=nijkbl(nbl)
!    nxmax=  maxval(nibl)
!    nzmax=  maxval(nkbl)
!    call reallocate(li,nxmax) ; call reallocate(lk,nzmax)

!!print*,'partitionning done'
!    call setind(1,nbl1)

!    if (allocated(num_typ2)) deallocate(num_typ2)
!    if (allocated(num_num2)) deallocate(num_num2)
!    if (allocated(loc2glob2)) deallocate(loc2glob2)
!    if (allocated(nijkbl2)) deallocate(nijkbl2)
!    if (allocated(loc2glob1)) deallocate(loc2glob1)
!    if (allocated(num_dr2)) deallocate(num_dr2)
!    if (allocated(num_bf)) deallocate( num_bf)
!    if (allocated(num_pr2)) deallocate( num_pr2)
!    if (allocated(num_bf2)) deallocate( num_bf2)
!    if (allocated(num_cf2)) deallocate( num_cf2)
!    if (allocated(num_ind2)) deallocate(num_ind2)
!    if (allocated(num_sdr2)) deallocate(num_sdr2)
!    if (allocated(num_sbf)) deallocate(num_sbf)
!    if (allocated(num_spr2)) deallocate(num_spr2)
!    if (allocated(num_sbf2)) deallocate(num_sbf2    )
!    if (allocated(num_scf2)) deallocate(num_scf2    )
!    if (allocated(num_cfi2)) deallocate(num_cfi2)
!    if (allocated(glob2loc)) deallocate(   glob2loc)

!    if (allocated(x2)) deallocate(x2)
!    if (allocated(y2)) deallocate(y2)
!    if (allocated(z2)) deallocate(z2)

!! todo reorder blocks with metis
!end subroutine partition

!subroutine num_split(nblock2,nblockt,nxyza,nprocs,nicv,njcv,nkcv)
!  implicit none
!    integer,intent(out) :: nblock2(nblockt)
!    integer,intent(in)  :: nicv(nblockt),njcv(nblockt),nkcv(nblockt)
!    integer,intent(in)  :: nblockt,nxyza,nprocs

!    integer             :: rsize,sblock(nblockt),i,j,k

!! todo
!! switch to the alternative version which permit to have ideal blocks size
!! need a criteria to avoid too small block, and need to manage the residual block
!! todo

!!   compute number of spliting of each blocks with the best equilibrium
!    do i=nblockt,nprocs-1                       ! split until nblockt>=nprocs
!      sblock=ceiling(nicv*njcv*nkcv*1./nblock2) ! compute the current size of blocks
!      j=maxloc(sblock,1)                        ! split the first bigest block
!      do k=j+1,nblockt
!        if(sblock(k)==sblock(j) &          ! if more than one bigest block
!        .and.nblock2(k)>nblock2(j)) &      ! split the most splitted
!           j=k
!      end do
!      nblock2(j)=nblock2(j)+1
!    end do


!!   compute number of spliting of each blocks with the ideal equilibrium
!!    rsize=nint(nxyza*1./nprocs)              ! ideal size of a block
!!    nblock2=ceiling(nicv*njcv*nkcv*1./rsize) ! number of split needed
!!    sblock=mod(nicv*njcv*nkcv,rsize)         ! size of the smallest block
!!    do i=1,nblockt
!!      if (sblock(i) <= something) &          ! allow for small imbalance in order to avoid too small blocks
!!           nblock2(i)=nblock2(i)-1
!!    end do

!end subroutine num_split

!subroutine triv_split(nblock2,nbl,nxyza,nicv,njcv,nkcv, &
!                      nblockd,num_cf2,nicv2,njcv2,nkcv2)
!  implicit none
!    integer,allocatable,intent(in)  :: nicv(:),njcv(:),nkcv(:)
!    integer,intent(in)              :: nxyza,nbl,nblock2

!    integer,allocatable,intent(out) :: nicv2(:,:,:),njcv2(:,:,:),nkcv2(:,:,:), num_cf2(:,:,:)
!    integer,intent(out)             :: nblockd(3)

!    integer             :: i,j,k,i1,j1,k1
!    integer,allocatable :: nibl2(:,:,:),njbl2(:,:,:),nkbl2(:,:,:),num_cft(:,:,:)

!!   trivial spliting : divide my block in nblock2 subblock
!!                      test all possiblities constisting in dividing
!!                      i times in the x direction, j times in the y direction and k times in the z direction
!    allocate(num_cf2(1,1,1)) ; num_cf2=nxyza*nblock2 ! useless initial big value
!    do k=1,nblock2
!    do i=1,nblock2
!    do j=1,nblock2
!      if(i*j*k==nblock2) then !           if we get the right number of blocks
!        allocate(njbl2(j,i,k),nibl2(j,i,k),nkbl2(j,i,k), num_cft(j,i,k))
!!       compute sizes of sub-blocks
!        do k1=1,k
!        do i1=1,i
!        do j1=1,j
!          njbl2(j1,i1,k1)=nint(j1*njcv(nbl)*1./j) - nint((j1-1.)*njcv(nbl)*1./j)
!          nibl2(j1,i1,k1)=nint(i1*nicv(nbl)*1./i) - nint((i1-1.)*nicv(nbl)*1./i)
!          nkbl2(j1,i1,k1)=nint(k1*nkcv(nbl)*1./k) - nint((k1-1.)*nkcv(nbl)*1./k)
!        end do
!        end do
!        end do
!        if (min(minval(nibl2),minval(njbl2),minval(nkbl2))>3) then ! if the splitting is acceptable
!!             compute sizes of new communication, must be over evaluated (including boundary condition)
!          num_cft=2*(njbl2*nibl2 + nibl2*nkbl2 + njbl2*nkbl2)
!!             choose the best splitting (less comm)
!          if(sum(num_cft)<sum(num_cf2)) then  !  todo : is sum better than maxval ?
!            nblockd=(/j,i,k/)
!            call reallocate(  njcv2,j,i,k) ;   njcv2=njbl2
!            call reallocate(  nicv2,j,i,k) ;   nicv2=nibl2
!            call reallocate(  nkcv2,j,i,k) ;   nkcv2=nkbl2
!            call reallocate(num_cf2,j,i,k) ; num_cf2=num_cft
!          end if
!        end if
!        deallocate(njbl2,nibl2,nkbl2,num_cft)
!      end if
!    end do
!    end do
!    end do
!end subroutine triv_split

!subroutine fict_interf(face,nibl2,njbl2,nkbl2,nbl,&
!                       i,j,k,inb,nblock2,nblockd, &
!                       inf,num_cf2,num_scf2,num_cfi2)
!  implicit none
!! compute the fictive interface
!! for trivial splitting!
!! walk on all the new fictive interfaces
!! of the same normal (face=1 : normal =x; face=2 : normal =-x, etc)
!! filling the "big" array num_cfi2

!  integer,intent(in)    :: nibl2(:,:,:),njbl2(:,:,:),nkbl2(:,:,:)
!  integer,intent(in)    :: i,j,k,face,nbl
!  integer,intent(in)    :: nblockd(3),inb,nblock2(:)

!  integer,intent(inout) :: inf,num_cf2(:,:,:),num_scf2(:),num_cfi2(:,:)

!  integer               :: i1,j1,k1,i2,j2,k2,i3,j3,k3,inp
!  integer               :: i4,j4,k4,i5,j5,k5,ii,jj,kk

!  i2=2 ; i3=nibl2(j,i,k)-1
!  j2=2 ; j3=njbl2(j,i,k)-1
!  k2=2 ; k3=nkbl2(j,i,k)-1
!  ii=0 ; jj=0 ; kk=0
!  i4=0 ; j4=0 ; k4=0
!  i5=0 ; j5=0 ; k5=0

!  select case(face) ! restrict walking on the surface of interest
!    case (1); i3=2 ; ii=-1 ; case (2); i2=nibl2(j,i,k)-1; ii=1
!    case (3); j3=2 ; jj=-1 ; case (4); j2=njbl2(j,i,k)-1; jj=1
!    case (5); k3=2 ; kk=-1 ; case (6); k2=nkbl2(j,i,k)-1; kk=1
!  end select

!  do k1=k2,k3
!  do i1=i2,i3 ! walk on the surface of interest
!  do j1=j2,j3
!    inp=j1 + &
!       (i1-1)*njbl2(j,i,k) + &                         ! new inp of the current point
!       (k1-1)*njbl2(j,i,k)*nibl2(j,i,k)
!!    inp1=loc2glob1(inp,inb,1)                          ! old inp of the current point

!!    glob2loc (inp1,nbl,1:3)=(/inp ,inb,-1-num_cfi2(inf,2)/)
!!    loc2glob1(inp ,inb,1:3)=(/inp1,nbl,-1-num_cfi2(inf,2)/)

!    select case(face)
!      case (1); i4=1  ; i5=nibl2(j+jj,i+ii,k+kk)-1-i1
!      case (2); i4=-1 ; i5=2-i1
!      case (3); j4=1  ; j5=njbl2(j+jj,i+ii,k+kk)-1-j1
!      case (4); j4=-1 ; j5=2-j1
!      case (5); k4=1  ; k5=nkbl2(j+jj,i+ii,k+kk)-1-k1
!      case (6); k4=-1 ; k5=2-k1
!    end select

!    inf=inf+1                                          ! count it globally
!    num_cf2(j,i,k)=num_cf2(j,i,k)+1                    ! count it for the current block
!    if(inb<nblock2(nbl)) num_scf2(inb+1:)=inf

!    num_cfi2(inf,2)=face                               ! face
!    num_cfi2(inf,1)=inp                                ! inp of the current point
!    num_cfi2(inf,3)=j1  -j4+ &
!                   (i1-1-i4)*njbl2(j,i,k) + &          ! inp of the point on the boundary
!                   (k1-1-k4)*njbl2(j,i,k)*nibl2(j,i,k)
!    num_cfi2(inf,4)=j1  +j5+ &
!                   (i1-1+i5)*njbl2(j+jj,i+ii,k+kk) + & ! inp of the point on the other side
!                   (k1-1+k5)*njbl2(j+jj,i+ii,k+kk)*nibl2(j+jj,i+ii,k+kk)
!    num_cfi2(inf,5)=j  +jj+ &
!                   (i-1+ii)*nblockd(1)+ &              ! inb of the other block
!                   (k-1+kk)*nblockd(2)*nblockd(1) + &
!                    sum(nblock2(1:rank))              ! using global numerotation
!  end do
!  end do
!  end do

!end subroutine fict_interf

