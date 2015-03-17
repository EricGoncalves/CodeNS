module mod_sch_weno5pond
  implicit none
contains
  subroutine sch_weno5pond(                                        &
       lm,ityprk,                                             &
       u,v,ff,                                                &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz,             &
       equat,                                                 &
       sn,lgsnlt,                                             &
       fxx,fyy,fzz,fxy,fxz,fyz,fex,fey,fez,                   &
       ps,                                                    &
       cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
!                                                                       
!***********************************************************************
!                                                                       
!_DA  DATE_C : janvier 2005 - Eric Goncalves / LEGI                     
!                                                                       
!     ACT                                                               
!_A   Schema WENO de Jiang&Chu, ordre 5 en maillage regulier.           
!_A   Formulation 2D avec schema de Roe. 
!_A   Prise en compte de la metrique du maillage. 
!                                                                       
!***********************************************************************
!                                                                       
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    implicit none
    integer          ::       i,     i1,   i1m1,   i1p1,     i2
    integer          ::    i2m1,   i2m2,     id,   iexp,   ind1
    integer          ::    ind2,isortie, ityprk,      j,     j1
    integer          ::    j1m1,   j1p1,     j2,   j2m1,   j2m2
    integer          ::      jd,      k,     k1,   k1m1,   k1p1
    integer          ::      k2,   k2m1,   k2m2,     kd,   kdir
    integer          ::  lgsnlt,     lm,      m,     m1,      n
    integer          ::     n0c,     n1,    nci,    ncj,    nid
    integer          ::    nijd,   ninc,    njd
    double precision ::    afc,   afd,   afe,   afg,   afh
    double precision ::    afi,    al,    am,  am2i,    ar
    double precision :: beta11,beta12,beta13,beta21,beta22
    double precision :: beta23,beta31,beta32,beta33,beta41
    double precision :: beta42,beta43,beta51,beta52,beta53
    double precision ::    bfc,   bfd,   bfe,   bfg,   bfh
    double precision ::    bfi,   c00,   c01,   c02,    c1
    double precision ::    c10,   c11,   c12,    c2,   c20
    double precision ::    c21,   c22,   c30,   c31,   c32
    double precision ::  cmui1, cmui2, cmuj1, cmuj2,  cnds
    double precision ::    cvi,   cvj,   df1,   df2,   df3
    double precision ::    df4,   df5,   dg1,   dg2,   dg3
    double precision ::    dg4,   dg5,   eps,    f1,   f11
    double precision ::    f12,   f13,    f2,   f21,   f22
    double precision ::    f23,    f3,   f31,   f32,   f33
    double precision ::     f4,   f41,   f42,   f43,    f5
    double precision ::    f51,   f52,   f53,   fc1,   fc2
    double precision ::    fc3,   fc4,   fc5,   fex,   fey
    double precision ::    fez,    ff,   fv2,   fv3,   fv4
    double precision ::    fv5,   fxx,   fxy,   fxz,   fyy
    double precision ::    fyz,   fzz,    g1,   g11,   g12
    double precision ::    g13,   g1m,   g1p,    g2,   g21
    double precision ::    g22,   g23,   g2m,   g2p,    g3
    double precision ::    g31,   g32,   g33,   g3m,   g3p
    double precision ::     g4,   g41,   g42,   g43,    g5
    double precision ::    g51,   g52,   g53,   gc1,   gc2
    double precision ::    gc3,   gc4,   gc5,    gd,   gd1
    double precision ::    gd2,   gv2,   gv3,   gv4,   gv5
    double precision ::     hl,    hm,    hr,    nx,    ny
    double precision ::     nz,   p11,   p12,   p13,   p14
    double precision ::    p15,   p21,   p22,   p23,   p24
    double precision ::    p25,   p31,   p32,   p33,   p34
    double precision ::    p35,   p41,   p42,   p43,   p44
    double precision ::    p45,   p51,   p52,   p53,   p54
    double precision ::    p55,    ps,   q11,   q12,   q13
    double precision ::    q14,   q15,   q1f, q1f1m, q1f1p
    double precision ::  q1f2m, q1f2p, q1f3p,   q21,   q22
    double precision ::    q23,   q24,   q25,   q2f, q2f1m
    double precision ::  q2f1p, q2f2m, q2f2p, q2f3p,   q31
    double precision ::    q32,   q33,   q34,   q35,   q3f
    double precision ::  q3f1m, q3f1p, q3f2m, q3f2p, q3f3p
    double precision ::    q41,   q42,   q43,   q44,   q45
    double precision ::    q4f, q4f1m, q4f1p, q4f2m, q4f2p
    double precision ::  q4f3p,   q51,   q52,   q53,   q54
    double precision ::    q55,   q5f, q5f1m, q5f1p, q5f2m
    double precision ::  q5f2p, q5f3p,   qcx,   qcy,   qcz
    double precision :: rhoami,rhoiam,  rhom, rhomi,   s11
    double precision ::    s12,   s13,   s14,   s21,   s22
    double precision ::    s23,   s24,   s31,   s32,   s33
    double precision ::    s34,   s41,   s42,   s43,   s44
    double precision ::    s51,   s52,   s53,   s54,    sn
    double precision ::     sw,   t11,   t12,   t13,   t14
    double precision ::    t15,   t16,   t21,   t22,   t23
    double precision ::    t24,   t25,   t26,   t31,   t32
    double precision ::    t33,   t34,   t35,   t36,   t41
    double precision ::    t42,   t43,   t44,   t45,   t46
    double precision ::    t51,   t52,   t53,   t54,   t55
    double precision ::    t56,  toxx,  toxy,  toxz,  toyy
    double precision ::   toyz,  tozz,     u,    ul,    um
    double precision ::     ur,     v,    v1,    v4,    v5
    double precision ::  vitm2,    vl,    vm,    vn,    vr
    double precision ::    w11,   w12,   w13,   w14,   w15
    double precision ::    w21,   w22,   w23,   w24,   w25
    double precision ::    w31,   w32,   w33,   w34,   w35
    double precision ::     wl,    wm,    wr,  ww11,  ww12
    double precision ::   ww13,  ww14,  ww15,  ww21,  ww22
    double precision ::   ww23,  ww24,  ww25,  ww31,  ww32
    double precision ::   ww33,  ww34,  ww35
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
    character(len=7 ) :: equat 
    dimension u(ip11,ip60),v(ip11,ip60),ff(ip11,ip60),ps(ip11) 
    dimension sn(lgsnlt,nind,ndir) 
    dimension toxx(ip12),toxy(ip12),toxz(ip12),                       &
         toyy(ip12),toyz(ip12),tozz(ip12),                       &
         qcx (ip12),qcy (ip12),qcz (ip12)           
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21),cvi(ip21),cvj(ip21)
    dimension fxx(ip00),fyy(ip00),fzz(ip00),fxy(ip00),fxz(ip00), &
         fyz(ip00),fex(ip00),fey(ip00),fez(ip00)
!



    n0c=npc(lm) 
    i1=ii1(lm) 
    i2=ii2(lm) 
    j1=jj1(lm) 
    j2=jj2(lm) 
    k1=kk1(lm) 
    k2=kk2(lm) 
!                                                                       
    nid = id2(lm)-id1(lm)+1 
    njd = jd2(lm)-jd1(lm)+1 
    nijd = nid*njd 
!                                                                       
    i1p1=i1+1 
    j1p1=j1+1 
    k1p1=k1+1 
    i2m1=i2-1 
    j2m1=j2-1 
    k2m1=k2-1 
    i2m2=i2-2 
    j2m2=j2-2 
    k2m2=k2-2 
    i1m1=i1-1 
    j1m1=j1-1 
    k1m1=k1-1 
!                                                                       
    nci = inc(1,0,0) 
    ncj = inc(0,1,0)                                                                      
!                                                                       
!     type de maillage=maillage (x,y)                                          
!     activation des sorties                                            
    isortie=0 
!                                                                       
!-----calcul des densites de flux convectifs -------------------------  
!                                                                       
    ind1 = indc(i1m1,j1m1,k1  ) 
    ind2 = indc(i2  ,j2  ,k2m1) 
    do n=ind1,ind2 
       m=n-n0c 
       u(n,1)=0. 
       u(n,2)=0. 
       u(n,3)=0. 
       u(n,4)=0. 
       u(n,5)=0. 
       fxx(m)=v(n,2)*(v(n,2)/v(n,1))+ps(n)-pinfl
       fxy(m)=v(n,3)*(v(n,2)/v(n,1))
       fxz(m)=v(n,4)*(v(n,2)/v(n,1))
       fyy(m)=v(n,3)*(v(n,3)/v(n,1))+ps(n)-pinfl
       fyz(m)=v(n,4)*(v(n,3)/v(n,1))
       fzz(m)=v(n,4)*(v(n,4)/v(n,1))+ps(n)-pinfl
       fex(m)=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
       fey(m)=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
       fez(m)=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
    enddo
!     coefficient pour calculs senseurs beta                            
    c1=13./12. 
    c2=0.25 
!     epsilon petit                                                     
!      eps=1.e-40
    eps=1.e-6
!                                                                       
!***********************************************************************
!                                                                       
!  Calcul du flux numerique par direction suivant les etapes successives
!    1) evaluation des matrices de passage et des valeurs               
!       propres associees a la matrice jacobienne A                     
!       (ces quantites sont evaluees a l'etat moyen de Roe)             
!    2) calculs des flux associes aux variables caracteristiques        
!    3) application de la procedure de reconstruction ENO scalaire      
!       a chaque composante                                             
!    4) calculs des flux convectifs a partir des flux reconstruits      
!    5) ajout des flux visqueux evalues avec schema centre              
!                                                                       
!***********************************************************************
!                                                                       
!------direction i------------------------------------------------------
!                                                                       
    kdir=1 
    ninc=nci 

    do k=k1,k2m1                                                     
       ind1 = indc(i1,j1  ,k)                                          
       ind2 = indc(i1,j2m1,k)                                          
       do n=ind1,ind2,ncj                                              
          m=n-n0c                                                        
          cvi(m-ninc)=cvi(m)
          cvi(m-2*ninc)=cvi(m)
       enddo
    enddo

    do k=k1,k2m1                                                     
       ind1 = indc(i2m1,j1  ,k)                                        
       ind2 = indc(i2m1,j2m1,k)                                        
       do n=ind1,ind2,ncj                                              
          m=n-n0c
          cvi(m+2*ninc)=cvi(m+ninc)
       enddo
    enddo

    do k=k1,k2m1 
       do j=j1,j2m1 
          ind1 = indc(i1  ,j,k) 
          ind2 = indc(i2m2,j,k) 
          do n=ind1,ind2 
             m=n-n0c 
             m1=m+ninc 
             n1=n+ninc 
!        vecteur normal unitaire a la face consideree (face i+1/2)      
             cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+                         &
                  sn(m1,kdir,2)*sn(m1,kdir,2)+                         &
                  sn(m1,kdir,3)*sn(m1,kdir,3))                         
             nx=sn(m1,kdir,1)/cnds 
             ny=sn(m1,kdir,2)/cnds 
             nz=sn(m1,kdir,3)/cnds 
!        calcul des etats gauche et droit                               
             ul=v(n,2)/v(n,1)
             vl=v(n,3)/v(n,1)
             wl=v(n,4)/v(n,1)
             ur=v(n1,2)/v(n1,1)
             vr=v(n1,3)/v(n1,1)
             wr=v(n1,4)/v(n1,1)
             al=sqrt(gam*ps(n )/v(n,1))
             ar=sqrt(gam*ps(n1)/v(n1,1))
             hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
             hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul de etat moyen de Roe
             gd=sqrt(v(n1,1)/v(n,1))
             gd1=1./(1.+gd)
             gd2=gd*gd1
             rhom=sqrt(v(n,1)*v(n1,1))
             rhomi=1./rhom
             um=gd1*ul+gd2*ur
             vm=gd1*vl+gd2*vr
             wm=gd1*wl+gd2*wr
             hm=gd1*hl+gd2*hr
             vitm2=0.5*(um**2+vm**2+wm**2)
             am=sqrt(abs(gam1*(hm-vitm2)))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
!        valeurs propres de la matrice jacobienne des flux
             v1=vn
             v4=vn+am
             v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
             q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
             q12=gam1*um*nx*am2i
             q13=gam1*vm*nx*am2i+nz*rhomi
             q14=gam1*wm*nx*am2i-ny*rhomi
             q15=-gam1*nx*am2i
             q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
             q22=gam1*um*ny*am2i-nz*rhomi
             q23=gam1*vm*ny*am2i
             q24=gam1*wm*ny*am2i+nx*rhomi
             q25=-gam1*ny*am2i
             q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
             q32=gam1*um*nz*am2i+ny*rhomi
             q33=gam1*vm*nz*am2i-nx*rhomi
             q34=gam1*wm*nz*am2i
             q35=-gam1*nz*am2i
             q41=gam1*vitm2*rhoami-vn*rhomi
             q42=nx*rhomi-gam1*um*rhoami
             q43=ny*rhomi-gam1*vm*rhoami
             q44=nz*rhomi-gam1*wm*rhoami
             q45=gam1*rhoami
             q51=gam1*vitm2*rhoami+vn*rhomi
             q52=-nx*rhomi-gam1*um*rhoami
             q53=-ny*rhomi-gam1*vm*rhoami
             q54=-nz*rhomi-gam1*wm*rhoami
             q55=gam1*rhoami
!
             p11=nx
             p12=ny
             p13=nz
             p14=0.5*rhoiam
             p15=0.5*rhoiam
             p21=um*nx
             p22=um*ny-rhom*nz
             p23=um*nz+rhom*ny
             p24=0.5*rhoiam*(um+nx*am)
             p25=0.5*rhoiam*(um-nx*am)
             p31=vm*nx+rhom*nz
             p32=vm*ny
             p33=vm*nz-rhom*nx
             p34=0.5*rhoiam*(vm+ny*am)
             p35=0.5*rhoiam*(vm-ny*am)
             p41=wm*nx-rhom*ny
             p42=wm*ny+rhom*nx
             p43=wm*nz
             p44=0.5*rhoiam*(wm+nz*am)
             p45=0.5*rhoiam*(wm-nz*am)
             p51=vitm2*nx+rhom*(vm*nz-wm*ny)
             p52=vitm2*ny+rhom*(wm*nx-um*nz)
             p53=vitm2*nz+rhom*(um*ny-vm*nx)
             p54=0.5*rhoiam*(hm+am*vn)
             p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-2) a (i+3)      
             q1f2m=q11*v(n-2*ninc,2)+q12*fxx(m-2*ninc)+q13*fxy(m-2*ninc)    &
                  +q14*fxz(m-2*ninc)+q15*fex(m-2*ninc)    
             q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc)      &
                  +q14*fxz(m-ninc)  +q15*fex(m-ninc)      
             q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m)           &
                  +q14*fxz(m)       +q15*fex(m)           
             q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc)      &
                  +q14*fxz(m+ninc)  +q15*fex(m+ninc)      
             q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc)    &
                  +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)    
             q1f3p=q11*v(n+3*ninc,2)+q12*fxx(m+3*ninc)+q13*fxy(m+3*ninc)    &
                  +q14*fxz(m+3*ninc)+q15*fex(m+3*ninc)    
!                                                                       
             q2f2m=q21*v(n-2*ninc,2)+q22*fxx(m-2*ninc)+q23*fxy(m-2*ninc)    &
                  +q24*fxz(m-2*ninc)+q25*fex(m-2*ninc)    
             q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc)      &
                  +q24*fxz(m-ninc)  +q25*fex(m-ninc)      
             q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m)           &
                  +q24*fxz(m)       +q25*fex(m)           
             q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc)      &
                  +q24*fxz(m+ninc)  +q25*fex(m+ninc)      
             q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc)    &
                  +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)    
             q2f3p=q21*v(n+3*ninc,2)+q22*fxx(m+3*ninc)+q23*fxy(m+3*ninc)    &
                  +q24*fxz(m+3*ninc)+q25*fex(m+3*ninc)    
!                                                                       
             q3f2m=q31*v(n-2*ninc,2)+q32*fxx(m-2*ninc)+q33*fxy(m-2*ninc)    &
                  +q34*fxz(m-2*ninc)+q35*fex(m-2*ninc)    
             q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc)      &
                  +q34*fxz(m-ninc)  +q35*fex(m-ninc)      
             q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m)           &
                  +q34*fxz(m)       +q35*fex(m)           
             q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc)      &
                  +q34*fxz(m+ninc)  +q35*fex(m+ninc)      
             q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc)    &
                  +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)    
             q3f3p=q31*v(n+3*ninc,2)+q32*fxx(m+3*ninc)+q33*fxy(m+3*ninc)    &
                  +q34*fxz(m+3*ninc)+q35*fex(m+3*ninc)    
!                                                                       
             q4f2m=q41*v(n-2*ninc,2)+q42*fxx(m-2*ninc)+q43*fxy(m-2*ninc)    &
                  +q44*fxz(m-2*ninc)+q45*fex(m-2*ninc)    
             q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc)      &
                  +q44*fxz(m-ninc)  +q45*fex(m-ninc)      
             q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m)           &
                  +q44*fxz(m)       +q45*fex(m)           
             q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc)      &
                  +q44*fxz(m+ninc)  +q45*fex(m+ninc)      
             q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc)    &
                  +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)    
             q4f3p=q41*v(n+3*ninc,2)+q42*fxx(m+3*ninc)+q43*fxy(m+3*ninc)    &
                  +q44*fxz(m+3*ninc)+q45*fex(m+3*ninc)    
!                                                                       
             q5f2m=q51*v(n-2*ninc,2)+q52*fxx(m-2*ninc)+q53*fxy(m-2*ninc)    &
                  +q54*fxz(m-2*ninc)+q55*fex(m-2*ninc)    
             q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc)      &
                  +q54*fxz(m-ninc)  +q55*fex(m-ninc)      
             q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m)           &
                  +q54*fxz(m)       +q55*fex(m)           
             q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc)      &
                  +q54*fxz(m+ninc)  +q55*fex(m+ninc)      
             q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc)    &
                  +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)    
             q5f3p=q51*v(n+3*ninc,2)+q52*fxx(m+3*ninc)+q53*fxy(m+3*ninc)    &
                  +q54*fxz(m+3*ninc)+q55*fex(m+3*ninc)    
!        coefficients Crj en maillage irregulier
             afg=cmui2(m1)*cvi(m1)
             afh=afg+cmui2(m+2*ninc)*cvi(m+2*ninc)
             afi=afh+cmui2(m+3*ninc)*cvi(m+3*ninc)
             afe=cmui2(m)*cvi(m)
             afd=afe+cmui2(m-ninc)*cvi(m-ninc)
             afc=afd+cmui2(m-2*ninc)*cvi(m-2*ninc)
             c00=afg*afh/((afg+afe)*(afh+afe))
             c02=-afe*afg/((afh+afe)*afh)
             c01=1.-c00-c02
             c10=-afg*afe/(afd*(afd+afg))
             c12=afe*afd/((afg+afe)*(afd+afg))
             c11=1.-c10-c12
             c20=afe*afd/(afc*(afc-afe))
             c22=1.+afe/afc+afe/afd      
             c21=1.-c20-c22
             c30=1.+afg/afh+afg/afi
             c32=afg*afh/(afi*(afi-afg))
             c31=1.-c30-c32
!         c00=1./3.
!         c02=-1./6.
!         c01=1.-c00-c02
!         c10=-1./6.
!         c12=1./3.
!         c11=1.-c10-c12
!         c20=1./3.
!         c22=11./6.
!         c21=1.-c20-c22
!         c30=11./6.
!         c32=1./3.
!         c31=1.-c30-c32
!        calcul des flux d'ordre 3 sur les 3 stencils       
             f11=0.5*(1.+sign(1.,v1))*(q1f2m*c20 +q1f1m*c21 +q1f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     
             f12=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     
             f13=0.5*(1.+sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1p*c22 +q1f2p*c21 +q1f3p*c20)     
!                                                                       
             f21=0.5*(1.+sign(1.,v1))*(q2f2m*c20 +q2f1m*c21 +q2f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     
             f22=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     
             f23=0.5*(1.+sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1p*c22 +q2f2p*c21 +q2f3p*c20)     
!                                                                       
             f31=0.5*(1.+sign(1.,v1))*(q3f2m*c20 +q3f1m*c21 +q3f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     
             f32=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     
             f33=0.5*(1.+sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1p*c22 +q3f2p*c21 +q3f3p*c20)     
!                                                                       
             f41=0.5*(1.+sign(1.,v4))*(q4f2m*c20 +q4f1m*c21 +q4f  *c22)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     
             f42=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     &
                  +0.5*(1.-sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     
             f43=0.5*(1.+sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1p*c22 +q4f2p*c21 +q4f3p*c20)     
!                                                                       
             f51=0.5*(1.+sign(1.,v5))*(q5f2m*c20 +q5f1m*c21 +q5f  *c22)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     
             f52=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     &
                  +0.5*(1.-sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     
             f53=0.5*(1.+sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1p*c22 +q5f2p*c21 +q5f3p*c20) 
!        calcul des senseurs beta (au carre)                            
             iexp=2 
!         iexp=1                                                        
             s11=(q1f3p-2.*q1f2p+q1f1p)**2 
             s12=(q1f2p-2.*q1f1p+q1f  )**2 
             s13=(q1f1p-2.*q1f  +q1f1m)**2 
             s14=(q1f  -2.*q1f1m+q1f2m)**2 
             t11=(q1f3p-4.*q1f2p+3.*q1f1p)**2 
             t12=(q1f2p-4.*q1f1p+3.*q1f  )**2 
             t13=(q1f2p-q1f  )**2 
             t14=(q1f1p-q1f1m)**2 
             t15=(3.*q1f1p-4.*q1f  +q1f1m)**2 
             t16=(3.*q1f  -4.*q1f1m+q1f2m)**2 
             beta11=(0.5*(1.+sign(1.,v1))*(c1*s14+c2*t16)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s13+c2*t15)+eps)**iexp        
             beta12=(0.5*(1.+sign(1.,v1))*(c1*s13+c2*t14)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s12+c2*t13)+eps)**iexp        
             beta13=(0.5*(1.+sign(1.,v1))*(c1*s12+c2*t12)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s11+c2*t11)+eps)**iexp        
!                                                                       
             s21=(q2f3p-2.*q2f2p+q2f1p)**2 
             s22=(q2f2p-2.*q2f1p+q2f  )**2 
             s23=(q2f1p-2.*q2f  +q2f1m)**2 
             s24=(q2f  -2.*q2f1m+q2f2m)**2 
             t21=(q2f3p-4.*q2f2p+3.*q2f1p)**2 
             t22=(q2f2p-4.*q2f1p+3.*q2f  )**2 
             t23=(q2f2p-q2f  )**2 
             t24=(q2f1p-q2f1m)**2 
             t25=(3.*q2f1p-4.*q2f  +q2f1m)**2 
             t26=(3.*q2f  -4.*q2f1m+q2f2m)**2 
             beta21=(0.5*(1.+sign(1.,v1))*(c1*s24+c2*t26)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s23+c2*t25)+eps)**iexp        
             beta22=(0.5*(1.+sign(1.,v1))*(c1*s23+c2*t24)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s22+c2*t23)+eps)**iexp        
             beta23=(0.5*(1.+sign(1.,v1))*(c1*s22+c2*t22)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s21+c2*t21)+eps)**iexp        
!                                                                       
             s31=(q3f3p-2.*q3f2p+q3f1p)**2 
             s32=(q3f2p-2.*q3f1p+q3f  )**2 
             s33=(q3f1p-2.*q3f  +q3f1m)**2 
             s34=(q3f  -2.*q3f1m+q3f2m)**2 
             t31=(q3f3p-4.*q3f2p+3.*q3f1p)**2 
             t32=(q3f2p-4.*q3f1p+3.*q3f  )**2 
             t33=(q3f2p-q3f  )**2 
             t34=(q3f1p-q3f1m)**2 
             t35=(3.*q3f1p-4.*q3f  +q3f1m)**2 
             t36=(3.*q3f  -4.*q3f1m+q3f2m)**2 
             beta31=(0.5*(1.+sign(1.,v1))*(c1*s34+c2*t36)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s33+c2*t35)+eps)**iexp        
             beta32=(0.5*(1.+sign(1.,v1))*(c1*s33+c2*t34)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s32+c2*t33)+eps)**iexp        
             beta33=(0.5*(1.+sign(1.,v1))*(c1*s32+c2*t32)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s31+c2*t31)+eps)**iexp        
!                                                                       
             s41=(q4f3p-2.*q4f2p+q4f1p)**2 
             s42=(q4f2p-2.*q4f1p+q4f  )**2 
             s43=(q4f1p-2.*q4f  +q4f1m)**2 
             s44=(q4f  -2.*q4f1m+q4f2m)**2 
             t41=(q4f3p-4.*q4f2p+3.*q4f1p)**2 
             t42=(q4f2p-4.*q4f1p+3.*q4f  )**2 
             t43=(q4f2p-q4f  )**2 
             t44=(q4f1p-q4f1m)**2 
             t45=(3.*q4f1p-4.*q4f  +q4f1m)**2 
             t46=(3.*q4f  -4.*q4f1m+q4f2m)**2 
             beta41=(0.5*(1.+sign(1.,v4))*(c1*s44+c2*t46)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s43+c2*t45)+eps)**iexp        
             beta42=(0.5*(1.+sign(1.,v4))*(c1*s43+c2*t44)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s42+c2*t43)+eps)**iexp        
             beta43=(0.5*(1.+sign(1.,v4))*(c1*s42+c2*t42)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s41+c2*t41)+eps)**iexp        
!                                                                       
             s51=(q5f3p-2.*q5f2p+q5f1p)**2 
             s52=(q5f2p-2.*q5f1p+q5f  )**2 
             s53=(q5f1p-2.*q5f  +q5f1m)**2 
             s54=(q5f  -2.*q5f1m+q5f2m)**2 
             t51=(q5f3p-4.*q5f2p+3.*q5f1p)**2 
             t52=(q5f2p-4.*q5f1p+3.*q5f  )**2 
             t53=(q5f2p-q5f  )**2 
             t54=(q5f1p-q5f1m)**2 
             t55=(3.*q5f1p-4.*q5f  +q5f1m)**2 
             t56=(3.*q5f  -4.*q5f1m+q5f2m)**2 
             beta51=(0.5*(1.+sign(1.,v5))*(c1*s54+c2*t56)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s53+c2*t55)+eps)**iexp        
             beta52=(0.5*(1.+sign(1.,v5))*(c1*s53+c2*t54)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s52+c2*t53)+eps)**iexp        
             beta53=(0.5*(1.+sign(1.,v5))*(c1*s52+c2*t52)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s51+c2*t51)+eps)**iexp   
!        coefficients gamma en maillage irregulier
             g1p=afg*afh/((afc+afg)*(afc+afh))
             g3p=afd*afc/((afd+afh)*(afc+afh))
             g2p=1.-g1p-g3p
             g1m=afh*afi/((afd+afh)*(afd+afi))
             g3m=afd*afe/((afd+afi)*(afe+afi))
             g2m=1.-g1m-g3m
!         g1p=0.1
!         g3p=0.3
!         g2p=1.-g1p-g3p
!         g1m=0.3
!         g3m=0.1
!         g2m=1.-g1m-g3m
!        calculs des poids wi    
             ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta11)
             ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta12)
             ww31=0.5*(1.+sign(1.,v1))*(g3p/beta13) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta13)
             sw=ww11+ww21+ww31
             w11=ww11/sw
             w21=ww21/sw
             w31=ww31/sw
!         ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
!              +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
!         ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
!              +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
!         ww31m=w31*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w31+w31**2)/(g3p**2+w31*(1.-2.*g3p)) &        
!              +w31*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w31+w31**2)/(g3m**2+w31*(1.-2.*g3m))
!         swm=ww11m+ww21m+ww31m 
!         w11=ww11m/swm 
!         w21=ww21m/swm 
!         w31=ww31m/swm
!
             ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta21)
             ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta22)
             ww32=0.5*(1.+sign(1.,v1))*(g3p/beta23) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta23)
             sw=ww12+ww22+ww32
             w12=ww12/sw
             w22=ww22/sw
             w32=ww32/sw
!         ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
!              +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
!         ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
!              +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
!         ww32m=w32*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w32+w32**2)/(g3p**2+w32*(1.-2.*g3p)) &        
!              +w32*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w32+w32**2)/(g3m**2+w32*(1.-2.*g3m))
!         swm=ww12m+ww22m+ww32m 
!         w12=ww12m/swm 
!         w22=ww22m/swm 
!         w32=ww32m/swm
!
             ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta31)
             ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta32)
             ww33=0.5*(1.+sign(1.,v1))*(g3p/beta33) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta33)
             sw=ww13+ww23+ww33
             w13=ww13/sw
             w23=ww23/sw
             w33=ww33/sw
!         ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
!              +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
!         ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &        
!              +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
!         ww33m=w33*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w33+w33**2)/(g3p**2+w33*(1.-2.*g3p)) &        
!              +w33*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w33+w33**2)/(g3m**2+w33*(1.-2.*g3m))
!         swm=ww13m+ww23m+ww33m 
!         w13=ww13m/swm 
!         w23=ww23m/swm 
!         w33=ww33m/swm
!
             ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                  +0.5*(1.-sign(1.,v4))*(g1m/beta41)
             ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                  +0.5*(1.-sign(1.,v4))*(g2m/beta42)
             ww34=0.5*(1.+sign(1.,v4))*(g3p/beta43) &
                  +0.5*(1.-sign(1.,v4))*(g3m/beta43)
             sw=ww14+ww24+ww34
             w14=ww14/sw
             w24=ww24/sw
             w34=ww34/sw
!         ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
!              +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
!         ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &    
!              +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
!         ww34m=w34*0.5*(1.+sign(1.,v4))*(g3p+g3p**2-3.*g3p*w34+w34**2)/(g3p**2+w34*(1.-2.*g3p)) &        
!              +w34*0.5*(1.-sign(1.,v4))*(g3m+g3m**2-3.*g3m*w34+w34**2)/(g3m**2+w34*(1.-2.*g3m))
!         swm=ww14m+ww24m+ww34m 
!         w14=ww14m/swm 
!         w24=ww24m/swm 
!         w34=ww34m/swm
!
             ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                  +0.5*(1.-sign(1.,v5))*(g1m/beta51)
             ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                  +0.5*(1.-sign(1.,v5))*(g2m/beta52)
             ww35=0.5*(1.+sign(1.,v5))*(g3p/beta53) &
                  +0.5*(1.-sign(1.,v5))*(g3m/beta53)
             sw=ww15+ww25+ww35
             w15=ww15/sw
             w25=ww25/sw
             w35=ww35/sw
!         ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
!              +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
!         ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &        
!              +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
!         ww35m=w35*0.5*(1.+sign(1.,v5))*(g3p+g3p**2-3.*g3p*w35+w35**2)/(g3p**2+w35*(1.-2.*g3p)) &        
!              +w35*0.5*(1.-sign(1.,v5))*(g3m+g3m**2-3.*g3m*w35+w35**2)/(g3m**2+w35*(1.-2.*g3m))
!         swm=ww15m+ww25m+ww35m
!         w15=ww15m/swm 
!         w25=ww25m/swm 
!         w35=ww35m/swm
!        calcul des flux convectifs projetes
             fc1=w11*f11+w21*f12+w31*f13
             fc2=w12*f21+w22*f22+w32*f23
             fc3=w13*f31+w23*f32+w33*f33
             fc4=w14*f41+w24*f42+w34*f43
             fc5=w15*f51+w25*f52+w35*f53
!        produit avec matrice P pour retour dans l'espace physique
             f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
             f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
             f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
             f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-------------------------------------------------------------------    
!        produit de Q avec les flux Euler aux points (i-2) a (i+3)      
             q1f2m=q11*v(n-2*ninc,3)+q12*fxy(m-2*ninc)+q13*fyy(m-2*ninc)    &
                  +q14*fyz(m-2*ninc)+q15*fey(m-2*ninc)    
             q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc)      &
                  +q14*fyz(m-ninc)  +q15*fey(m-ninc)      
             q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m)           &
                  +q14*fyz(m)       +q15*fey(m)           
             q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc)      &
                  +q14*fyz(m+ninc)  +q15*fey(m+ninc)      
             q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc)    &
                  +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)    
             q1f3p=q11*v(n+3*ninc,3)+q12*fxy(m+3*ninc)+q13*fyy(m+3*ninc)    &
                  +q14*fyz(m+3*ninc)+q15*fey(m+3*ninc)    
!                                                                       
             q2f2m=q21*v(n-2*ninc,3)+q22*fxy(m-2*ninc)+q23*fyy(m-2*ninc)    &
                  +q24*fyz(m-2*ninc)+q25*fey(m-2*ninc)    
             q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc)      &
                  +q24*fyz(m-ninc)  +q25*fey(m-ninc)      
             q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m)           &
                  +q24*fyz(m)       +q25*fey(m)           
             q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc)      &
                  +q24*fyz(m+ninc)  +q25*fey(m+ninc)      
             q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc)    &
                  +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)    
             q2f3p=q21*v(n+3*ninc,3)+q22*fxy(m+3*ninc)+q23*fyy(m+3*ninc)    &
                  +q24*fyz(m+3*ninc)+q25*fey(m+3*ninc)    
!                                                                       
             q3f2m=q31*v(n-2*ninc,3)+q32*fxy(m-2*ninc)+q33*fyy(m-2*ninc)    &
                  +q34*fyz(m-2*ninc)+q35*fey(m-2*ninc)    
             q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc)      &
                  +q34*fyz(m-ninc)  +q35*fey(m-ninc)      
             q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m)           &
                  +q34*fyz(m)       +q35*fey(m)           
             q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc)      &
                  +q34*fyz(m+ninc)  +q35*fey(m+ninc)      
             q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc)    &
                  +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)    
             q3f3p=q31*v(n+3*ninc,3)+q32*fxy(m+3*ninc)+q33*fyy(m+3*ninc)    &
                  +q34*fyz(m+3*ninc)+q35*fey(m+3*ninc)    
!                                                                       
             q4f2m=q41*v(n-2*ninc,3)+q42*fxy(m-2*ninc)+q43*fyy(m-2*ninc)    &
                  +q44*fyz(m-2*ninc)+q45*fey(m-2*ninc)    
             q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc)      &
                  +q44*fyz(m-ninc)  +q45*fey(m-ninc)      
             q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m)           &
                  +q44*fyz(m)       +q45*fey(m)           
             q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc)      &
                  +q44*fyz(m+ninc)  +q45*fey(m+ninc)      
             q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc)    &
                  +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)    
             q4f3p=q41*v(n+3*ninc,3)+q42*fxy(m+3*ninc)+q43*fyy(m+3*ninc)    &
                  +q44*fyz(m+3*ninc)+q45*fey(m+3*ninc)    
!                                                                       
             q5f2m=q51*v(n-2*ninc,3)+q52*fxy(m-2*ninc)+q53*fyy(m-2*ninc)    &
                  +q54*fyz(m-2*ninc)+q55*fey(m-2*ninc)    
             q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc)      &
                  +q54*fyz(m-ninc)  +q55*fey(m-ninc)      
             q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m)           &
                  +q54*fyz(m)       +q55*fey(m)           
             q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc)      &
                  +q54*fyz(m+ninc)  +q55*fey(m+ninc)      
             q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc)    &
                  +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)    
             q5f3p=q51*v(n+3*ninc,3)+q52*fxy(m+3*ninc)+q53*fyy(m+3*ninc)    &
                  +q54*fyz(m+3*ninc)+q55*fey(m+3*ninc)    
!        calcul des flux d'ordre 3 sur les 3 stencils                   
             g11=0.5*(1.+sign(1.,v1))*(q1f2m*c20 +q1f1m*c21 +q1f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     
             g12=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     
             g13=0.5*(1.+sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1p*c22 +q1f2p*c21 +q1f3p*c20)     
!                                                                       
             g21=0.5*(1.+sign(1.,v1))*(q2f2m*c20 +q2f1m*c21 +q2f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     
             g22=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     
             g23=0.5*(1.+sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1p*c22 +q2f2p*c21 +q2f3p*c20)     
!                                                                       
             g31=0.5*(1.+sign(1.,v1))*(q3f2m*c20 +q3f1m*c21 +q3f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     
             g32=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     
             g33=0.5*(1.+sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1p*c22 +q3f2p*c21 +q3f3p*c20)     
!                                                                       
             g41=0.5*(1.+sign(1.,v4))*(q4f2m*c20 +q4f1m*c21 +q4f  *c22)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     
             g42=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     &
                  +0.5*(1.-sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     
             g43=0.5*(1.+sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1p*c22 +q4f2p*c21 +q4f3p*c20)     
!                                                                       
             g51=0.5*(1.+sign(1.,v5))*(q5f2m*c20 +q5f1m*c21 +q5f  *c22)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     
             g52=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     &
                  +0.5*(1.-sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     
             g53=0.5*(1.+sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1p*c22 +q5f2p*c21 +q5f3p*c20)     
!        calcul des senseurs beta (au carre)                            
             iexp=2 
!         iexp=1                                                        
             s11=(q1f3p-2.*q1f2p+q1f1p)**2 
             s12=(q1f2p-2.*q1f1p+q1f  )**2 
             s13=(q1f1p-2.*q1f  +q1f1m)**2 
             s14=(q1f  -2.*q1f1m+q1f2m)**2 
             t11=(q1f3p-4.*q1f2p+3.*q1f1p)**2 
             t12=(q1f2p-4.*q1f1p+3.*q1f  )**2 
             t13=(q1f2p-q1f  )**2 
             t14=(q1f1p-q1f1m)**2 
             t15=(3.*q1f1p-4.*q1f  +q1f1m)**2 
             t16=(3.*q1f  -4.*q1f1m+q1f2m)**2 
             beta11=(0.5*(1.+sign(1.,v1))*(c1*s14+c2*t16)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s13+c2*t15)+eps)**iexp        
             beta12=(0.5*(1.+sign(1.,v1))*(c1*s13+c2*t14)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s12+c2*t13)+eps)**iexp        
             beta13=(0.5*(1.+sign(1.,v1))*(c1*s12+c2*t12)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s11+c2*t11)+eps)**iexp        
!                                                                       
             s21=(q2f3p-2.*q2f2p+q2f1p)**2 
             s22=(q2f2p-2.*q2f1p+q2f  )**2 
             s23=(q2f1p-2.*q2f  +q2f1m)**2 
             s24=(q2f  -2.*q2f1m+q2f2m)**2 
             t21=(q2f3p-4.*q2f2p+3.*q2f1p)**2 
             t22=(q2f2p-4.*q2f1p+3.*q2f  )**2 
             t23=(q2f2p-q2f  )**2 
             t24=(q2f1p-q2f1m)**2 
             t25=(3.*q2f1p-4.*q2f  +q2f1m)**2 
             t26=(3.*q2f  -4.*q2f1m+q2f2m)**2 
             beta21=(0.5*(1.+sign(1.,v1))*(c1*s24+c2*t26)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s23+c2*t25)+eps)**iexp        
             beta22=(0.5*(1.+sign(1.,v1))*(c1*s23+c2*t24)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s22+c2*t23)+eps)**iexp        
             beta23=(0.5*(1.+sign(1.,v1))*(c1*s22+c2*t22)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s21+c2*t21)+eps)**iexp        
!                                                                       
             s31=(q3f3p-2.*q3f2p+q3f1p)**2 
             s32=(q3f2p-2.*q3f1p+q3f  )**2 
             s33=(q3f1p-2.*q3f  +q3f1m)**2 
             s34=(q3f  -2.*q3f1m+q3f2m)**2 
             t31=(q3f3p-4.*q3f2p+3.*q3f1p)**2 
             t32=(q3f2p-4.*q3f1p+3.*q3f  )**2 
             t33=(q3f2p-q3f  )**2 
             t34=(q3f1p-q3f1m)**2 
             t35=(3.*q3f1p-4.*q3f  +q3f1m)**2 
             t36=(3.*q3f  -4.*q3f1m+q3f2m)**2 
             beta31=(0.5*(1.+sign(1.,v1))*(c1*s34+c2*t36)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s33+c2*t35)+eps)**iexp        
             beta32=(0.5*(1.+sign(1.,v1))*(c1*s33+c2*t34)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s32+c2*t33)+eps)**iexp        
             beta33=(0.5*(1.+sign(1.,v1))*(c1*s32+c2*t32)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s31+c2*t31)+eps)**iexp        
!                                                                       
             s41=(q4f3p-2.*q4f2p+q4f1p)**2 
             s42=(q4f2p-2.*q4f1p+q4f  )**2 
             s43=(q4f1p-2.*q4f  +q4f1m)**2 
             s44=(q4f  -2.*q4f1m+q4f2m)**2 
             t41=(q4f3p-4.*q4f2p+3.*q4f1p)**2 
             t42=(q4f2p-4.*q4f1p+3.*q4f  )**2 
             t43=(q4f2p-q4f  )**2 
             t44=(q4f1p-q4f1m)**2 
             t45=(3.*q4f1p-4.*q4f  +q4f1m)**2 
             t46=(3.*q4f  -4.*q4f1m+q4f2m)**2 
             beta41=(0.5*(1.+sign(1.,v4))*(c1*s44+c2*t46)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s43+c2*t45)+eps)**iexp        
             beta42=(0.5*(1.+sign(1.,v4))*(c1*s43+c2*t44)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s42+c2*t43)+eps)**iexp        
             beta43=(0.5*(1.+sign(1.,v4))*(c1*s42+c2*t42)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s41+c2*t41)+eps)**iexp        
!                                                                       
             s51=(q5f3p-2.*q5f2p+q5f1p)**2 
             s52=(q5f2p-2.*q5f1p+q5f  )**2 
             s53=(q5f1p-2.*q5f  +q5f1m)**2 
             s54=(q5f  -2.*q5f1m+q5f2m)**2 
             t51=(q5f3p-4.*q5f2p+3.*q5f1p)**2 
             t52=(q5f2p-4.*q5f1p+3.*q5f  )**2 
             t53=(q5f2p-q5f  )**2 
             t54=(q5f1p-q5f1m)**2 
             t55=(3.*q5f1p-4.*q5f  +q5f1m)**2 
             t56=(3.*q5f  -4.*q5f1m+q5f2m)**2 
             beta51=(0.5*(1.+sign(1.,v5))*(c1*s54+c2*t56)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s53+c2*t55)+eps)**iexp        
             beta52=(0.5*(1.+sign(1.,v5))*(c1*s53+c2*t54)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s52+c2*t53)+eps)**iexp        
             beta53=(0.5*(1.+sign(1.,v5))*(c1*s52+c2*t52)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s51+c2*t51)+eps)**iexp  
!        calculs des poids wi    
             ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta11)
             ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta12)
             ww31=0.5*(1.+sign(1.,v1))*(g3p/beta13) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta13)
             sw=ww11+ww21+ww31
             w11=ww11/sw
             w21=ww21/sw
             w31=ww31/sw
!         ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
!              +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
!         ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
!              +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
!         ww31m=w31*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w31+w31**2)/(g3p**2+w31*(1.-2.*g3p)) &        
!              +w31*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w31+w31**2)/(g3m**2+w31*(1.-2.*g3m))
!         swm=ww11m+ww21m+ww31m 
!         w11=ww11m/swm 
!         w21=ww21m/swm 
!         w31=ww31m/swm
!
             ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta21)
             ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta22)
             ww32=0.5*(1.+sign(1.,v1))*(g3p/beta23) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta23)
             sw=ww12+ww22+ww32
             w12=ww12/sw
             w22=ww22/sw
             w32=ww32/sw
!         ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
!              +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
!         ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
!              +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
!         ww32m=w32*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w32+w32**2)/(g3p**2+w32*(1.-2.*g3p)) &        
!              +w32*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w32+w32**2)/(g3m**2+w32*(1.-2.*g3m))
!         swm=ww12m+ww22m+ww32m 
!         w12=ww12m/swm 
!         w22=ww22m/swm 
!         w32=ww32m/swm
!
             ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta31)
             ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta32)
             ww33=0.5*(1.+sign(1.,v1))*(g3p/beta33) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta33)
             sw=ww13+ww23+ww33
             w13=ww13/sw
             w23=ww23/sw
             w33=ww33/sw
!         ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
!              +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
!         ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &        
!              +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
!         ww33m=w33*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w33+w33**2)/(g3p**2+w33*(1.-2.*g3p)) &        
!              +w33*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w33+w33**2)/(g3m**2+w33*(1.-2.*g3m))
!         swm=ww13m+ww23m+ww33m 
!         w13=ww13m/swm 
!         w23=ww23m/swm 
!         w33=ww33m/swm
!
             ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                  +0.5*(1.-sign(1.,v4))*(g1m/beta41)
             ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                  +0.5*(1.-sign(1.,v4))*(g2m/beta42)
             ww34=0.5*(1.+sign(1.,v4))*(g3p/beta43) &
                  +0.5*(1.-sign(1.,v4))*(g3m/beta43)
             sw=ww14+ww24+ww34
             w14=ww14/sw
             w24=ww24/sw
             w34=ww34/sw
!         ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
!              +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
!         ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &    
!              +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
!         ww34m=w34*0.5*(1.+sign(1.,v4))*(g3p+g3p**2-3.*g3p*w34+w34**2)/(g3p**2+w34*(1.-2.*g3p)) &        
!              +w34*0.5*(1.-sign(1.,v4))*(g3m+g3m**2-3.*g3m*w34+w34**2)/(g3m**2+w34*(1.-2.*g3m))
!         swm=ww14m+ww24m+ww34m 
!         w14=ww14m/swm 
!         w24=ww24m/swm 
!         w34=ww34m/swm
!
             ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                  +0.5*(1.-sign(1.,v5))*(g1m/beta51)
             ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                  +0.5*(1.-sign(1.,v5))*(g2m/beta52)
             ww35=0.5*(1.+sign(1.,v5))*(g3p/beta53) &
                  +0.5*(1.-sign(1.,v5))*(g3m/beta53)
             sw=ww15+ww25+ww35
             w15=ww15/sw
             w25=ww25/sw
             w35=ww35/sw
!         ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
!              +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
!         ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &        
!              +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
!         ww35m=w35*0.5*(1.+sign(1.,v5))*(g3p+g3p**2-3.*g3p*w35+w35**2)/(g3p**2+w35*(1.-2.*g3p)) &        
!              +w35*0.5*(1.-sign(1.,v5))*(g3m+g3m**2-3.*g3m*w35+w35**2)/(g3m**2+w35*(1.-2.*g3m))
!         swm=ww15m+ww25m+ww35m
!         w15=ww15m/swm 
!         w25=ww25m/swm 
!         w35=ww35m/swm
!        calcul des flux convectifs projetes
             gc1=w11*g11+w21*g12+w31*g13
             gc2=w12*g21+w22*g22+w32*g23
             gc3=w13*g31+w23*g32+w33*g33
             gc4=w14*g41+w24*g42+w34*g43
             gc5=w15*g51+w25*g52+w35*g53
!        produit avec matrice P pour retour dans l'espace physique
             g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
             g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
             g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
             g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
             df1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
             df2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
             df3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
             df5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        calcul des flux visqueux (multiplies par -2)
             fv2=(cmui2(m1)*toxx(n)+cmui1(m1)*toxx(n1))*sn(m1,kdir,1) &
                  +(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,2)
             fv3=(cmui2(m1)*toxy(n)+cmui1(m1)*toxy(n1))*sn(m1,kdir,1) &
                  +(cmui2(m1)*toyy(n)+cmui1(m1)*toyy(n1))*sn(m1,kdir,2)  
             fv5=(cmui2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                  +cmui1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                  +(cmui2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                  +cmui1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
             u(n1,1)=u(n1,1)-df1 
             u(n1,2)=u(n1,2)-df2+0.5*fv2 
             u(n1,3)=u(n1,3)-df3+0.5*fv3 
             u(n1,5)=u(n1,5)-df5+0.5*fv5 
             u(n,1)=u(n,1)+df1 
             u(n,2)=u(n,2)+df2-0.5*fv2 
             u(n,3)=u(n,3)+df3-0.5*fv3 
             u(n,5)=u(n,5)+df5-0.5*fv5 
          enddo
       enddo
    enddo
!
!------direction j----------------------------------------------
!
    kdir=2
    ninc=ncj

    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
       do n=ind1,ind2
          m=n-n0c
          cvj(m-ninc)=cvj(m)
          cvj(m-2*ninc)=cvj(m)
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2m1,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2
          m=n-n0c
          cvj(m+2*ninc)=cvj(m+ninc)
       enddo
    enddo

    do k=k1,k2m1 
       do j=j1,j2m2 
          ind1 = indc(i1  ,j,k) 
          ind2 = indc(i2m1,j,k) 
          do n=ind1,ind2 
             m=n-n0c 
             m1=m+ninc 
             n1=n+ninc 
!        vecteur normal unitaire a la face consideree (face j+1/2)      
             cnds=sqrt(sn(m1,kdir,1)*sn(m1,kdir,1)+ &
                  sn(m1,kdir,2)*sn(m1,kdir,2)+ &
                  sn(m1,kdir,3)*sn(m1,kdir,3))
             nx=sn(m1,kdir,1)/cnds
             ny=sn(m1,kdir,2)/cnds
             nz=sn(m1,kdir,3)/cnds
!        calcul des etats gauche et droit
             ul=v(n,2)/v(n,1)
             vl=v(n,3)/v(n,1)
             wl=v(n,4)/v(n,1)
             ur=v(n1,2)/v(n1,1)
             vr=v(n1,3)/v(n1,1)
             wr=v(n1,4)/v(n1,1)
             al=sqrt(gam*ps(n )/v(n,1))
             ar=sqrt(gam*ps(n1)/v(n1,1))
             hl=al*al/gam1+0.5*(ul**2+vl**2+wl**2)
             hr=ar*ar/gam1+0.5*(ur**2+vr**2+wr**2)
!        calcul des etats moyens de Roe
             gd=sqrt(v(n1,1)/v(n,1))
             gd1=1./(1.+gd)
             gd2=gd*gd1
             rhom=sqrt(v(n,1)*v(n1,1))
             rhomi=1./rhom
             um=gd1*ul+gd2*ur
             vm=gd1*vl+gd2*vr
             wm=gd1*wl+gd2*wr
             hm=gd1*hl+gd2*hr
             vitm2=0.5*(um**2+vm**2+wm**2)
             am=sqrt(abs(gam1*(hm-vitm2)))
             am2i=1./(am*am)
             vn=um*nx+vm*ny+wm*nz
             rhoiam=rhom/am
             rhoami=am2i/rhoiam
!        valeurs propres
             v1=vn
             v4=vn+am
             v5=vn-am
!        calcul des matrices de passage a gauche Q et a droite P
             q11=(1.-gam1*vitm2*am2i)*nx-(vm*nz-wm*ny)*rhomi
             q12=gam1*um*nx*am2i
             q13=gam1*vm*nx*am2i+nz*rhomi
             q14=gam1*wm*nx*am2i-ny*rhomi
             q15=-gam1*nx*am2i
             q21=(1.-gam1*vitm2*am2i)*ny-(wm*nx-um*nz)*rhomi
             q22=gam1*um*ny*am2i-nz*rhomi
             q23=gam1*vm*ny*am2i
             q24=gam1*wm*ny*am2i+nx*rhomi
             q25=-gam1*ny*am2i
             q31=(1.-gam1*vitm2*am2i)*nz-(um*ny-vm*nx)*rhomi
             q32=gam1*um*nz*am2i+ny*rhomi
             q33=gam1*vm*nz*am2i-nx*rhomi
             q34=gam1*wm*nz*am2i
             q35=-gam1*nz*am2i
             q41=gam1*vitm2*rhoami-vn*rhomi
             q42=nx*rhomi-gam1*um*rhoami
             q43=ny*rhomi-gam1*vm*rhoami
             q44=nz*rhomi-gam1*wm*rhoami
             q45=gam1*rhoami
             q51=gam1*vitm2*rhoami+vn*rhomi
             q52=-nx*rhomi-gam1*um*rhoami
             q53=-ny*rhomi-gam1*vm*rhoami
             q54=-nz*rhomi-gam1*wm*rhoami
             q55=gam1*rhoami
!
             p11=nx
             p12=ny
             p13=nz
             p14=0.5*rhoiam
             p15=0.5*rhoiam
             p21=um*nx
             p22=um*ny-rhom*nz
             p23=um*nz+rhom*ny
             p24=0.5*rhoiam*(um+nx*am)
             p25=0.5*rhoiam*(um-nx*am)
             p31=vm*nx+rhom*nz
             p32=vm*ny
             p33=vm*nz-rhom*nx
             p34=0.5*rhoiam*(vm+ny*am)
             p35=0.5*rhoiam*(vm-ny*am)
             p41=wm*nx-rhom*ny
             p42=wm*ny+rhom*nx
             p43=wm*nz
             p44=0.5*rhoiam*(wm+nz*am)
             p45=0.5*rhoiam*(wm-nz*am)
             p51=vitm2*nx+rhom*(vm*nz-wm*ny)
             p52=vitm2*ny+rhom*(wm*nx-um*nz)
             p53=vitm2*nz+rhom*(um*ny-vm*nx)
             p54=0.5*rhoiam*(hm+am*vn)
             p55=0.5*rhoiam*(hm-am*vn)
!        produit de Q avec les flux Euler aux points (i-2) a (i+3)      
             q1f2m=q11*v(n-2*ninc,2)+q12*fxx(m-2*ninc)+q13*fxy(m-2*ninc)    &
                  +q14*fxz(m-2*ninc)+q15*fex(m-2*ninc)    
             q1f1m=q11*v(n-ninc  ,2)+q12*fxx(m-ninc)  +q13*fxy(m-ninc)      &
                  +q14*fxz(m-ninc)  +q15*fex(m-ninc)      
             q1f  =q11*v(n       ,2)+q12*fxx(m)       +q13*fxy(m)           &
                  +q14*fxz(m)       +q15*fex(m)           
             q1f1p=q11*v(n+ninc  ,2)+q12*fxx(m+ninc)  +q13*fxy(m+ninc)      &
                  +q14*fxz(m+ninc)  +q15*fex(m+ninc)      
             q1f2p=q11*v(n+2*ninc,2)+q12*fxx(m+2*ninc)+q13*fxy(m+2*ninc)    &
                  +q14*fxz(m+2*ninc)+q15*fex(m+2*ninc)    
             q1f3p=q11*v(n+3*ninc,2)+q12*fxx(m+3*ninc)+q13*fxy(m+3*ninc)    &
                  +q14*fxz(m+3*ninc)+q15*fex(m+3*ninc)    
!                                                                       
             q2f2m=q21*v(n-2*ninc,2)+q22*fxx(m-2*ninc)+q23*fxy(m-2*ninc)    &
                  +q24*fxz(m-2*ninc)+q25*fex(m-2*ninc)    
             q2f1m=q21*v(n-ninc  ,2)+q22*fxx(m-ninc)  +q23*fxy(m-ninc)      &
                  +q24*fxz(m-ninc)  +q25*fex(m-ninc)      
             q2f  =q21*v(n       ,2)+q22*fxx(m)       +q23*fxy(m)           &
                  +q24*fxz(m)       +q25*fex(m)           
             q2f1p=q21*v(n+ninc  ,2)+q22*fxx(m+ninc)  +q23*fxy(m+ninc)      &
                  +q24*fxz(m+ninc)  +q25*fex(m+ninc)      
             q2f2p=q21*v(n+2*ninc,2)+q22*fxx(m+2*ninc)+q23*fxy(m+2*ninc)    &
                  +q24*fxz(m+2*ninc)+q25*fex(m+2*ninc)    
             q2f3p=q21*v(n+3*ninc,2)+q22*fxx(m+3*ninc)+q23*fxy(m+3*ninc)    &
                  +q24*fxz(m+3*ninc)+q25*fex(m+3*ninc)    
!                                                                       
             q3f2m=q31*v(n-2*ninc,2)+q32*fxx(m-2*ninc)+q33*fxy(m-2*ninc)    &
                  +q34*fxz(m-2*ninc)+q35*fex(m-2*ninc)    
             q3f1m=q31*v(n-ninc  ,2)+q32*fxx(m-ninc)  +q33*fxy(m-ninc)      &
                  +q34*fxz(m-ninc)  +q35*fex(m-ninc)      
             q3f  =q31*v(n       ,2)+q32*fxx(m)       +q33*fxy(m)           &
                  +q34*fxz(m)       +q35*fex(m)           
             q3f1p=q31*v(n+ninc  ,2)+q32*fxx(m+ninc)  +q33*fxy(m+ninc)      &
                  +q34*fxz(m+ninc)  +q35*fex(m+ninc)      
             q3f2p=q31*v(n+2*ninc,2)+q32*fxx(m+2*ninc)+q33*fxy(m+2*ninc)    &
                  +q34*fxz(m+2*ninc)+q35*fex(m+2*ninc)    
             q3f3p=q31*v(n+3*ninc,2)+q32*fxx(m+3*ninc)+q33*fxy(m+3*ninc)    &
                  +q34*fxz(m+3*ninc)+q35*fex(m+3*ninc)    
!                                                                       
             q4f2m=q41*v(n-2*ninc,2)+q42*fxx(m-2*ninc)+q43*fxy(m-2*ninc)    &
                  +q44*fxz(m-2*ninc)+q45*fex(m-2*ninc)    
             q4f1m=q41*v(n-ninc  ,2)+q42*fxx(m-ninc)  +q43*fxy(m-ninc)      &
                  +q44*fxz(m-ninc)  +q45*fex(m-ninc)      
             q4f  =q41*v(n       ,2)+q42*fxx(m)       +q43*fxy(m)           &
                  +q44*fxz(m)       +q45*fex(m)           
             q4f1p=q41*v(n+ninc  ,2)+q42*fxx(m+ninc)  +q43*fxy(m+ninc)      &
                  +q44*fxz(m+ninc)  +q45*fex(m+ninc)      
             q4f2p=q41*v(n+2*ninc,2)+q42*fxx(m+2*ninc)+q43*fxy(m+2*ninc)    &
                  +q44*fxz(m+2*ninc)+q45*fex(m+2*ninc)    
             q4f3p=q41*v(n+3*ninc,2)+q42*fxx(m+3*ninc)+q43*fxy(m+3*ninc)    &
                  +q44*fxz(m+3*ninc)+q45*fex(m+3*ninc)    
!                                                                       
             q5f2m=q51*v(n-2*ninc,2)+q52*fxx(m-2*ninc)+q53*fxy(m-2*ninc)    &
                  +q54*fxz(m-2*ninc)+q55*fex(m-2*ninc)    
             q5f1m=q51*v(n-ninc  ,2)+q52*fxx(m-ninc)  +q53*fxy(m-ninc)      &
                  +q54*fxz(m-ninc)  +q55*fex(m-ninc)      
             q5f  =q51*v(n       ,2)+q52*fxx(m)       +q53*fxy(m)           &
                  +q54*fxz(m)       +q55*fex(m)           
             q5f1p=q51*v(n+ninc  ,2)+q52*fxx(m+ninc)  +q53*fxy(m+ninc)      &
                  +q54*fxz(m+ninc)  +q55*fex(m+ninc)      
             q5f2p=q51*v(n+2*ninc,2)+q52*fxx(m+2*ninc)+q53*fxy(m+2*ninc)    &
                  +q54*fxz(m+2*ninc)+q55*fex(m+2*ninc)    
             q5f3p=q51*v(n+3*ninc,2)+q52*fxx(m+3*ninc)+q53*fxy(m+3*ninc)    &
                  +q54*fxz(m+3*ninc)+q55*fex(m+3*ninc)   
!        coefficients Crj en maillage irregulier
             bfg=cmuj2(m1)*cvj(m1)
             bfh=bfg+cmuj2(m+2*ninc)*cvj(m+2*ninc)
             bfi=bfh+cmuj2(m+3*ninc)*cvj(m+3*ninc)
             bfe=cmuj2(m)*cvj(m)
             bfd=bfe+cmuj2(m-ninc)*cvj(m-ninc)
             bfc=bfd+cmuj2(m-2*ninc)*cvj(m-2*ninc)
             c00=bfg*bfh/((bfg+bfe)*(bfh+bfe))
             c02=-bfe*bfg/((bfh+bfe)*bfh)
             c01=1.-c00-c02
             c10=-bfg*bfe/(bfd*(bfd+bfg))
             c12=bfe*bfd/((bfg+bfe)*(bfd+bfg))
             c11=1.-c10-c12
             c20=bfe*bfd/(bfc*(bfc-bfe))
             c22=1.+bfe/bfc+bfe/bfd      
             c21=1.-c20-c22
             c30=1.+bfg/bfh+bfg/bfi
             c32=bfg*bfh/(bfi*(bfi-bfg))
             c31=1.-c30-c32
!         c00=1./3.
!         c02=-1./6.
!         c01=1.-c00-c02
!         c10=-1./6.
!         c12=1./3.
!         c11=1.-c10-c12
!         c20=1./3.
!         c22=11./6.
!         c21=1.-c20-c22
!         c30=11./6.
!         c32=1./3.
!         c31=1.-c30-c32
!        calcul des flux d'ordre 3 sur les 3 stencils       
             f11=0.5*(1.+sign(1.,v1))*(q1f2m*c20 +q1f1m*c21 +q1f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     
             f12=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     
             f13=0.5*(1.+sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1p*c22 +q1f2p*c21 +q1f3p*c20)     
!                                                                       
             f21=0.5*(1.+sign(1.,v1))*(q2f2m*c20 +q2f1m*c21 +q2f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     
             f22=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     
             f23=0.5*(1.+sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1p*c22 +q2f2p*c21 +q2f3p*c20)     
!                                                                       
             f31=0.5*(1.+sign(1.,v1))*(q3f2m*c20 +q3f1m*c21 +q3f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     
             f32=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     
             f33=0.5*(1.+sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1p*c22 +q3f2p*c21 +q3f3p*c20)     
!                                                                       
             f41=0.5*(1.+sign(1.,v4))*(q4f2m*c20 +q4f1m*c21 +q4f  *c22)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     
             f42=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     &
                  +0.5*(1.-sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     
             f43=0.5*(1.+sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1p*c22 +q4f2p*c21 +q4f3p*c20)     
!                                                                       
             f51=0.5*(1.+sign(1.,v5))*(q5f2m*c20 +q5f1m*c21 +q5f  *c22)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     
             f52=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     &
                  +0.5*(1.-sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     
             f53=0.5*(1.+sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1p*c22 +q5f2p*c21 +q5f3p*c20)     
!        calcul des senseurs beta (au carre)                            
             iexp=2 
!         iexp=1                                                        
             s11=(q1f3p-2.*q1f2p+q1f1p)**2 
             s12=(q1f2p-2.*q1f1p+q1f  )**2 
             s13=(q1f1p-2.*q1f  +q1f1m)**2 
             s14=(q1f  -2.*q1f1m+q1f2m)**2 
             t11=(q1f3p-4.*q1f2p+3.*q1f1p)**2 
             t12=(q1f2p-4.*q1f1p+3.*q1f  )**2 
             t13=(q1f2p-q1f  )**2 
             t14=(q1f1p-q1f1m)**2 
             t15=(3.*q1f1p-4.*q1f  +q1f1m)**2 
             t16=(3.*q1f  -4.*q1f1m+q1f2m)**2 
             beta11=(0.5*(1.+sign(1.,v1))*(c1*s14+c2*t16)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s13+c2*t15)+eps)**iexp        
             beta12=(0.5*(1.+sign(1.,v1))*(c1*s13+c2*t14)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s12+c2*t13)+eps)**iexp        
             beta13=(0.5*(1.+sign(1.,v1))*(c1*s12+c2*t12)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s11+c2*t11)+eps)**iexp        
!                                                                       
             s21=(q2f3p-2.*q2f2p+q2f1p)**2 
             s22=(q2f2p-2.*q2f1p+q2f  )**2 
             s23=(q2f1p-2.*q2f  +q2f1m)**2 
             s24=(q2f  -2.*q2f1m+q2f2m)**2 
             t21=(q2f3p-4.*q2f2p+3.*q2f1p)**2 
             t22=(q2f2p-4.*q2f1p+3.*q2f  )**2 
             t23=(q2f2p-q2f  )**2 
             t24=(q2f1p-q2f1m)**2 
             t25=(3.*q2f1p-4.*q2f  +q2f1m)**2 
             t26=(3.*q2f  -4.*q2f1m+q2f2m)**2 
             beta21=(0.5*(1.+sign(1.,v1))*(c1*s24+c2*t26)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s23+c2*t25)+eps)**iexp        
             beta22=(0.5*(1.+sign(1.,v1))*(c1*s23+c2*t24)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s22+c2*t23)+eps)**iexp        
             beta23=(0.5*(1.+sign(1.,v1))*(c1*s22+c2*t22)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s21+c2*t21)+eps)**iexp        
!                                                                       
             s31=(q3f3p-2.*q3f2p+q3f1p)**2 
             s32=(q3f2p-2.*q3f1p+q3f  )**2 
             s33=(q3f1p-2.*q3f  +q3f1m)**2 
             s34=(q3f  -2.*q3f1m+q3f2m)**2 
             t31=(q3f3p-4.*q3f2p+3.*q3f1p)**2 
             t32=(q3f2p-4.*q3f1p+3.*q3f  )**2 
             t33=(q3f2p-q3f  )**2 
             t34=(q3f1p-q3f1m)**2 
             t35=(3.*q3f1p-4.*q3f  +q3f1m)**2 
             t36=(3.*q3f  -4.*q3f1m+q3f2m)**2 
             beta31=(0.5*(1.+sign(1.,v1))*(c1*s34+c2*t36)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s33+c2*t35)+eps)**iexp        
             beta32=(0.5*(1.+sign(1.,v1))*(c1*s33+c2*t34)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s32+c2*t33)+eps)**iexp        
             beta33=(0.5*(1.+sign(1.,v1))*(c1*s32+c2*t32)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s31+c2*t31)+eps)**iexp        
!                                                                       
             s41=(q4f3p-2.*q4f2p+q4f1p)**2 
             s42=(q4f2p-2.*q4f1p+q4f  )**2 
             s43=(q4f1p-2.*q4f  +q4f1m)**2 
             s44=(q4f  -2.*q4f1m+q4f2m)**2 
             t41=(q4f3p-4.*q4f2p+3.*q4f1p)**2 
             t42=(q4f2p-4.*q4f1p+3.*q4f  )**2 
             t43=(q4f2p-q4f  )**2 
             t44=(q4f1p-q4f1m)**2 
             t45=(3.*q4f1p-4.*q4f  +q4f1m)**2 
             t46=(3.*q4f  -4.*q4f1m+q4f2m)**2 
             beta41=(0.5*(1.+sign(1.,v4))*(c1*s44+c2*t46)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s43+c2*t45)+eps)**iexp        
             beta42=(0.5*(1.+sign(1.,v4))*(c1*s43+c2*t44)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s42+c2*t43)+eps)**iexp        
             beta43=(0.5*(1.+sign(1.,v4))*(c1*s42+c2*t42)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s41+c2*t41)+eps)**iexp        
!                                                                       
             s51=(q5f3p-2.*q5f2p+q5f1p)**2 
             s52=(q5f2p-2.*q5f1p+q5f  )**2 
             s53=(q5f1p-2.*q5f  +q5f1m)**2 
             s54=(q5f  -2.*q5f1m+q5f2m)**2 
             t51=(q5f3p-4.*q5f2p+3.*q5f1p)**2 
             t52=(q5f2p-4.*q5f1p+3.*q5f  )**2 
             t53=(q5f2p-q5f  )**2 
             t54=(q5f1p-q5f1m)**2 
             t55=(3.*q5f1p-4.*q5f  +q5f1m)**2 
             t56=(3.*q5f  -4.*q5f1m+q5f2m)**2 
             beta51=(0.5*(1.+sign(1.,v5))*(c1*s54+c2*t56)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s53+c2*t55)+eps)**iexp        
             beta52=(0.5*(1.+sign(1.,v5))*(c1*s53+c2*t54)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s52+c2*t53)+eps)**iexp        
             beta53=(0.5*(1.+sign(1.,v5))*(c1*s52+c2*t52)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s51+c2*t51)+eps)**iexp     
!        coefficients gamma en maillage irregulier
             g1p=bfg*bfh/((bfc+bfg)*(bfc+bfh))
             g3p=bfd*bfc/((bfd+bfh)*(bfc+bfh))
             g2p=1.-g1p-g3p
             g1m=bfh*bfi/((bfd+bfh)*(bfd+bfi))
             g3m=bfd*bfe/((bfd+bfi)*(bfe+bfi))
             g2m=1.-g1m-g3m
!         g1p=0.1
!         g3p=0.3
!         g2p=1.-g1p-g3p
!         g1m=0.3
!         g3m=0.1
!         g2m=1.-g1m-g3m
!        calculs des poids wi    
             ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta11)
             ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta12)
             ww31=0.5*(1.+sign(1.,v1))*(g3p/beta13) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta13)
             sw=ww11+ww21+ww31
             w11=ww11/sw
             w21=ww21/sw
             w31=ww31/sw
!         ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
!              +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
!         ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
!              +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
!         ww31m=w31*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w31+w31**2)/(g3p**2+w31*(1.-2.*g3p)) &        
!              +w31*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w31+w31**2)/(g3m**2+w31*(1.-2.*g3m))
!         swm=ww11m+ww21m+ww31m 
!         w11=ww11m/swm 
!         w21=ww21m/swm 
!         w31=ww31m/swm
!
             ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta21)
             ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta22)
             ww32=0.5*(1.+sign(1.,v1))*(g3p/beta23) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta23)
             sw=ww12+ww22+ww32
             w12=ww12/sw
             w22=ww22/sw
             w32=ww32/sw
!         ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
!              +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
!         ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
!              +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
!         ww32m=w32*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w32+w32**2)/(g3p**2+w32*(1.-2.*g3p)) &        
!              +w32*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w32+w32**2)/(g3m**2+w32*(1.-2.*g3m))
!         swm=ww12m+ww22m+ww32m 
!         w12=ww12m/swm 
!         w22=ww22m/swm 
!         w32=ww32m/swm
!
             ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta31)
             ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta32)
             ww33=0.5*(1.+sign(1.,v1))*(g3p/beta33) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta33)
             sw=ww13+ww23+ww33
             w13=ww13/sw
             w23=ww23/sw
             w33=ww33/sw
!         ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
!              +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
!         ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &        
!              +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
!         ww33m=w33*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w33+w33**2)/(g3p**2+w33*(1.-2.*g3p)) &        
!              +w33*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w33+w33**2)/(g3m**2+w33*(1.-2.*g3m))
!         swm=ww13m+ww23m+ww33m 
!         w13=ww13m/swm 
!         w23=ww23m/swm 
!         w33=ww33m/swm
!
             ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                  +0.5*(1.-sign(1.,v4))*(g1m/beta41)
             ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                  +0.5*(1.-sign(1.,v4))*(g2m/beta42)
             ww34=0.5*(1.+sign(1.,v4))*(g3p/beta43) &
                  +0.5*(1.-sign(1.,v4))*(g3m/beta43)
             sw=ww14+ww24+ww34
             w14=ww14/sw
             w24=ww24/sw
             w34=ww34/sw
!         ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
!              +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
!         ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &    
!              +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
!         ww34m=w34*0.5*(1.+sign(1.,v4))*(g3p+g3p**2-3.*g3p*w34+w34**2)/(g3p**2+w34*(1.-2.*g3p)) &        
!              +w34*0.5*(1.-sign(1.,v4))*(g3m+g3m**2-3.*g3m*w34+w34**2)/(g3m**2+w34*(1.-2.*g3m))
!         swm=ww14m+ww24m+ww34m 
!         w14=ww14m/swm 
!         w24=ww24m/swm 
!         w34=ww34m/swm
!
             ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                  +0.5*(1.-sign(1.,v5))*(g1m/beta51)
             ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                  +0.5*(1.-sign(1.,v5))*(g2m/beta52)
             ww35=0.5*(1.+sign(1.,v5))*(g3p/beta53) &
                  +0.5*(1.-sign(1.,v5))*(g3m/beta53)
             sw=ww15+ww25+ww35
             w15=ww15/sw
             w25=ww25/sw
             w35=ww35/sw
!         ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
!              +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
!         ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &        
!              +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
!         ww35m=w35*0.5*(1.+sign(1.,v5))*(g3p+g3p**2-3.*g3p*w35+w35**2)/(g3p**2+w35*(1.-2.*g3p)) &        
!              +w35*0.5*(1.-sign(1.,v5))*(g3m+g3m**2-3.*g3m*w35+w35**2)/(g3m**2+w35*(1.-2.*g3m))
!         swm=ww15m+ww25m+ww35m
!         w15=ww15m/swm 
!         w25=ww25m/swm 
!         w35=ww35m/swm
!        calcul des flux convectifs projetes
             fc1=w11*f11+w21*f12+w31*f13
             fc2=w12*f21+w22*f22+w32*f23
             fc3=w13*f31+w23*f32+w33*f33
             fc4=w14*f41+w24*f42+w34*f43
             fc5=w15*f51+w25*f52+w35*f53
!        produit avec matrice P pour retour dans l'espace physique
             f1=fc1*p11+fc2*p12+fc3*p13+fc4*p14+fc5*p15
             f2=fc1*p21+fc2*p22+fc3*p23+fc4*p24+fc5*p25
             f3=fc1*p31+fc2*p32+fc3*p33+fc4*p34+fc5*p35
             f5=fc1*p51+fc2*p52+fc3*p53+fc4*p54+fc5*p55
!-------------------------------------------------------------------    
!        produit de Q avec les flux Euler aux points (i-2) a (i+3)      
             q1f2m=q11*v(n-2*ninc,3)+q12*fxy(m-2*ninc)+q13*fyy(m-2*ninc)    &
                  +q14*fyz(m-2*ninc)+q15*fey(m-2*ninc)    
             q1f1m=q11*v(n-ninc  ,3)+q12*fxy(m-ninc)  +q13*fyy(m-ninc)      &
                  +q14*fyz(m-ninc)  +q15*fey(m-ninc)      
             q1f  =q11*v(n       ,3)+q12*fxy(m)       +q13*fyy(m)           &
                  +q14*fyz(m)       +q15*fey(m)           
             q1f1p=q11*v(n+ninc  ,3)+q12*fxy(m+ninc)  +q13*fyy(m+ninc)      &
                  +q14*fyz(m+ninc)  +q15*fey(m+ninc)      
             q1f2p=q11*v(n+2*ninc,3)+q12*fxy(m+2*ninc)+q13*fyy(m+2*ninc)    &
                  +q14*fyz(m+2*ninc)+q15*fey(m+2*ninc)    
             q1f3p=q11*v(n+3*ninc,3)+q12*fxy(m+3*ninc)+q13*fyy(m+3*ninc)    &
                  +q14*fyz(m+3*ninc)+q15*fey(m+3*ninc)    
!                                                                       
             q2f2m=q21*v(n-2*ninc,3)+q22*fxy(m-2*ninc)+q23*fyy(m-2*ninc)    &
                  +q24*fyz(m-2*ninc)+q25*fey(m-2*ninc)    
             q2f1m=q21*v(n-ninc  ,3)+q22*fxy(m-ninc)  +q23*fyy(m-ninc)      &
                  +q24*fyz(m-ninc)  +q25*fey(m-ninc)      
             q2f  =q21*v(n       ,3)+q22*fxy(m)       +q23*fyy(m)           &
                  +q24*fyz(m)       +q25*fey(m)           
             q2f1p=q21*v(n+ninc  ,3)+q22*fxy(m+ninc)  +q23*fyy(m+ninc)      &
                  +q24*fyz(m+ninc)  +q25*fey(m+ninc)      
             q2f2p=q21*v(n+2*ninc,3)+q22*fxy(m+2*ninc)+q23*fyy(m+2*ninc)    &
                  +q24*fyz(m+2*ninc)+q25*fey(m+2*ninc)    
             q2f3p=q21*v(n+3*ninc,3)+q22*fxy(m+3*ninc)+q23*fyy(m+3*ninc)    &
                  +q24*fyz(m+3*ninc)+q25*fey(m+3*ninc)    
!                                                                       
             q3f2m=q31*v(n-2*ninc,3)+q32*fxy(m-2*ninc)+q33*fyy(m-2*ninc)    &
                  +q34*fyz(m-2*ninc)+q35*fey(m-2*ninc)    
             q3f1m=q31*v(n-ninc  ,3)+q32*fxy(m-ninc)  +q33*fyy(m-ninc)      &
                  +q34*fyz(m-ninc)  +q35*fey(m-ninc)      
             q3f  =q31*v(n       ,3)+q32*fxy(m)       +q33*fyy(m)           &
                  +q34*fyz(m)       +q35*fey(m)           
             q3f1p=q31*v(n+ninc  ,3)+q32*fxy(m+ninc)  +q33*fyy(m+ninc)      &
                  +q34*fyz(m+ninc)  +q35*fey(m+ninc)      
             q3f2p=q31*v(n+2*ninc,3)+q32*fxy(m+2*ninc)+q33*fyy(m+2*ninc)    &
                  +q34*fyz(m+2*ninc)+q35*fey(m+2*ninc)    
             q3f3p=q31*v(n+3*ninc,3)+q32*fxy(m+3*ninc)+q33*fyy(m+3*ninc)    &
                  +q34*fyz(m+3*ninc)+q35*fey(m+3*ninc)    
!                                                                       
             q4f2m=q41*v(n-2*ninc,3)+q42*fxy(m-2*ninc)+q43*fyy(m-2*ninc)    &
                  +q44*fyz(m-2*ninc)+q45*fey(m-2*ninc)    
             q4f1m=q41*v(n-ninc  ,3)+q42*fxy(m-ninc)  +q43*fyy(m-ninc)      &
                  +q44*fyz(m-ninc)  +q45*fey(m-ninc)      
             q4f  =q41*v(n       ,3)+q42*fxy(m)       +q43*fyy(m)           &
                  +q44*fyz(m)       +q45*fey(m)           
             q4f1p=q41*v(n+ninc  ,3)+q42*fxy(m+ninc)  +q43*fyy(m+ninc)      &
                  +q44*fyz(m+ninc)  +q45*fey(m+ninc)      
             q4f2p=q41*v(n+2*ninc,3)+q42*fxy(m+2*ninc)+q43*fyy(m+2*ninc)    &
                  +q44*fyz(m+2*ninc)+q45*fey(m+2*ninc)    
             q4f3p=q41*v(n+3*ninc,3)+q42*fxy(m+3*ninc)+q43*fyy(m+3*ninc)    &
                  +q44*fyz(m+3*ninc)+q45*fey(m+3*ninc)    
!                                                                       
             q5f2m=q51*v(n-2*ninc,3)+q52*fxy(m-2*ninc)+q53*fyy(m-2*ninc)    &
                  +q54*fyz(m-2*ninc)+q55*fey(m-2*ninc)    
             q5f1m=q51*v(n-ninc  ,3)+q52*fxy(m-ninc)  +q53*fyy(m-ninc)      &
                  +q54*fyz(m-ninc)  +q55*fey(m-ninc)      
             q5f  =q51*v(n       ,3)+q52*fxy(m)       +q53*fyy(m)           &
                  +q54*fyz(m)       +q55*fey(m)           
             q5f1p=q51*v(n+ninc  ,3)+q52*fxy(m+ninc)  +q53*fyy(m+ninc)      &
                  +q54*fyz(m+ninc)  +q55*fey(m+ninc)      
             q5f2p=q51*v(n+2*ninc,3)+q52*fxy(m+2*ninc)+q53*fyy(m+2*ninc)    &
                  +q54*fyz(m+2*ninc)+q55*fey(m+2*ninc)    
             q5f3p=q51*v(n+3*ninc,3)+q52*fxy(m+3*ninc)+q53*fyy(m+3*ninc)    &
                  +q54*fyz(m+3*ninc)+q55*fey(m+3*ninc)    
!        calcul des flux d'ordre 3 sur les 3 stencils                   
             g11=0.5*(1.+sign(1.,v1))*(q1f2m*c20 +q1f1m*c21 +q1f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     
             g12=0.5*(1.+sign(1.,v1))*(q1f1m*c10 +q1f  *c11 +q1f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     
             g13=0.5*(1.+sign(1.,v1))*(q1f  *c20 +q1f1p*c11 +q1f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q1f1p*c22 +q1f2p*c21 +q1f3p*c20)     
!                                                                       
             g21=0.5*(1.+sign(1.,v1))*(q2f2m*c20 +q2f1m*c21 +q2f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     
             g22=0.5*(1.+sign(1.,v1))*(q2f1m*c10 +q2f  *c11 +q2f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     
             g23=0.5*(1.+sign(1.,v1))*(q2f  *c20 +q2f1p*c11 +q2f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q2f1p*c22 +q2f2p*c21 +q2f3p*c20)     
!                                                                       
             g31=0.5*(1.+sign(1.,v1))*(q3f2m*c20 +q3f1m*c21 +q3f  *c22)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     
             g32=0.5*(1.+sign(1.,v1))*(q3f1m*c10 +q3f  *c11 +q3f1p*c20)     &
                  +0.5*(1.-sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     
             g33=0.5*(1.+sign(1.,v1))*(q3f  *c20 +q3f1p*c11 +q3f2p*c10)     &
                  +0.5*(1.-sign(1.,v1))*(q3f1p*c22 +q3f2p*c21 +q3f3p*c20)     
!                                                                       
             g41=0.5*(1.+sign(1.,v4))*(q4f2m*c20 +q4f1m*c21 +q4f  *c22)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     
             g42=0.5*(1.+sign(1.,v4))*(q4f1m*c10 +q4f  *c11 +q4f1p*c20)     &
                  +0.5*(1.-sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     
             g43=0.5*(1.+sign(1.,v4))*(q4f  *c20 +q4f1p*c11 +q4f2p*c10)     &
                  +0.5*(1.-sign(1.,v4))*(q4f1p*c22 +q4f2p*c21 +q4f3p*c20)     
!                                                                       
             g51=0.5*(1.+sign(1.,v5))*(q5f2m*c20 +q5f1m*c21 +q5f  *c22)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     
             g52=0.5*(1.+sign(1.,v5))*(q5f1m*c10 +q5f  *c11 +q5f1p*c20)     &
                  +0.5*(1.-sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     
             g53=0.5*(1.+sign(1.,v5))*(q5f  *c20 +q5f1p*c11 +q5f2p*c10)     &
                  +0.5*(1.-sign(1.,v5))*(q5f1p*c22 +q5f2p*c21 +q5f3p*c20)     
!        calcul des senseurs beta (au carre)                            
             iexp=2 
!         iexp=1                                                        
             s11=(q1f3p-2.*q1f2p+q1f1p)**2 
             s12=(q1f2p-2.*q1f1p+q1f  )**2 
             s13=(q1f1p-2.*q1f  +q1f1m)**2 
             s14=(q1f  -2.*q1f1m+q1f2m)**2 
             t11=(q1f3p-4.*q1f2p+3.*q1f1p)**2 
             t12=(q1f2p-4.*q1f1p+3.*q1f  )**2 
             t13=(q1f2p-q1f  )**2 
             t14=(q1f1p-q1f1m)**2 
             t15=(3.*q1f1p-4.*q1f  +q1f1m)**2 
             t16=(3.*q1f  -4.*q1f1m+q1f2m)**2 
             beta11=(0.5*(1.+sign(1.,v1))*(c1*s14+c2*t16)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s13+c2*t15)+eps)**iexp        
             beta12=(0.5*(1.+sign(1.,v1))*(c1*s13+c2*t14)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s12+c2*t13)+eps)**iexp        
             beta13=(0.5*(1.+sign(1.,v1))*(c1*s12+c2*t12)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s11+c2*t11)+eps)**iexp        
!                                                                       
             s21=(q2f3p-2.*q2f2p+q2f1p)**2 
             s22=(q2f2p-2.*q2f1p+q2f  )**2 
             s23=(q2f1p-2.*q2f  +q2f1m)**2 
             s24=(q2f  -2.*q2f1m+q2f2m)**2 
             t21=(q2f3p-4.*q2f2p+3.*q2f1p)**2 
             t22=(q2f2p-4.*q2f1p+3.*q2f  )**2 
             t23=(q2f2p-q2f  )**2 
             t24=(q2f1p-q2f1m)**2 
             t25=(3.*q2f1p-4.*q2f  +q2f1m)**2 
             t26=(3.*q2f  -4.*q2f1m+q2f2m)**2 
             beta21=(0.5*(1.+sign(1.,v1))*(c1*s24+c2*t26)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s23+c2*t25)+eps)**iexp        
             beta22=(0.5*(1.+sign(1.,v1))*(c1*s23+c2*t24)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s22+c2*t23)+eps)**iexp        
             beta23=(0.5*(1.+sign(1.,v1))*(c1*s22+c2*t22)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s21+c2*t21)+eps)**iexp        
!                                                                       
             s31=(q3f3p-2.*q3f2p+q3f1p)**2 
             s32=(q3f2p-2.*q3f1p+q3f  )**2 
             s33=(q3f1p-2.*q3f  +q3f1m)**2 
             s34=(q3f  -2.*q3f1m+q3f2m)**2 
             t31=(q3f3p-4.*q3f2p+3.*q3f1p)**2 
             t32=(q3f2p-4.*q3f1p+3.*q3f  )**2 
             t33=(q3f2p-q3f  )**2 
             t34=(q3f1p-q3f1m)**2 
             t35=(3.*q3f1p-4.*q3f  +q3f1m)**2 
             t36=(3.*q3f  -4.*q3f1m+q3f2m)**2 
             beta31=(0.5*(1.+sign(1.,v1))*(c1*s34+c2*t36)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s33+c2*t35)+eps)**iexp        
             beta32=(0.5*(1.+sign(1.,v1))*(c1*s33+c2*t34)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s32+c2*t33)+eps)**iexp        
             beta33=(0.5*(1.+sign(1.,v1))*(c1*s32+c2*t32)                   &
                  +0.5*(1.-sign(1.,v1))*(c1*s31+c2*t31)+eps)**iexp        
!                                                                       
             s41=(q4f3p-2.*q4f2p+q4f1p)**2 
             s42=(q4f2p-2.*q4f1p+q4f  )**2 
             s43=(q4f1p-2.*q4f  +q4f1m)**2 
             s44=(q4f  -2.*q4f1m+q4f2m)**2 
             t41=(q4f3p-4.*q4f2p+3.*q4f1p)**2 
             t42=(q4f2p-4.*q4f1p+3.*q4f  )**2 
             t43=(q4f2p-q4f  )**2 
             t44=(q4f1p-q4f1m)**2 
             t45=(3.*q4f1p-4.*q4f  +q4f1m)**2 
             t46=(3.*q4f  -4.*q4f1m+q4f2m)**2 
             beta41=(0.5*(1.+sign(1.,v4))*(c1*s44+c2*t46)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s43+c2*t45)+eps)**iexp        
             beta42=(0.5*(1.+sign(1.,v4))*(c1*s43+c2*t44)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s42+c2*t43)+eps)**iexp        
             beta43=(0.5*(1.+sign(1.,v4))*(c1*s42+c2*t42)                   &
                  +0.5*(1.-sign(1.,v4))*(c1*s41+c2*t41)+eps)**iexp        
!                                                                       
             s51=(q5f3p-2.*q5f2p+q5f1p)**2 
             s52=(q5f2p-2.*q5f1p+q5f  )**2 
             s53=(q5f1p-2.*q5f  +q5f1m)**2 
             s54=(q5f  -2.*q5f1m+q5f2m)**2 
             t51=(q5f3p-4.*q5f2p+3.*q5f1p)**2 
             t52=(q5f2p-4.*q5f1p+3.*q5f  )**2 
             t53=(q5f2p-q5f  )**2 
             t54=(q5f1p-q5f1m)**2 
             t55=(3.*q5f1p-4.*q5f  +q5f1m)**2 
             t56=(3.*q5f  -4.*q5f1m+q5f2m)**2 
             beta51=(0.5*(1.+sign(1.,v5))*(c1*s54+c2*t56)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s53+c2*t55)+eps)**iexp        
             beta52=(0.5*(1.+sign(1.,v5))*(c1*s53+c2*t54)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s52+c2*t53)+eps)**iexp        
             beta53=(0.5*(1.+sign(1.,v5))*(c1*s52+c2*t52)                   &
                  +0.5*(1.-sign(1.,v5))*(c1*s51+c2*t51)+eps)**iexp   
!        calculs des poids wi    
             ww11=0.5*(1.+sign(1.,v1))*(g1p/beta11) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta11)
             ww21=0.5*(1.+sign(1.,v1))*(g2p/beta12) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta12)
             ww31=0.5*(1.+sign(1.,v1))*(g3p/beta13) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta13)
             sw=ww11+ww21+ww31
             w11=ww11/sw
             w21=ww21/sw
             w31=ww31/sw
!         ww11m=w11*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w11+w11**2)/(g1p**2+w11*(1.-2.*g1p)) &
!              +w11*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w11+w11**2)/(g1m**2+w11*(1.-2.*g1m))
!         ww21m=w21*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w21+w21**2)/(g2p**2+w21*(1.-2.*g2p)) &
!              +w21*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w21+w21**2)/(g2m**2+w21*(1.-2.*g2m))
!         ww31m=w31*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w31+w31**2)/(g3p**2+w31*(1.-2.*g3p)) &        
!              +w31*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w31+w31**2)/(g3m**2+w31*(1.-2.*g3m))
!         swm=ww11m+ww21m+ww31m 
!         w11=ww11m/swm 
!         w21=ww21m/swm 
!         w31=ww31m/swm
!
             ww12=0.5*(1.+sign(1.,v1))*(g1p/beta21) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta21)
             ww22=0.5*(1.+sign(1.,v1))*(g2p/beta22) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta22)
             ww32=0.5*(1.+sign(1.,v1))*(g3p/beta23) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta23)
             sw=ww12+ww22+ww32
             w12=ww12/sw
             w22=ww22/sw
             w32=ww32/sw
!         ww12m=w12*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w12+w12**2)/(g1p**2+w12*(1.-2.*g1p)) &
!              +w12*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w12+w12**2)/(g1m**2+w12*(1.-2.*g1m))
!         ww22m=w22*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w22+w22**2)/(g2p**2+w22*(1.-2.*g2p)) &
!              +w22*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w22+w22**2)/(g2m**2+w22*(1.-2.*g2m))
!         ww32m=w32*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w32+w32**2)/(g3p**2+w32*(1.-2.*g3p)) &        
!              +w32*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w32+w32**2)/(g3m**2+w32*(1.-2.*g3m))
!         swm=ww12m+ww22m+ww32m 
!         w12=ww12m/swm 
!         w22=ww22m/swm 
!         w32=ww32m/swm
!
             ww13=0.5*(1.+sign(1.,v1))*(g1p/beta31) &
                  +0.5*(1.-sign(1.,v1))*(g1m/beta31)
             ww23=0.5*(1.+sign(1.,v1))*(g2p/beta32) &
                  +0.5*(1.-sign(1.,v1))*(g2m/beta32)
             ww33=0.5*(1.+sign(1.,v1))*(g3p/beta33) &
                  +0.5*(1.-sign(1.,v1))*(g3m/beta33)
             sw=ww13+ww23+ww33
             w13=ww13/sw
             w23=ww23/sw
             w33=ww33/sw
!         ww13m=w13*0.5*(1.+sign(1.,v1))*(g1p+g1p**2-3.*g1p*w13+w13**2)/(g1p**2+w13*(1.-2.*g1p)) &
!              +w13*0.5*(1.-sign(1.,v1))*(g1m+g1m**2-3.*g1m*w13+w13**2)/(g1m**2+w13*(1.-2.*g1m))
!         ww23m=w23*0.5*(1.+sign(1.,v1))*(g2p+g2p**2-3.*g2p*w23+w23**2)/(g2p**2+w23*(1.-2.*g2p)) &        
!              +w23*0.5*(1.-sign(1.,v1))*(g2m+g2m**2-3.*g2m*w23+w23**2)/(g2m**2+w23*(1.-2.*g2m))
!         ww33m=w33*0.5*(1.+sign(1.,v1))*(g3p+g3p**2-3.*g3p*w33+w33**2)/(g3p**2+w33*(1.-2.*g3p)) &        
!              +w33*0.5*(1.-sign(1.,v1))*(g3m+g3m**2-3.*g3m*w33+w33**2)/(g3m**2+w33*(1.-2.*g3m))
!         swm=ww13m+ww23m+ww33m 
!         w13=ww13m/swm 
!         w23=ww23m/swm 
!         w33=ww33m/swm
!
             ww14=0.5*(1.+sign(1.,v4))*(g1p/beta41) &
                  +0.5*(1.-sign(1.,v4))*(g1m/beta41)
             ww24=0.5*(1.+sign(1.,v4))*(g2p/beta42) &
                  +0.5*(1.-sign(1.,v4))*(g2m/beta42)
             ww34=0.5*(1.+sign(1.,v4))*(g3p/beta43) &
                  +0.5*(1.-sign(1.,v4))*(g3m/beta43)
             sw=ww14+ww24+ww34
             w14=ww14/sw
             w24=ww24/sw
             w34=ww34/sw
!         ww14m=w14*0.5*(1.+sign(1.,v4))*(g1p+g1p**2-3.*g1p*w14+w14**2)/(g1p**2+w14*(1.-2.*g1p)) &
!              +w14*0.5*(1.-sign(1.,v4))*(g1m+g1m**2-3.*g1m*w14+w14**2)/(g1m**2+w14*(1.-2.*g1m))
!         ww24m=w24*0.5*(1.+sign(1.,v4))*(g2p+g2p**2-3.*g2p*w24+w24**2)/(g2p**2+w24*(1.-2.*g2p)) &    
!              +w24*0.5*(1.-sign(1.,v4))*(g2m+g2m**2-3.*g2m*w24+w24**2)/(g2m**2+w24*(1.-2.*g2m))
!         ww34m=w34*0.5*(1.+sign(1.,v4))*(g3p+g3p**2-3.*g3p*w34+w34**2)/(g3p**2+w34*(1.-2.*g3p)) &        
!              +w34*0.5*(1.-sign(1.,v4))*(g3m+g3m**2-3.*g3m*w34+w34**2)/(g3m**2+w34*(1.-2.*g3m))
!         swm=ww14m+ww24m+ww34m 
!         w14=ww14m/swm 
!         w24=ww24m/swm 
!         w34=ww34m/swm
!
             ww15=0.5*(1.+sign(1.,v5))*(g1p/beta51) &
                  +0.5*(1.-sign(1.,v5))*(g1m/beta51)
             ww25=0.5*(1.+sign(1.,v5))*(g2p/beta52) &
                  +0.5*(1.-sign(1.,v5))*(g2m/beta52)
             ww35=0.5*(1.+sign(1.,v5))*(g3p/beta53) &
                  +0.5*(1.-sign(1.,v5))*(g3m/beta53)
             sw=ww15+ww25+ww35
             w15=ww15/sw
             w25=ww25/sw
             w35=ww35/sw
!         ww15m=w15*0.5*(1.+sign(1.,v5))*(g1p+g1p**2-3.*g1p*w15+w15**2)/(g1p**2+w15*(1.-2.*g1p)) &
!              +w15*0.5*(1.-sign(1.,v5))*(g1m+g1m**2-3.*g1m*w15+w15**2)/(g1m**2+w15*(1.-2.*g1m))
!         ww25m=w25*0.5*(1.+sign(1.,v5))*(g2p+g2p**2-3.*g2p*w25+w25**2)/(g2p**2+w25*(1.-2.*g2p)) &        
!              +w25*0.5*(1.-sign(1.,v5))*(g2m+g2m**2-3.*g2m*w25+w25**2)/(g2m**2+w25*(1.-2.*g2m))
!         ww35m=w35*0.5*(1.+sign(1.,v5))*(g3p+g3p**2-3.*g3p*w35+w35**2)/(g3p**2+w35*(1.-2.*g3p)) &        
!              +w35*0.5*(1.-sign(1.,v5))*(g3m+g3m**2-3.*g3m*w35+w35**2)/(g3m**2+w35*(1.-2.*g3m))
!         swm=ww15m+ww25m+ww35m
!         w15=ww15m/swm 
!         w25=ww25m/swm 
!         w35=ww35m/swm
!        calcul des flux convectifs projetes
             gc1=w11*g11+w21*g12+w31*g13
             gc2=w12*g21+w22*g22+w32*g23
             gc3=w13*g31+w23*g32+w33*g33
             gc4=w14*g41+w24*g42+w34*g43
             gc5=w15*g51+w25*g52+w35*g53
!        produit avec matrice P pour retour dans l'espace physique
             g1=gc1*p11+gc2*p12+gc3*p13+gc4*p14+gc5*p15
             g2=gc1*p21+gc2*p22+gc3*p23+gc4*p24+gc5*p25
             g3=gc1*p31+gc2*p32+gc3*p33+gc4*p34+gc5*p35
             g5=gc1*p51+gc2*p52+gc3*p53+gc4*p54+gc5*p55
!        calcul du flux numerique et bilan de flux
             dg1=f1*sn(m1,kdir,1)+g1*sn(m1,kdir,2)
             dg2=f2*sn(m1,kdir,1)+g2*sn(m1,kdir,2)
             dg3=f3*sn(m1,kdir,1)+g3*sn(m1,kdir,2)
             dg5=f5*sn(m1,kdir,1)+g5*sn(m1,kdir,2)
!        calcul des flux visqueux (multiplies par -2)                   
             gv2=(cmuj2(m1)*toxx(n)+cmuj1(m1)*toxx(n1))*sn(m1,kdir,1) &
                  +(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,2)
             gv3=(cmuj2(m1)*toxy(n)+cmuj1(m1)*toxy(n1))*sn(m1,kdir,1) &
                  +(cmuj2(m1)*toyy(n)+cmuj1(m1)*toyy(n1))*sn(m1,kdir,2)
             gv5=(cmuj2(m1)*(toxx(n )*ul+toxy(n )*vl+qcx(n )) &
                  +cmuj1(m1)*(toxx(n1)*ur+toxy(n1)*vr+qcx(n1)))*sn(m1,kdir,1) &
                  +(cmuj2(m1)*(toxy(n )*ul+toyy(n )*vl+qcy(n )) &
                  +cmuj1(m1)*(toxy(n1)*ur+toyy(n1)*vr+qcy(n1)))*sn(m1,kdir,2)
             u(n1,1)=u(n1,1)-dg1 
             u(n1,2)=u(n1,2)-dg2+0.5*gv2 
             u(n1,3)=u(n1,3)-dg3+0.5*gv3 
             u(n1,5)=u(n1,5)-dg5+0.5*gv5 
             u(n,1)=u(n,1)+dg1 
             u(n,2)=u(n,2)+dg2-0.5*gv2 
             u(n,3)=u(n,3)+dg3-0.5*gv3 
             u(n,5)=u(n,5)+dg5-0.5*gv5 
          enddo
       enddo
    enddo
!                                                                       
!-----traitement des bords------------------------------------------    
!                                                                       
    kdir=1 
    ninc=nci 
!                                                                       
    do k=k1,k2m1
       ind1 = indc(i1,j1  ,k)
       ind2 = indc(i1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
!       flux a la facette frontiere
          f1=v(n-ninc,2)*sn(m,kdir,1) &
               + v(n-ninc,3)*sn(m,kdir,2)
          f2=(fxx(m-ninc)-toxx(n-ninc))*sn(m,kdir,1) &
               +(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,2)
          f3=(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,1) &
               +(fyy(m-ninc)-toyy(n-ninc))*sn(m,kdir,2)
          f5=(fex(m-ninc)-(toxx(n-ninc)*v(n-ninc,2)+toxy(n-ninc)* &
               v(n-ninc,3)+ toxz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcx(n-ninc))*sn(m,kdir,1) &
               +(fey(m-ninc)-(toxy(n-ninc)*v(n-ninc,2)+toyy(n-ninc)* &
               v(n-ninc,3)+ toyz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcy(n-ninc))*sn(m,kdir,2)
          u(n,1)=u(n,1)-f1
          u(n,2)=u(n,2)-f2
          u(n,3)=u(n,3)-f3
          u(n,5)=u(n,5)-f5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i2m1,j1  ,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2,ncj
          m=n-n0c
          m1=m+ninc
          n1=n+ninc
!       flux a la facette frontiere
          f1=v(n1,2)*sn(m1,kdir,1) &
               + v(n1,3)*sn(m1,kdir,2)
          f2=(fxx(m1)-toxx(n1))*sn(m1,kdir,1) &
               +(fxy(m1)-toxy(n1))*sn(m1,kdir,2)
          f3=(fxy(m1)-toxy(n1))*sn(m1,kdir,1) &
               +(fyy(m1)-toyy(n1))*sn(m1,kdir,2)
          f5=(fex(m1)-(toxx(n1)*v(n1,2)+toxy(n1)*v(n1,3) &
               + toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1))*sn(m1,kdir,1) &
               +(fey(m1)-(toxy(n1)*v(n1,2)+toyy(n1)*v(n1,3) &
               + toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1))*sn(m1,kdir,2)
          u(n,1)=u(n,1)+f1
          u(n,2)=u(n,2)+f2
          u(n,3)=u(n,3)+f3
          u(n,5)=u(n,5)+f5
       enddo
    enddo
!
    kdir=2
    ninc=ncj
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j1,k)
       ind2 = indc(i2m1,j1,k)
       do n=ind1,ind2
          m=n-n0c
!       flux a la facette frontiere
          g1=v(n-ninc,2)*sn(m,kdir,1) &
               + v(n-ninc,3)*sn(m,kdir,2)
          g2=(fxx(m-ninc)-toxx(n-ninc))*sn(m,kdir,1) &
               +(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,2)
          g3=(fxy(m-ninc)-toxy(n-ninc))*sn(m,kdir,1) &
               +(fyy(m-ninc)-toyy(n-ninc))*sn(m,kdir,2)
          g5=(fex(m-ninc)-(toxx(n-ninc)*v(n-ninc,2)+toxy(n-ninc)* &
               v(n-ninc,3)+ toxz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcx(n-ninc))*sn(m,kdir,1) &
               +(fey(m-ninc)-(toxy(n-ninc)*v(n-ninc,2)+toyy(n-ninc)* &
               v(n-ninc,3)+ toyz(n-ninc)*v(n-ninc,4))/v(n-ninc,1) &
               -qcy(n-ninc))*sn(m,kdir,2)
          u(n,1)=u(n,1)-g1
          u(n,2)=u(n,2)-g2
          u(n,3)=u(n,3)-g3
          u(n,5)=u(n,5)-g5
       enddo
    enddo
!
    do k=k1,k2m1
       ind1 = indc(i1  ,j2m1,k)
       ind2 = indc(i2m1,j2m1,k)
       do n=ind1,ind2
          m=n-n0c
          m1=m+ninc
          n1=n+ninc
!       flux a la facette frontiere
          g1=v(n1,2)*sn(m1,kdir,1) &
               + v(n1,3)*sn(m1,kdir,2)
          g2=(fxx(m1)-toxx(n1))*sn(m1,kdir,1) &
               +(fxy(m1)-toxy(n1))*sn(m1,kdir,2)
          g3=(fxy(m1)-toxy(n1))*sn(m1,kdir,1) &
               +(fyy(m1)-toyy(n1))*sn(m1,kdir,2)
          g5=(fex(m1)-(toxx(n1)*v(n1,2)+toxy(n1)*v(n1,3) &
               + toxz(n1)*v(n1,4))/v(n1,1)-qcx(n1))*sn(m1,kdir,1) &
               +(fey(m1)-(toxy(n1)*v(n1,2)+toyy(n1)*v(n1,3) &
               + toyz(n1)*v(n1,4))/v(n1,1)-qcy(n1))*sn(m1,kdir,2)
          u(n,1)=u(n,1)+g1
          u(n,2)=u(n,2)+g2
          u(n,3)=u(n,3)+g3
          u(n,5)=u(n,5)+g5
       enddo
    enddo
!
    if(isortie.eq.1) then
       write(6,'("===>sch_weno5_pond: increment explicite")')
       k=1
       i=13
!       i=160 
       do j=j1,j2m1
          n=indc(i,j,k)
          m=n-n0c
          write(6,'(i4,i6,4(1pe12.4))') &
               j,n,u(n,1),u(n,2),u(n,3),u(n,5)
       enddo
    endif
!
!-----calcul de la forcing function pour le multigrille--------------
!
    if(ityprk.ne.0) then
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                ff(n,1) = ff(n,1) - u(n,1) 
                ff(n,2) = ff(n,2) - u(n,2)
                ff(n,3) = ff(n,3) - u(n,3)
                ff(n,4) = ff(n,4) - u(n,4)
                ff(n,5) = ff(n,5) - u(n,5)
             enddo
          enddo
       enddo
    endif

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine sch_weno5pond
end module mod_sch_weno5pond
