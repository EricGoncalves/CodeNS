module para_fige
implicit none
      integer :: ndir,nind,lz,lg,lt,mtb,mtt,neqt,nsta,lsta,ista,nobj
      integer :: nmx,lgcmdx
  parameter(  ndir=3     )
  parameter(  nind=3     )
  parameter(    lz=50    )
  parameter(    lg=6     )
  parameter(    lt=lz*lg )
  parameter(   mtb=600   )
  parameter(   mtt=mtb*lg)
  parameter(  neqt=7     )
  parameter(  nsta=50    )
  parameter(  lsta=7     )
  parameter(  ista=2     )
  parameter(  nobj=600   )
  parameter(   nmx=500   )
  parameter(lgcmdx=1316  )
end module para_fige
!
module para_var
implicit none
      integer :: ndimub,ndimctf,ndimnts,ndimntu,kdimg,kdimv,kdimk
      integer :: mdimub,mdimtbf,mdimtnf,mdimtcf,mdimtrf,nvar
      integer :: ip00,ip11,ip12,ip13,ip60,ip40,ip21,ip31,ip41
      integer :: ip42,ip43,ip44
      double precision :: ccg,cfg,cng
  parameter(ndimub =120000)
  parameter(ndimctf=100000)
  parameter(ndimnts=100000)
  parameter(ndimntu=1)
  parameter(kdimg  =1)
  parameter(kdimv  =1)
  parameter(kdimk  =1)
  parameter(mdimub =400)
  parameter(mdimtbf=900)
  parameter(mdimtnf=900)
  parameter(mdimtcf=160)
  parameter(mdimtrf=1)
  parameter(nvar   =7)
!parameter(    ccg=1./3.) !3D
!parameter(    cng=1./3.) !3D
  parameter(    ccg=1./2.) 
  parameter(    cng=1./2.)
  parameter(    cfg=1.)
  parameter(ip00=ndimub)
  parameter(ip11=int((1.+kdimg*ccg)*ndimctf))
  parameter(ip12=kdimv*(ip11-1)+1)
  parameter(ip13=kdimk*(ip11-1)+1)
  parameter(ip21=int((1.+kdimg*cng)*ndimnts+ndimntu))
  parameter(ip31=int(1.+3.*(1.+kdimg*cng)*ndimnts))
  parameter(ip40=mdimub)
  parameter(ip41=int((1.+kdimg*cfg)*mdimtbf))
  parameter(ip42=int((1.+kdimg*cfg)*mdimtnf))
  parameter(ip43=int((1.+kdimg*cfg)*mdimtcf))
  parameter(ip44=int((1.+kdimg*cfg)*mdimtrf))
  parameter(ip60=nvar)
end module para_var
!
module boundary
  use para_fige
implicit none
  character(len=4) :: cl
  character(len=2) :: indfl
  integer iminb,imaxb,jminb,jmaxb,kminb,kmaxb
  integer mpn,nfbn
  integer mpr,nfbr,ndrr,srotr,crotr
  integer nbd,lbd,nbdko,lbdko
  integer nfba,nbdc,kexl
  integer mmb,mpb,nba,ndlb,nfei
  integer mpc,nfbc,ndcc,mdnc,mper
  real bc
  dimension nbdc(mtb),bc(mtb,ista*lsta)
  dimension iminb(mtt),imaxb(mtt),jminb(mtt),jmaxb(mtt),kminb(mtt),kmaxb(mtt)
  dimension mpn(mtt),nfbn(mtb)
  dimension mpr(mtt),nfbr(mtb),ndrr(mtb),srotr(mtb),crotr(mtb)
  dimension lbd(mtt),nfba(mtb),lbdko(mtt)
  dimension cl(mtb),indfl(mtb)
  dimension mmb(mtt),mpb(mtt),nba(mtb),ndlb(mtb),nfei(mtb)
  dimension mpc(mtt),nfbc(mtb),ndcc(mtb),mdnc(mtt),mper(mtt)
end module boundary
!
module maillage
  use para_fige
implicit none
  integer nptot,npn,nnn,npc,nnc,npfb,nnfb
  integer kvn,kcaldis,kecrdis,klecdis
  integer lzx,lgx,mtbx,mtnx,mtcx,mtrx,mtax,ndimubx,ndimctbx,ndimntbx, &
          mdimubx,mdimtbx,mdimtnx,mdimtcx,mdimtrx,neqtx
  integer id1,ii1,ii2,id2,jd1,jj1,jj2,jd2,kd1,kk1,kk2,kd2
  integer nbdrat,npbrat,lbdrat
  dimension npn(lt),nnn(lt),npc(lt),nnc(lt),npfb(lt),nnfb(lt)
  dimension id1(lt),ii1(lt),ii2(lt),id2(lt),jd1(lt),jj1(lt), &
            jj2(lt),jd2(lt),kd1(lt),kk1(lt),kk2(lt),kd2(lt)
  dimension nbdrat(lz),npbrat(lz),lbdrat(mtb)
end module maillage
!
module definition
  use para_fige
implicit none
  integer klomg
  real roa1,aa1,ta1,pa1,ha1
  real perio,ptrans,protat,omg
  real ronz,anz,tnz,dnz,pnz,rnz
  real varst
  dimension varst(nsta,lsta)
end module definition
!
module chainecarac
implicit none
  character(len=24) :: c0,c1,c2,c3
  character(len=32) :: cb,cc,cd,ch,ci,cf,cm,cr,cs
  data c0/'VALEUR NON SIGNIFICATIVE'/
  data c1/'valeur par defaut       '/
  data c2/'                        '/
  data c3/'valeur anterieure       '/
  data cb/'                                '/
  data cc/' character trop long!           '/
  data cd/' pas de domaine deja cree!      '/
  data ch/' character attendu!             '/
  data ci/' entier attendu!                '/
  data cf/' pas de frontiere deja creee!   '/
  data cm/' liste d''entiers attendue!     '/
  data cr/' reel attendu!                  '/
  data cs/' suite de commande attendue!    '/
  character(len=7)  :: equat,equatt
  character(len=80) :: titrt1
  character(len=4)  :: config
end module chainecarac
!
module constantes
implicit none
  integer intmx,linx  
  data intmx/999999/
  data linx/132/
  real reelmx,reelmn
  data reelmx/999999999./
  data reelmn/1.e-30/
  real pis2,raddeg,degrad
  data pis2/1.570796327/
  data raddeg/57.29577951/
  data degrad/0.01745329252/
end module constantes
!
module proprieteflu
implicit none
   real gam,gam1,gam2,gam3,gam4,gam5,rd
   real pr,prt,reynz
   real rgp,cp,cv
   real pinfl,ql
end module proprieteflu
!
module kcle
  use para_fige
implicit none
  integer klzx,klgx,kmtbx,kmtnx,kmtcx,kmtrx,kmtax,kndimubx,kndimctbx, &
          kndimntbx,kmdimubx,kmdimtbx, &
          kmdimtnx,kmdimtcx,kmdimtrx,kneqtx
  integer kkdualns,ktol,ktolke,kniter,knitur
  integer kischema,kmuscl,kilim,kxk
  integer kkprec,kcte,kkvisq
  integer kgam,krd,kpr,kprt,kreynz
  integer kpinfl,kql
  integer knba,kkvn,knfi,kkfmg,kkcg
  integer kncyresi,kncysave,kncyexpl,kdiscsv
  integer kequat,kconfig,kkexl,ktitrt1
  integer kklomg,komg,kperio
  integer kronz,kanz,ktnz,kdnz,kpnz,krnz
  integer kkdtl,kicychr0,kncychro,kdt1min,keta
  integer kki2,kki4,kkmf,klmax,knumt,kncycle
  integer kicytur0,kncyturb,kpctvort
  integer kvarst
  dimension keta(lt)
  dimension kki2(lt),kki4(lt)
  dimension kkmf(lt),klmax(lt)
  dimension kpctvort(lt)
  dimension kvarst(nsta,lsta)
  dimension kncycle(lg)
end module kcle
!
module schemanum
  use para_fige
implicit none
   character(len=4) :: discsv
  integer muscl,ilim,ischema,lacou
  integer kdualns,kdualto
  integer niter,nitur,mgl
  integer ktrans,kprec,kvisq,klroe
  integer ncyresi,ncysave,ncyexpl
  integer nfi,kfmg,kcg
  integer kmf,lmax,numt,ncycle
  integer kdtl,icychr0,ncychro
  real rki2t,rki4t,xk,cte,epsroe
  real resno1,resite,reske1,reskeite
  real rm0,al0,be0,tol,tolke
  real vrtmac,vrtalp,vrtcz,vrtlre,vrtxre,vrtzre
  real x0,y0,z0,freq,cga
  real ki2,ki4,eta,dt1min
  dimension ki2(lt),ki4(lt)
  dimension kmf(lt),lmax(lt)
  dimension ncycle(lg) 
  dimension eta(lt)
end module schemanum
!
module modeleturb
  use para_fige
implicit none
  integer keasm,keinit,kesst,kfracom,kditur,kwsst,ktransi
  integer ncycrac,naprng,ncytuke0,ncycke
  integer imxclko,komsst,icytur0,ncyturb,lparoi
  integer kcutke,kfludis,ksecmb,kcmut,kclkep,kinke,kparoi,kutau
  real utaumin,rapvisq
  real cb1,sigma,cb2,kappa,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4
  real rokinf,roeinf,epsk,epse,rkplus
  real cmu,cke1,cke2,alfak,alfae
  real rtrac,drtrac,vkar,cllog,yp0
  real sigme1,sigma1,beta1,wsig1,betae,okappa
  real sigme2,sigma2,beta2,wsig2
  real allfae0,allfa0,rrk,romeg,rbeta,bheta,sigmeb,sigmab,bethae
  real sigmk,sigmw,sigmd,beta,betas
  real ccmu,cc1,cc2,ce2,cgl,ceta,sigk,sige
  real cmukl,cklb1,ckle2,sigmak,sigmal,xkappa
  real epspid,epstaud,epsvord
  real pctvort
  dimension pctvort(lt)
end module modeleturb
!
module sortiefichier
  use para_fige
implicit none
  integer lsortie,nfreq
  integer kimp,kfa,lec,imp,out,sec,sor1,sor2,sor3
  integer kdgv,kdgc,kdac,kdgcf,kdacf,kres
  integer don1,inia1,sorf1,sorf2
  integer inig1,kdav,kfi,kfb,kfn,kfc,kfr
  data kimp/3/
  data lec,imp,out,sec,sor1,sor2,kfa,sor3/11,12,13,14,15,16,17,18/
  data kdgv,kdav,kdgc,kdac/21,22,23,24/
!  data inig1,kfi,kfb,kfn,kfc,kfr,kdgcf,kdacf,kres/31,32,33,34,35,36,37,38,39/
  data don1,inia1,sorf1,sorf2/41,42,43,44/
  integer kvglo,nbfll,nmfint
  real xref,yref,zref,sref,xlref
  real alpha0,beta0,p0spi0,q0spi0,v0
  dimension nmfint(mtb)
end module sortiefichier
!
!
program solve
!
!***********************************************************************
!
!     ACT
!_A   Logiciel de calcul par domaines d'ecoulements tridimensionnels
!_A   de fluide parfait par resolution des equations d'Euler
!_A   ou de fluide visqueux par resolution des equations de
!_A   resolution des equations de Navier-Stokes moyennees 
!_A   completees par un modele de turbulence.
!
!     VAL
!_V    gaz parfait et gaz raide.
!
!     INP
!_I    out        : com int          ; unite logiq, moyennes des residus
!_I    sec        : com int          ; unite logiq, pts a p ou ro negatifs
!_I    imp        : com int          ; unite logiq, sorties de controle
!_I    kdgv       : com int          ; unite logiq, coordonnees des noeuds
!_I    kdav       : com int          ; unite logiq, var aero aux noeuds
!_I    kdgc       : com int          ; unite logiq, coordonnees des centres
!_I    kdac       : com int          ; unite logiq, var aero aux centres
!_I    kdgcf      : com int          ; unite logiq, coordonnees des centres
!_I                                     et centres des facettes front
!_I    kdacf      : com int          ; unite logiq, var aero aux centres
!_I                                     et centres des facettes front
!_I    kfi        : com int          ; unite logiq, indices pts frontieres
!_I    kfb        : com int          ; unite logiq, tableaux de base front
!_I    kfn        : com int          ; unite logiq, tableaux normales
!_I    kfc        : com int          ; unite logiq, tableaux front coinc
!_I    kfr        : com int          ; unite logiq, tableaux front rec
!_I    kres       : com int          ; unite logiq, residus en chaque pt
!_I    reelmx     : com real         ; nombre reel grand
!
!     I/O
!
!     LOC
!_L    ki2        : com real(lt    ) ; coef de dissipation ordre 2    
!_L    ki4        : com real(lt    ) ; coef de dissipation ordre 4           
!_L    titrt1     : com char         ; titre du calcul                     
!_L    cb         : com char         ; mess interpretation, err indeterminee
!_L    cc         : com char         ; mess interpretation, char trop long   
!_L    cd         : com char         ; mess interpretation, pas de dom cree  
!_L    ch         : com char         ; mess interpretation, donnee non char  
!_L    ci         : com char         ; mess interpretation, donnee non entier
!_L    cf         : com char         ; mess interpretation, pas de fr creee  
!_L    cm         : com char         ; mess interpretation, donnee non liste d'entiers
!_L    cr         : com char         ; mess interpretation, donnee non reel 
!_L    cs         : com char         ; mess interpretation, suite de commande manquante  
!_L    equat      : com char         ; type d'equations modelisant l'ecoulement
!_L    cl         : com char(mtb   ) ; type de cond lim a appliquer
!_L    config     : com char         ; type de config geometrique du calcul
!_L    discsv     : com char         ; changement de discretisation (centre/noeud) 
!_L    td         : com char(lz    ) ; type de domaine (struct/non struct))
!_L    indfl      : com char(mtb   ) ; type de plan de la frontiere 
!_L    pis2       : com real         ; pi divise par 2 
!_L    raddeg     : com real         ; coef de transf rad en deg   
!_L    degrad     : com real         ; coef de transf deg en rad    
!_L    lzx        : com int          ; nbr total de domaines
!_L    lgx        : com int          ; nombre de niveaux de grille utilises
!_L    mtbx       : com int          ; nbr total de frontieres
!_L    mtnx       : com int          ; nbr total de frontieres a normales stockes
!_L    mtcx       : com int          ; nbr total de frontieres coincidentes
!_L    mtrx       : com int          ; nbr total de frontieres recouvertes
!_L    mtax       : com int          ; nbr total de frontieres autres     
!_L    ndimubx    : com int          ; nbr de cellules du plus grd domaine (pts fictifs inclus)  
!_L    ndimctbx   : com int          ; nbr de cellules de tts les domaines (pts fictifs inclus)
!_L    ndimntbx   : com int          ; nbr de noeuds de tts les domaines (pts fictifs inclus) 
!_L    mdimubx    : com int          ; nbr de pts de la plus grde front    
!_L    mdimtbx    : com int          ; nbr de pts de ttes les front        
!_L    mdimtnx    : com int          ; nbr de pts de ttes les front a normales stockees 
!_L    mdimtcx    : com int          ; nbr de pts de ttes les front coincidentes   
!_L    mdimtrx    : com int          ; nbr de pts de ttes les front recouvertes   
!_L    kimp       : com int          ; niveau de sortie sur unite logi imp
!_L    perio      : com real         ; periodicite geometrique en angle ou distance selon config
!_L    ptrans     : com real         ; distance pour periodicite    
!_L    protat     : com real         ; angle(rad) pour periodicite     
!_L    klomg      : com int          ; cle pour rotation du repere relatif  
!_L    omg        : com real         ; vitesse rotation du repere relatif
!_L    numt       : com int          ; cycle courant du calcul           
!_L    ncycle     : com int          ; nbr tot de cycles de l'execution courante 
!_L    npn        : com int (lt    ) ; pointeur fin de dom precedent dans tab tous noeuds
!_L    nnn        : com int (lt    ) ; nombre de noeuds du dom (dont fic.)  
!_L    npc        : com int (lt    ) ; pointeur fin de dom precedent dans tab toutes cellules
!_L    nnc        : com int (lt    ) ; nombre de cellules du dom (dont fic.) 
!_L    npfb       : com int (lt    ) ; pointeur fin de dom precedent dans tab toutes facettes
!_L    nnfb       : com int (lt    ) ; nombre de facettes du dom (dont fic.) 
!_L    id1        : com int (lt    ) ; indice min en i fictif         
!_L    ii1        : com int (lt    ) ; indice min en i reel 
!_L    ii2        : com int (lt    ) ; indice max en i reel
!_L    id2        : com int (lt    ) ; indice max en i fictif
!_L    jd1        : com int (lt    ) ; indice min en j fictif
!_L    jj1        : com int (lt    ) ; indice min en j reel
!_L    jj2        : com int (lt    ) ; indice max en j reel
!_L    jd2        : com int (lt    ) ; indice max en j fictif
!_L    kd1        : com int (lt    ) ; indice min en k fictif
!_L    kk1        : com int (lt    ) ; indice min en k reel
!_L    kk2        : com int (lt    ) ; indice max en k reel
!_L    kd2        : com int (lt    ) ; indice max en k fictif
!_L    npfq       : com int (lt    ) ; pointeur fin de dom precedent dans tab toutes facettes quad
!_L    nnfq       : com int (lt    ) ; nbr facettes quad du dom (dont fic.)  
!_L    nnfql      : com int (lt    ) ; nbr facettes quad limites du dom   
!_L    nnfqc      : com int (lt    ) ; nbr facettes quad coincidentes du dom 
!_L    npft       : com int (lt    ) ; pointeur fin de dom precedent dans tab toutes facettes triang 
!_L    nnft       : com int (lt    ) ; nbr facettes trian du dom (dont fic.) 
!_L    nnftl      : com int (lt    ) ; nbr facettes triang limites du dom   
!_L    nnftc      : com int (lt    ) ; nbr facettes triang coinc du dom      
!_L    nbd        : com int          ; nombre de frontieres a traiter       
!_L    lbd        : com int (mtt   ) ; numero de front a traiter             
!_L    mmb        : com int (mtt   ) ; nombre de pts d'une frontiere         
!_L    mpb        : com int (mtt   ) ; pointeur fin de front precedente dans tableaux de base des front.
!_L    nba        : com int (mtb   ) ; rang de traitement d'une front       
!_L    ndlb       : com int (mtb   ) ; numero dom contenant la frontiere    
!_L    nfei       : com int (mtb   ) ; numero de base interne d'une front en fct du numero externe 
!_L    iminb      : com int (mtt   ) ; indice min en i d'une front
!_L    imaxb      : com int (mtt   ) ; indice max en i d'une front
!_L    jminb      : com int (mtt   ) ; indice min en j d'une front
!_L    jmaxb      : com int (mtt   ) ; indice max en j d'une front
!_L    kminb      : com int (mtt   ) ; indice min en k d'une front
!_L    kmaxb      : com int (mtt   ) ; indice max en k d'une front
!_L    mpr        : com int (mtt   ) ; pointeur fin de front precedente dans tab front recouvertes 
!_L    nfbr       : com int (mtb   ) ; numero dans numerotation interne d'une frontiere recouverte
!_L    ndrr       : com int (mtb   ) ; numero de domaine recouvrant     
!_L    srotr      : com real(mtb   ) ; rotation amenant une front periodique dans sa position recouverte, sinus  
!_L    crotr      : com real(mtb   ) ; rotation amenant une front periodique dans sa position recouverte, cosinus 
!_L    mpc        : com int (mtt   ) ; pointeur fin de front precedente dans tableaux front coinc    
!_L    nfbc       : com int (mtb   ) ; numero dans numerotation interne d'une frontiere coincidente
!_L    ndcc       : com int (mtb   ) ; numero du dom coicident     
!_L    mdnc       : com int (mtt   ) ; saut d'ind entre pt front fictif et pt interieur pour la front coinc
!_L    mpn        : com int (mtt   ) ; pointeur fin de front precedente dans tab front a normales stockees
!_L    nfbn       : com int (mtb   ) ; numero dans numerotation interne
!
!_L    nfba       : com int (mtb   ) ; numero dans numerotation interne d'une frontiere autre   
!_L    gam        : com real         ; rapport des chaleurs specifiques
!_L    gam1       : com real         ; rap chal spec -1                   
!_L    gam2       : com real         ; (rap chal spec -1)/2                 
!_L    gam3       : com real         ; 1/rap chal spec                     
!_L    gam4       : com real         ; 1/(rap chal spec -1                  
!_L    gam5       : com real         ; rap chal spec/(rap chal spec -1)     
!_L    rgp        : com real         ; constante des gaz parfaits adim       
!_L    cp         : com real         ; chal spec a pres cste adim            
!_L    cv         : com real         ; chal spec a vol cst adim             
!_L    pr         : com real         ; nombre de Prandtl                   
!_L    prt        : com real         ; nombre de Prandtl turbulent           
!_L    reynz      : com real         ; nombre de Reynolds calcule avec les grandeurs d'adimensionnement,
!_L                                    pour definir la loi de Sutherland     
!_L    nfi        : com int          ; nbr de rangees de pts fictifs         
!_L    kdit       : com int          ; type dissipation artificielle        
!_L    kfmg       : com int          ; choix technique  multigrille         
!_L    kcg        : com int          ; cle calcul metrique grilles grossieres;
!_L    tnz        : com real         ; etat pour adimensionnement, temperature
!_L    ronz       : com real         ; etat pour adimensionnement, masse volumique
!_L    anz        : com real         ; etat pour adimensionnement, vitesse du son d'arret
!_L    dnz        : com real         ; etat pour adimensionnement longueur 
!_L    roa1       : com real         ; etat de reference utilisateur adimensionne, masse volumique d'arret
!_L    aa1        : com real         ; etat de reference utilisateur adimensionne, vitesse du son d'arret
!_L    ta1        : com real         ; etat de reference utilisateur adimensionne, temperature d'arret
!_L    pa1        : com real         ; pression d'arret de l'etat de reference utilisateur adimensionne
!_L    ha1        : com real         ; enthalpie d'arret de l'etat de reference utilisateur adimensionne
!_L    icytur0    : com int          ; nbr de cycl en deb de calcul au cours desquelles mut n'est pas mis a jour
!_L    ncyturb    : com int          ; freq en it de mise a jour de mut
!_L    pctvort    : com real(lt    ) ; pourcentage de tourbillon pour calcul d'epaisseur de couche lim 
!_L    kdtl       : com int          ; cle d'utilisation pas de temps local
!_L    icychr0    : com int          ; nbr de cycle en deb de calcul au cours desquelles le pas de temps 
!_L                                    est mis a jour a chaque it
!_L    ncychro    : com int          ; freq en it de mise a jour du pas de temps
!_L    dt1min     : com real         ; pas de temps constant desire
!_L    eta        : com real(lt    ) ; nombre de CFL
!_L    ncyresi    : com int          ; freq en it de calcul des residus
!_L    ncysave    : com int          ; freq en it de sauvegarde des var aero
!_L    kmf        : com int (lt    ) ; cle phase implicite
!_L    kvn        : com int          ; cle controle mailles deversees
!_L    kexl       : com int          ; cle traitement frontieres avant exploitation
!
!     COM
!_C
!_C   definitions :
!_C   -----------
!_C   sous-domaine         : partie du domaine de calcul maillee avec un maillage structure i,j,k .
!_C   frontiere            : ensemble de points appartenant a certaines des six surfaces limites des   
!_C                          sous-domaines, en tout point duquel on applique une meme condition a la limite .
!_C   frontiere coincidente: partie d' une des six surfaces limites d'un sous-domaine , qui coincide avec 
!_C                          une autre partie d'une des six surfaces limites d'un sous-domaine.
!_C   frontiere non coincidente: partie d' une des six surfaces limites d'un sous-domaine, qui est adjacente
!_C                              a une autre partie d'une des six surfaces limites d'un sous-domaine.
!_C   frontiere recouverte : frontiere limite d'un sous-domaine dont les
!_C                          points sont recouverts par un autre sous-domaine
!_C
!_C   dimensions :
!_C   ----------
!_C     les dimensions sont definies dans des instructions 'parameter' .
!_C     selon le cas de calcul il faut ajuster:
!_C
!_C     ndimub  : nbr max de pts d'un dom (pts fictifs inclus)
!_C     ndimctf : nbr max de cellules de tous les dom (pts fictifs inclus)
!_C     ndimnts : nbr max de noeuds de ts les dom struc (pts fictifs inclus)
!_C     ndimntu : nbr max de noeuds de ts les dom non-struc (pts fictifs inclus)
!_C     kdimg   : 1 pour multigrille si non 0
!_C     kdimv   : 1 pour Navier-Stokes si non 0
!_C     kdimk   : 1 pour k-epsilon si non 0
!_C     mdimub  : nbr max de facettes d'une front
!_C     mdimtbf : nbr max de facettes de ttes les front
!_C     mdimtnf : nbr max de facettes de ttes les front a normales stockees
!_C     mdimtcf : nbr max de facettes de ttes les front coincidentes
!_C     mdimtrf : nbr max de facettes de ttes les front recouvertes
!_C     neqt    : nbr max d'equations
!_C     ccg     : coef de calc du nbr de centres   :   1./3. (3d) ou 4./9. (2d)
!_C     cng     : coef de calc du nbr de noeuds    :   1./3. (3d) ou 4./9. (2d)
!_C     cfg     : coef de calc du nbr de pts       :   1.    (2d) ou 1.    (1d)
!_C
!_C     les dimensions suivantes sont figees dans le logiciel:
!_C
!_C     ndir    : nbr max de directions d'espace   :     3
!_C     nind    : nbr max d'indices en structure   :     3
!_C     lz      : nbr max de domaines              :    50
!_C     lg      : nbr max de niveaux de grille     :     6
!_C     lt      : nbr max de grilles               : lz*lg
!_C     mtb     : nbr max de frontieres            :   600
!_C     mtt     : nbr max de surfaces              :mtb*lg
!_C     nsta    : nbr max d'etats                  :    50
!_C     nobj    : nbr max d'objets                 :   600
!_C     nmx     : nbr max de mots                  :    50
!_C     lgcmdx  : longueur max d'une commande      :  1316
!_C
!_C   normalisation des variables :
!_C   ---------------------------
!_C     les coordonnees des points des maillages et les variables aerodynamiques
!_C     doivent etre definies par des nombres de l'ordre de l' unite.
!_C
!_C   fichiers :
!_C   --------
!_C     lec   - 11 - flec       -     formatte - donnees du calcul
!_C     imp   - 12 - fimp       -     formatte - sorties de controle
!_C     out   - 13 - fout       -     formatte - moyennes des residus et cx,xy,cz
!_C     sec   - 14 - fsec       -     formatte - pts a p ou ro negatifs
!_C
!_C     kdgv  - 21 - fgv        - non formatte - coordonnees des noeuds
!_C     kdav  - 22 - fav        - non formatte - variables aerodynamiques aux noeuds
!_C     kdgc  - 23 - fgc        - non formatte - coordonnees des centres de mailles
!_C     kdac  - 24 - fac        - non formatte - variables aerodynamiques
!_C                                              aux centres des mailles
!_C     kdgcf - 25 - fgcf       - non formatte - coordonnees des centres de
!_C                                              mailles et de facettes frontieres
!_C     kdacf - 26 - facf       - non formatte - variables aerodynamiques aux centres
!_C                                              des mailles etdes facettes frontieres
!_C     kfi   - 31 - fi         - non formatte - indice dans le domaine
!_C     kfb   - 32 - fb         - non formatte - increment vers l'interieur
!_C     kfn   - 33 - fn         - non formatte - normales frontieres
!_C     kfc   - 34 - fc         - non formatte - indice du point coincident
!_C     kfr   - 35 - fr         - non formatte - indice maille recouvrante
!_C                                              et coefficients d'interpolations.
!_C     kres  - 41 - fres       - non formatte - residus aux centres des mailles
!
!***********************************************************************
!
!-----parameters --------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use maillage
      use definition
      use chainecarac
      use sortiefichier
implicit none
!
!-----------------------------------------------------------------------
!
      integer           :: iyplus
      integer           :: img,imot,l,mfbi,mfc,mfn,mfr,mnc,mnpar,mnr,ncbd,ncin
      integer           :: ncyc,nmot
      real              :: mu,mut,nxn,nyn,nzn
      double precision  :: aam,bceqt,cfke,cmui1,cmui2,cmuj1,cmuj2,cmuk1
      double precision  :: cmuk2,cson,cvi,cvj,cvk,d0x,d0y,d0z,dist,dt,exs1
      double precision  :: exs2,fgam,pres,pression,ptdual,qcx,qcy,qcz,qtx,qty,qtz,r,res,roam
      double precision  :: rod,roed,roud,rovd,rowd,rpi,rti,sn,tam,tm1,tm10,tm11,tm12,tm13
      double precision  :: tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tn1,tn10,tn2,tn3,tn4,tn5,tn6,tn7
      double precision  :: tn8,tn9,tnte1,tnte2,tnte3,tnte4,toxx,toxy,toxz,toyy,toyz,tozz,tp
      double precision  :: utau,v,vdual,vdual1,vdual2,vol,x,xnr,y,ynr,z,znr,ztemp
      character(len=32) :: comment,mot(nmx)
!
      dimension imot(nmx)
      dimension v(ip11,ip60)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60),ptdual(ip11,ip60)
      dimension r(ip11),dt(ip11),vol(ip11),ztemp(ip11),pression(ip11),cson(ip11)
      dimension mu(ip12),mut(ip12),dist(ip12), &
                toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12), &
                tozz(ip12),qcx(ip12),qcy(ip12),qcz(ip12),mnpar(ip12)
      dimension cfke(ip13)
      dimension sn(ip31*ndir)
      dimension x(ip21),y(ip21),z(ip21)
      dimension cvi(ip21),cvj(ip21),cvk(ip21), &
                cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21),cmuk1(ip21),cmuk2(ip21)
      dimension rpi(ip40),rti(ip40),d0x(ip40),d0y(ip40),d0z(ip40), &
                qtx(ip40),qty(ip40),qtz(ip40),res(ip40),tp(ip40), &
                rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension ncbd(ip41),ncin(ip41)
      dimension bceqt(ip41,neqt)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),fgam(ip42),utau(ip42)
      dimension mnc(ip43)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
!     tableaux de travail
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40), &
                tm5(ip40),tm6(ip40),tm7(ip40),tm8(ip40), &
                tm9(ip40),tm10(ip40),tm11(ip40),tm12(ip40),tm13(ip40)
!    tm14(ip40),tm15(ip40),tm16(ip40), &
!                tm17(ip40),tm18(ip40),tm19(ip40),tm20(ip40), &
!                tm21(ip40),tm22(ip40),tm23(ip40),tm24(ip40), &
!                tm25(ip40),tm26(ip40),tm27(ip40),tm28(ip40), &
!                tm29(ip40),tm30(ip40),tm31(ip40),tm32(ip40), &
!                tm33(ip40),tm34(ip40),tm35(ip40),tm36(ip40), &
!                tm37(ip40),tm38(ip40),tm39(ip40),tm40(ip40), &
!                tm41(ip40),tm42(ip40),tm43(ip40),tm44(ip40)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8(ip00),tn9(ip00),tn10(ip00)
      dimension tnte1(ip11,ip60),tnte2(ip11,ip60),tnte3(ip11,ip60),tnte4(ip11,ip60)
!
      data exs1,exs2/1.,0./
!      data exr1,exr2/2.,-1./
!
!     reservation des unites logiques generales
!
      open(lec  ,file='flec')
      open(imp  ,file='fimp')
      open(out  ,file='fout')
      open(sec  ,file='fsec')
      open(sor1 ,file='smoy')
      open(sor2 ,file='pres')
      open(sor3 ,file='resro')
      open(kfa  ,file='fcla',form='formatted')
      open(kdgv ,file='fgv' ,form='unformatted')
      open(kdgc ,file='fgc' ,form='unformatted')
      open(kdac ,file='fac' ,form='unformatted')
!      open(kdav ,file='fav' ,form='unformatted')
!      open(kdgcf,file='fgcf',form='unformatted')
!      open(kdacf,file='facf',form='unformatted')
!      open(kfi  ,file='fi'  ,form='unformatted')
!      open(kfb  ,file='fb'  ,form='unformatted')
!      open(kfn  ,file='fn'  ,form='unformatted')
!      open(kfc  ,file='fc'  ,form='unformatted')
!      open(kfr  ,file='fr'  ,form='unformatted')
!      open(kres ,file='fres',form='unformatted')
!
!     lecture fichier "fatdon" des donnees des modeles 2 equations
      call atlecdon
!
!     initialisations
!
       call inimem(          &
                 dt,v,mu,mut, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 sn, &
                 vol, &
                 ptdual,vdual,vdual1,vdual2, &
                 cvi,cvj,cvk, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                 pression,ztemp,cson, &
                 tnte1,tnte3,tnte4, &
                 ncyc, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 comment)
!
!  lecture et interpretation des donnees  ******************************
!
      mot=""
      do while(mot(1)(1:3).ne.'end')
!
      call rdcmd(mot,imot,nmot)
!
!--   #
      if((imot(1).eq.1).and.(mot(1)(1:1).eq.'#')) then
        continue
!--   CALL
      elseif((imot(1).eq.4).and.(mot(1)(1:4).eq.'call')) then
!--   CALL UTDON
        if((imot(2).eq.9).and.(mot(2)(1:9).eq.'utdon_gen')) then
            call utdon_gen( &
                 config,cl,x,y,z,omg, &
                 ncbd,v, &
                 nxn,nyn,nzn)
!
!--   CALL UTSOR
        elseif((imot(2).eq.5).and.(mot(2)(1:5).eq.'utsor')) then
          if(equat(1:2).eq.'ns') then
            call utsorfr( &
                 ncbd,ncin,v,mu,mut, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 x,y,z,nxn,nyn,nzn, &
                 pression,ztemp,cson)
!
!           calcul et ecriture de y+ pour la premiere maille
!            hauteur demi-maille adjacente aux parois
            iyplus=1
            if(iyplus.eq.1) then
            call  met_yplus( &
                 ncbd,ncin,v,mu,dist, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 x,y,z,nxn,nyn,nzn)
            endif
!
!           integration des epaisseurs de couche limite
!
            call met_intep3( &
                 ncbd,ncin,v, &
                 sn,vol, &
                 dist,mnpar,mu, &
                 tn1,tn2,tn3,tn4,tn5, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 x,y,z,nxn,nyn,nzn, &
                 pression,cson,ztemp, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          else
!           calcul Euler
!
!            call utsor( &
!                 ncbd,v,mut, &
!                 x,y,z,nxn,nyn,nzn)
          endif
!
        else
            call synterr(mot,imot,2,cb)
        endif
!--   COMPUTE
      elseif  ((imot(1).eq.7).and.(mot(1)(1:7).eq.'compute')) then
!--   COMPUTE FLOW
        if((imot(2).eq.4).and.(mot(2)(1:4).eq.'flow')) then
!
            call c_cpfw( &
                 mot,imot,nmot, &
                 ncyc, &
                 x,y,z,r,exs1,exs2,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 mu,mut,dist,cfke, &
                 mnpar,fgam,utau, &
                 v,dt, &
                 ptdual,vdual,vdual1,vdual2, &
                 tnte1,tnte2,tnte3,tnte4, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
                 tm11,tm12,tm13, &
                 ncin, &
                 mnc, &
                 ncbd,mnr,xnr,ynr,znr, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,pres,tp, &
                 rod,roud,rovd,rowd,roed, &
                 pression,ztemp,cson, &
                 cvi,cvj,cvk, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!--   COMPUTE BOUNDARY
        else if((imot(2).eq.8).and.(mot(2)(1:8).eq.'boundary')) then
            call c_cpbd( &
                 mot,imot,nmot, &
                 ncin,nxn,nyn,nzn,ncbd, &
                 sn,vol,v,mut, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,x,y,z,omg, &
                 pres,tp,rod,roud,rovd,rowd,roed, &
                 mnr,xnr,ynr,znr,mnc, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13,pression,ztemp,cson)
!
        else
            call synterr(mot,imot,2,cb)
        endif
!--   CREATE
      elseif((imot(1).eq.6).and.(mot(1)(1:6).eq.'create')) then
!--   CREATE DOM
        if((imot(2).eq.3).and.(mot(2)(1:3).eq.'dom')) then
!--   CREATE DOM ST
          if((imot(3).eq.2).and.(mot(3)(1:2).eq.'st')) then
            call c_crdms(mot,imot,nmot)
          else
            call synterr(mot,imot,3,cb)
          endif
!--   CREATE BOUNDARY
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'boundary')) then
!--   CREATE BOUNDARY ST
          if((imot(3).eq.2).and.(mot(3)(1:2).eq.'st')) then
            call c_crbds(mot,imot,nmot,ncbd)
          else
            call synterr(mot,imot,3,cb)
          endif
        else
            call synterr(mot,imot,2,cb)
        endif
!--   DEFINE
      elseif((imot(1).eq.6).and.(mot(1)(1:6).eq.'define')) then
!--   DEFINE TITLE
        if((imot(2).eq.5).and.(mot(2)(1:5).eq.'title')) then
            call c_dftl1(mot,imot,nmot)
!--   DEFINE GEOMETRY
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'geometry')) then
            call c_dfgm(mot,imot,nmot)
!--   DEFINE FLOW
        elseif((imot(2).eq.4).and.(mot(2)(1:4).eq.'flow')) then
            call c_dffw(mot,imot,nmot)
!--   DEFINE PHYSICS
        elseif((imot(2).eq.7).and.(mot(2)(1:7).eq.'physics')) then
            call c_dfph(mot,imot,nmot)
!--   DEFINE NUMERICS
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'numerics')) then
            call c_dfnm(mot,imot,nmot)
!--   DEFINE STATE
        elseif((imot(2).eq.5).and.(mot(2)(1:5).eq.'state')) then
            call c_dfst(mot,imot,nmot)
!--   DEFINE NORMALIZ
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'normaliz')) then
            call c_dfnzst(mot,imot,nmot)
!
!     definition d'un etat amont
      roam=ronz
      aam=anz
      tam=tnz
!
      call c_dfst0(roam,aam,tam)
!
!--   normalisation etat
!
      call c_nzst(roam,aam,tam)
!
!--   DEFINE PM_DTD
        elseif((imot(2).eq.6).and.(mot(2)(1:6).eq.'pm_dtd')) then
            call c_dfpmdtd(mot,imot,nmot)
!--   DEFINE PM_DTG
        elseif((imot(2).eq.6).and.(mot(2)(1:6).eq.'pm_dtg')) then
            call c_dfpmdtg(mot,imot,nmot)
!--   DEFINE PM_TURBN
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'pm_turbn')) then
            call c_dfpmtbn(mot,imot,nmot)
!--   DEFINE PM_TURBP
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'pm_turbp')) then
!--   DEFINE PM_TURBP KEPS
          if((imot(3).eq.4).and.(mot(3)(1:4).eq.'keps')) then
!--   DEFINE PM_TURBP KEPS GENERAL
            if((imot(4).eq.7).and.(mot(4)(1:7).eq.'general')) then
            call c_dfpmtbkeg(mot,imot,nmot)
            else
                call synterr(mot,imot,4,cb)
            endif
          else
              call synterr(mot,imot,3,cb)
          endif
!--   DEFINE PM_NUMD
        elseif((imot(2).eq.7).and.(mot(2)(1:7).eq.'pm_numd')) then
            call c_dfpmdsd(mot,imot,nmot)
!--   DEFINE PM_NUMI
        elseif((imot(2).eq.7).and.(mot(2)(1:7).eq.'pm_numi')) then
            call c_dfpmimd(mot,imot,nmot)
!--   DEFINE PM_CFG
        elseif((imot(2).eq.6).and.(mot(2)(1:6).eq.'pm_cfg')) then
            call c_dfpmcfg(mot,imot,nmot)
        else
            call synterr(mot,imot,2,cb)
        endif
!--   DISPLAY
      elseif((imot(1).eq.7).and.(mot(1)(1:7).eq.'display')) then
!--   DISPLAY BOUNDARY
        if((imot(2).eq.8).and.(mot(2)(1:8).eq.'boundary')) then
            call c_dpbd( &
                 mot,imot,nmot, &
                 ncbd,ncin, &
                 nxn,nyn,nzn, &
                 mnc, &
                 mnr,xnr,ynr,znr, &
                 tm1,tm2,tm3, &
                 tm4,tm5,tm6)
!--   DISPLAY DIMENSION
        elseif((imot(2).eq.9).and.(mot(2)(1:9).eq.'dimension')) then
            call c_dpdim(mot,imot,nmot)
        endif
!--   END
      elseif((imot(1).eq.3).and.(mot(1)(1:3).eq.'end')) then
         call c_end(mot,imot,nmot)
!--   INIT
      elseif((imot(1).eq.4).and.(mot(1)(1:4).eq.'init')) then
!--   INIT DOM
        if((imot(2).eq.3).and.(mot(2)(1:3).eq.'dom')) then
!--   INIT DOM XYZ
          if((imot(3).eq.3).and.(mot(3)(1:3).eq.'xyz')) then
            call c_ingr( &
                 mot,imot,nmot, &
                 x,y,z)
!--   INIT DOM CSVMUT
          elseif((imot(3).eq.6).and.(mot(3)(1:6).eq.'csvmut')) then
            call c_infw( &
                 mot,imot,nmot, &
                 x,y,z,v,mut,tnte1,utau, &
                 vdual,vdual1,vdual2)
!--   INIT DOM NUMT
          else if((imot(3).eq.4).and.(mot(3)(1:4).eq.'numt')) then
            call c_intn(mot,imot,nmot)
          else
            call synterr(mot,imot,3,cb)
          endif
!--   INIT BOUNDARY
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'boundary')) then
!--   INIT BOUNDARY BASIC
          if((imot(3).eq.5).and.(mot(3)(1:5).eq.'basic')) then
            call c_inbdb( &
                 mot,imot,nmot, &
                 ncbd,ncin,bceqt)
!--   INIT BOUNDARY NORM
          elseif((imot(3).eq.4).and.(mot(3)(1:4).eq.'norm')) then
            call c_inbdn( &
                 mot,imot,nmot, &
                 x,y,z, &
                 sn, &
                 ncbd,nxn,nyn,nzn, &
                 tn1,tn2,tn3,tn4,tn5,tn6, &
                 tn7,tn8,tn9)
!--   INIT BOUNDARY COIN
          elseif((imot(3).eq.4).and.(mot(3)(1:4).eq.'coin')) then
            call c_inbdc( &
                 mot,imot,nmot, &
                 exs1,exs2, &
                 x,y,z, &
                 ncbd,ncin,mnc)
!--   INIT BOUNDARY NONCOIN
          elseif((imot(3).eq.7).and.(mot(3)(1:7).eq.'noncoin')) then
!            call c_inbdr( &
!                 mot,imot,nmot, &
!                 exr1,exr2,exs1,exs2, &
!                 x,y,z, &
!                 ncbd,mnr,xnr,ynr,znr, &
!                 tm1 ,tm2 ,tm3 ,tm4 ,tm5 ,tm6 ,tm7 ,tm8 ,tm9 ,tm10, &
!                 tm11,tm12,tm13,tm14,tm15,tm16,tm17,tm18,tm19,tm20, &
!                 tm21,tm22,tm23,tm24,tm25,tm26,tm27,tm28,tm29,tm30, &
!                 tm31,tm32,tm33,tm34,tm35,tm36,tm37,tm38,tm39,tm40, &
!                 tm41,tm42,tm43,tm44, &
!                 tn1,tn2,tn3)
          else
            call synterr(mot,imot,3,cb)
          endif
        else
            call synterr(mot,imot,2,cb)
        endif
!--   SAVE
      elseif((imot(1).eq.4).and.(mot(1)(1:4).eq.'save')) then
!--   SAVE DOM
        if((imot(2).eq.3).and.(mot(2)(1:3).eq.'dom')) then
!--   SAVE DOM XYZ
          if((imot(3).eq.3).and.(mot(3)(1:3).eq.'xyz')) then
            do l=1,lzx
            call c_svgr( &
                 mot,imot,nmot, &
                 l,x,y,z, &
                 tn1,tn2,tn3)
            enddo
!--   SAVE DOM CSVMUT
          elseif((imot(3).eq.6).and.(mot(3)(1:6).eq.'csvmut')) then
            do l=1,lzx
            call c_svfw( &
                 mot,imot,nmot, &
                 l,v,mut,utau, &
                 ncin,ncbd, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
            enddo
          else
            call synterr(mot,imot,3,cb)
          endif
!--   SAVE DUAL
         elseif((imot(2).eq.4).and.(mot(2)(1:4).eq.'dual')) then
          do l=1,lzx
           call svdual(l,vdual,vdual1,vdual2)
          enddo
!--   SAVE BOUNDARY
        elseif((imot(2).eq.8).and.(mot(2)(1:8).eq.'boundary')) then
!--   SAVE BOUNDARY 'CREATION'
          if(nmot.eq.2) then
           do mfbi=1,mtbx
            call c_svbd( &
                 mot,imot,nmot, &
                 mfbi, &
                 ncbd)
           enddo
!--   SAVE BOUNDARY BASIC
          else if((imot(3).eq.5).and.(mot(3)(1:5).eq.'basic')) then
           do mfbi=1,mtbx
            call c_svbdb( &
                 mot,imot,nmot, &
                 mfbi, &
                 ncin)
           enddo
!--   SAVE BOUNDARY NORM
          elseif((imot(3).eq.4).and.(mot(3)(1:4).eq.'norm')) then
           do mfn=1,mtnx
           mfbi=nfbn(mfn)
            call c_svbdn( &
                 mot,imot,nmot, &
                 mfbi, &
                 nxn,nyn,nzn)
           enddo
!--   SAVE BOUNDARY COIN
          elseif((imot(3).eq.4).and.(mot(3)(1:4).eq.'coin')) then
           do mfc=1,mtcx
            mfbi=nfbc(mfc)
            call c_svbdc( &
                 mot,imot,nmot, &
                 mfbi, &
                 mnc)
      enddo
!--   SAVE BOUNDARY NONCOIN
          elseif((imot(3).eq.7).and.(mot(3)(1:7).eq.'noncoin')) then
           do mfr=1,mtrx
            mfbi=nfbr(mfr)
            do img=1,lgx
!             call c_svbdr( &
!                 mot,imot,nmot, &
!                 mfbi,img, &
!                 mnr,xnr,ynr,znr)
            enddo
           enddo
          else
            call synterr(mot,imot,3,cb)
          endif
        else
            call synterr(mot,imot,2,cb)
        endif
!--   SET
      elseif((imot(1).eq.3).and.(mot(1)(1:3).eq.'set')) then
!--   SET ENV_CPFW
        if((imot(2).eq.8).and.(mot(2)(1:8).eq.'env_cpfw')) then
          call c_secpfw(mot,imot,nmot)
        else
          call synterr(mot,imot,2,cb)
        endif
      else
         call synterr(mot,imot,1,cb)
      endif
!
      enddo
!
      end
