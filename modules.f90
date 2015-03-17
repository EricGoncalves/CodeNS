module para_fige
  implicit none
  integer          ::   ista,    lg,lgcmdx,  lsta,    lt
  integer          ::     lz,   mtb,   mtt,  ndir,  neqt
  integer          ::   nind,   nmx,  nobj,  nsta
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
  integer          ::    ip00,   ip11,   ip12,   ip13,   ip21
  integer          ::    ip31,   ip40,   ip41,   ip42,   ip43
  integer          ::    ip44,   ip60,  kdimg,  kdimk,  kdimv
  integer          :: mdimtbf,mdimtcf,mdimtnf,mdimtrf, mdimub
  integer          :: ndimctf,ndimnts,ndimntu, ndimub,   nvar
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
  parameter(ip11=nint((1.+kdimg*ccg)*ndimctf))
  parameter(ip12=kdimv*(ip11-1)+1)
  parameter(ip13=kdimk*(ip11-1)+1)
  parameter(ip21=nint((1.+kdimg*cng)*ndimnts+ndimntu))
  parameter(ip31=nint(1.+3.*(1.+kdimg*cng)*ndimnts))
  parameter(ip40=mdimub)
  parameter(ip41=nint((1.+kdimg*cfg)*mdimtbf))
  parameter(ip42=nint((1.+kdimg*cfg)*mdimtnf))
  parameter(ip43=nint((1.+kdimg*cfg)*mdimtcf))
  parameter(ip44=nint((1.+kdimg*cfg)*mdimtrf))
  parameter(ip60=nvar)
end module para_var
!
module boundary
  use para_fige
implicit none
  integer          :: crotr,imaxb,iminb,jmaxb,jminb
  integer          ::  kexl,kmaxb,kminb,  lbd,lbdko
  integer          ::  mdnc,  mmb,  mpb,  mpc, mper
  integer          ::   mpn,  mpr,  nba,  nbd, nbdc
  integer          :: nbdko, ndcc, ndlb, ndrr, nfba
  integer          ::  nfbc, nfbn, nfbr, nfei,srotr
  double precision :: bc
  character(len=4) :: cl
  character(len=2) :: indfl
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
  integer          ::      id1,     id2,     ii1,     ii2,     jd1
  integer          ::      jd2,     jj1,     jj2, kcaldis,     kd1
  integer          ::      kd2, kecrdis,     kk1,     kk2, klecdis
  integer          ::      kvn,  lbdrat,     lgx,     lzx, mdimtbx
  integer          ::  mdimtcx, mdimtnx, mdimtrx, mdimubx,    mtax
  integer          ::     mtbx,    mtcx,    mtnx,    mtrx,  nbdrat
  integer          :: ndimctbx,ndimntbx, ndimubx,   neqtx,     nnc
  integer          ::     nnfb,     nnn,  npbrat,     npc,    npfb
  integer          ::      npn,   nptot
  dimension npn(lt),nnn(lt),npc(lt),nnc(lt),npfb(lt),nnfb(lt)
  dimension id1(lt),ii1(lt),ii2(lt),id2(lt),jd1(lt),jj1(lt), &
            jj2(lt),jd2(lt),kd1(lt),kk1(lt),kk2(lt),kd2(lt)
  dimension nbdrat(lz),npbrat(lz),lbdrat(mtb)
end module maillage
!
module definition
  use para_fige
implicit none
  integer          :: klomg
  double precision ::    aa1,   anz,   dnz,   ha1,   omg
  double precision ::    pa1, perio,   pnz,protat,ptrans
  double precision ::    rnz,  roa1,  ronz,   ta1,   tnz
  double precision ::  varst
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
  integer          :: intmx, linx
  double precision :: degrad,  pis2,raddeg,reelmn,reelmx
  data linx/132/
!  data intmx/999999/
!  data reelmx/999999999./
!  data reelmn/1.e-30/
!  data pis2/1.570796327/
!  data raddeg/57.29577951/
!  data degrad/0.01745329252/
end module constantes
!
module proprieteflu
implicit none
  double precision ::    cp,   cv,  gam, gam1, gam2
  double precision ::  gam3, gam4, gam5,pinfl,   pr
  double precision ::   prt,   ql,   rd,reynz,  rgp
end module proprieteflu
!
module kcle
  use para_fige
implicit none
  integer          ::      kanz,  kconfig,     kcte,  kdiscsv,     kdnz
  integer          ::   kdt1min,   kequat,     keta,     kgam, kicychr0
  integer          ::  kicytur0,    kilim, kischema,     kkcg,    kkdtl
  integer          ::  kkdualns,    kkexl,    kkfmg,     kki2,     kki4
  integer          ::    kklomg,     kkmf,   kkprec,   kkvisq,     kkvn
  integer          ::      klgx,    klmax,     klzx, kmdimtbx, kmdimtcx
  integer          ::  kmdimtnx, kmdimtrx, kmdimubx,    kmtax,    kmtbx
  integer          ::     kmtcx,    kmtnx,    kmtrx,   kmuscl,     knba
  integer          ::  kncychro,  kncycle, kncyexpl, kncyresi, kncysave
  integer          ::  kncyturb,kndimctbx,kndimntbx, kndimubx,   kneqtx
  integer          ::      knfi,   kniter,   knitur,    knumt,     komg
  integer          ::  kpctvort,   kperio,   kpinfl,     kpnz,      kpr
  integer          ::      kprt,      kql,      krd,   kreynz,     krnz
  integer          ::     kronz,  ktitrt1,     ktnz,     ktol,   ktolke
  integer          ::    kvarst,      kxk
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
  integer          :: icychr0,   ilim,ischema,    kcg,   kdtl
  integer          :: kdualns,kdualto,   kfmg,  klroe,    kmf
  integer          ::   kprec, ktrans,  kvisq,  lacou,   lmax
  integer          ::     mgl,  muscl,ncychro, ncycle,ncyexpl
  integer          :: ncyresi,ncysave,    nfi,  niter,  nitur
  integer          ::    numt
  double precision ::      al0,     be0,     cga,     cte,  dt1min
  double precision ::   epsroe,     eta,    freq,     ki2,     ki4
  double precision ::   resite,  reske1,reskeite,  resno1,   rki2t
  double precision ::    rki4t,     rm0,     tol,   tolke,  vrtalp
  double precision ::    vrtcz,  vrtlre,  vrtmac,  vrtxre,  vrtzre
  double precision ::       x0,      xk,      y0,      z0
   character(len=4) :: discsv
  dimension ki2(lt),ki4(lt)
  dimension kmf(lt),lmax(lt)
  dimension ncycle(lg) 
  dimension eta(lt)
end module schemanum
!
module modeleturb
  use para_fige
implicit none
  integer          ::  icytur0, imxclko,  kclkep,   kcmut,  kcutke
  integer          ::   kditur,   keasm,  keinit,   kesst, kfludis
  integer          ::  kfracom,   kinke,  komsst,  kparoi,  ksecmb
  integer          ::  ktransi,   kutau,   kwsst,  lparoi,  naprng
  integer          ::   ncycke, ncycrac,ncytuke0, ncyturb
  double precision ::   alfae,  alfak, allfa0,allfae0,   beta
  double precision ::   beta1,  beta2,  betae,  betas, bethae
  double precision ::   bheta,    cb1,    cb2,    cc1,    cc2
  double precision ::    ccmu,    ce2,   ceta,    cgl,   cke1
  double precision ::    cke2,  cklb1,  ckle2,  cllog,    cmu
  double precision ::   cmukl,    ct1,    ct2,    ct3,    ct4
  double precision ::     cv1,    cw1,    cw2,    cw3, drtrac
  double precision ::    epse,   epsk, epspid,epstaud,epsvord
  double precision ::   kappa, okappa,pctvort,rapvisq,  rbeta
  double precision ::  rkplus, roeinf, rokinf,  romeg,    rrk
  double precision ::   rtrac,   sige,   sigk,  sigma, sigma1
  double precision ::  sigma2, sigmab, sigmak, sigmal,  sigmd
  double precision ::  sigme1, sigme2, sigmeb,  sigmk,  sigmw
  double precision :: utaumin,   vkar,  wsig1,  wsig2, xkappa
  double precision ::     yp0
  dimension pctvort(lt)
end module modeleturb
!
module sortiefichier
  use para_fige
implicit none
  integer          ::    don1,    imp,  inia1,  inig1,   kdac
  integer          ::   kdacf,   kdav,   kdgc,  kdgcf,   kdgv
  integer          ::     kfa,    kfb,    kfc,    kfi,    kfn
  integer          ::     kfr,   kimp,   kres,  kvglo,    lec
  integer          :: lsortie,  nbfll,  nfreq, nmfint,    out
  integer          ::     sec,   sor1,   sor2,   sor3,  sorf1
  integer          ::   sorf2
  double precision :: alpha0, beta0,p0spi0,q0spi0,  sref
  double precision ::     v0, xlref,  xref,  yref,  zref
  data kimp/3/
  data lec,imp,out,sec,sor1,sor2,kfa,sor3/11,12,13,14,15,16,17,18/
  data kdgv,kdav,kdgc,kdac/21,22,23,24/
!  data inig1,kfi,kfb,kfn,kfc,kfr,kdgcf,kdacf,kres/31,32,33,34,35,36,37,38,39/
  data don1,inia1,sorf1,sorf2/41,42,43,44/
  dimension nmfint(mtb)
end module sortiefichier
