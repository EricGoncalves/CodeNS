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
  integer linx ,intmx
  real reelmx,reelmn
  real pis2,raddeg,degrad
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
