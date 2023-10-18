
libname h "d:\HRbog";
data pbc3;
  set h.pbc3;
  run;

data pbc3; set pbc3;
log2bili=log2(bili);
run;

/* bootstrappes - 120123: 200 ændret til 1000*/

data bootpbc;
do sampnum = 1 to 1000; /* nboot=1000*/
do i = 1 to 349; /*nobs=349*/
x=round(ranuni(0)*349); /*nobs=349*/
set pbc3
point=x;
output;
end;
end;
stop;
run;

/* Arealer under trappekurver */

%macro areastepby(data,byvar,beh,grp,tid,y,tau);

data select; set &data;
where &beh=&grp;
run;

data select; set select;

by &byvar;
retain mu oldt oldy;

if first.&byvar then do oldt=0; oldy=1; mu=0;  end;

if &tid>&tau then do;
&tid=&tau; &y=oldy; end;

if not first.&byvar then mu+oldy*(&tid-oldt);
if last.&byvar then do;
if &tid<&tau then mu+(&tau-&tid)*&y; end;
oldy=&y; oldt=&tid;
run;

data sidste; set select;

by  &byvar;
if last.&byvar;
run;

proc print data=sidste; run;

%mend areastepby;


proc phreg data=bootpbc;
by sampnum;
model followup*status(0)=;
strata tment;
baseline out=survdat survival=km / method=pl;
run;

/* Restricted mean v.hj. af MACRO Tab 4.1 */

%areastepby(survdat,sampnum,tment,0,followup,km,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(survdat,sampnum,tment,1,followup,km,3);

proc means data=sidste stddev mean;
var mu;
run;

/* predikeret surv given covariates */

data cov;
tment=0; alb=38; log2bili=log2(45); output;
tment=1; alb=38; log2bili=log2(45); output;
run;

proc phreg data=bootpbc;
by sampnum;
class tment (ref='0');
model followup*status(0)=tment alb log2bili/rl;
baseline out=predsurv survival=surv covariates=cov/ method=breslow;
run;

%areastepby(predsurv,sampnum,tment,0,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(predsurv,sampnum,tment,1,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

data cov;
tment=0; alb=20; log2bili=log2(90); output;
tment=1; alb=20; log2bili=log2(90); output;
run;

proc phreg data=bootpbc;
by sampnum;
class tment (ref='0');
model followup*status(0)=tment alb log2bili/rl;
baseline out=predsurv2 survival=surv covariates=cov/ method=breslow;
run;

%areastepby(predsurv2,sampnum,tment,0,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(predsurv2,sampnum,tment,1,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

/* g-formel */

proc phreg data=bootpbc;
by sampnum;
class tment (ref='0');
model followup*status(0)=tment alb log2bili/rl;
baseline out=gsurv survival=surv stderr=se/ method=breslow diradj group=tment;
run;

%areastepby(gsurv,sampnum,tment,0,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(gsurv,sampnum,tment,1,followup,surv,3);

proc means data=sidste stddev mean;
var mu;
run;

/* Find S(2 yr) og beregn bootstrap mean og SD */

data s20; set gsurv;
where tment=0;
by sampnum;
retain told sold seold;
if followup>=2 and told<= 2 then do;
output;
end;
told=followup; sold=surv; seold=se;
run;

data s20; set s20;
var=seold*seold;
s0=sold;
proc means mean stddev;
var var sold;
run;

data s21; set gsurv;
where tment=1;
by sampnum;
retain told sold seold;
if followup>=2 and told<= 2 then do;
output;
end;
told=followup; sold=surv; seold=se;
run;

data s21; set s21;
var=seold*seold;
s1=sold;
proc means mean stddev;
var var sold;
run;

data s2; merge s20 s21; by sampnum;
diff=s1-s0;
proc means mean stddev;
var diff;
run;


/* Comp. risks. years lost */

proc phreg data=bootpbc;
by sampnum;
model followup*status(0)=/eventcode=1;
strata tment;
baseline out=cuminc1 cif=cif1 stdcif=std1;
run;



/* Tab 4.3 */

%areastepby(cuminc1,sampnum,tment,0,followup,cif1,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(cuminc1,sampnum,tment,1,followup,cif1,3);

proc means data=sidste stddev mean;
var mu;
run;

proc phreg data=bootpbc;
by sampnum;
model followup*status(0)=/eventcode=2;
strata tment;
baseline out=cuminc2 cif=cif2;
run;

/* Tab 4.3 */

%areastepby(cuminc2,sampnum,tment,0,followup,cif2,3);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby(cuminc2,sampnum,tment,1,followup,cif2,3);

proc means data=sidste stddev mean;
var mu;
run;

/* Fine-Gray g-formel t=2 */

data cov;
input sex tment age alb log2bili;
datalines;
0 0 0 0 0
;
run;

/* På originale data, baseline til tid 2 aflæses til 0.7358654259  for transpl og 0.0003138795 for død,
    betaer aflæses ligeledes */

proc phreg data=pbc3 outest=beta1orig;
model followup*status(0)=sex tment age log2bili alb/eventcode=1;
baseline out=cuminc1orig covariates=cov cif=cif1;
run;


data pbc3g1; set pbc3;
cold=0.7358654259;
betatreat=-0.408624125; betasex=0.0916082925; betaage=-0.075026871; betaalb=-0.06967979;
betalog2bili=0.6192339539;
pred0=1-(1-cold)**exp(betatreat*0+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
pred1=1-(1-cold)**exp(betatreat*1+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
run;

proc means data=pbc3g1 mean;
var pred0 pred1;
run;



proc phreg data=pbc3 outest=beta2orig;
model followup*status(0)=sex tment age log2bili alb/eventcode=2;
baseline out=cuminc2orig covariates=cov cif=cif2;
run;


data pbc3g2; set pbc3;
cold=0.0003138795;
betatreat=-0.352897932; betasex=0.414825845; betaage=0.0874743944; betaalb=-0.061217334;
betalog2bili=0.6156589227;
pred0=1-(1-cold)**exp(betatreat*0+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
pred1=1-(1-cold)**exp(betatreat*1+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
run;

proc means data=pbc3g2 mean;
var pred0 pred1;
run;

/* Nu på boot data */

proc phreg data=bootpbc outest=beta1;
by sampnum;
model followup*status(0)=sex tment age log2bili alb/eventcode=1;
baseline out=cuminc1 covariates=cov cif=cif1;
run;


/* Find F1(2 yr) og beregn bootstrap mean og SD */

data c21; set cuminc1;
by sampnum;
retain told cold;
if followup>=2 and told<= 2 then do;
output;
end;
told=followup; cold=cif1;
run;

data c21rens; set c21;
drop tment sex alb age log2bili told cif1 followup;
run;

data betafin1; set beta1;
betatreat=tment; betasex=sex; betaage=age; betaalb=alb; betalog2bili=log2bili;
keep sampnum betatreat betasex betaage betaalb betalog2bili;
run;

data gformel1; merge betafin1 c21rens; by sampnum; run;

data gfinal1; merge gformel1 bootpbc; by sampnum; run;

data gfinal1; set gfinal1;
pred0=1-(1-cold)**exp(betatreat*0+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
pred1=1-(1-cold)**exp(betatreat*1+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
diff=pred1-pred0;
run;

proc means data=gfinal1;
by sampnum;
var pred0;
output out=gresults10 mean=xbar0;
run;

proc means data=gresults10 mean stddev;
var xbar0;
run;

/* I 9 bootstrap datasæt er prediktionen=1 */

proc means data=gresults10 mean stddev;
where xbar0<1;
var xbar0;
run;

proc means data=gfinal1;
by sampnum;
var pred1;
output out=gresults11 mean=xbar1;
run;

proc means data=gresults11 mean stddev;
var xbar1;
run;


proc means data=gresults11 mean stddev;
where xbar1<1;
var xbar1;
run;

/* treatment diff in 'good' data sets */

proc means data=gfinal1;
by sampnum;
var diff;
output out=gresults1d mean=xbard;
run;

proc means data=gresults1d mean stddev;
where xbard ne 0;
var xbard;
run;

/* nu død uden transpl. */

proc phreg data=bootpbc outest=beta2;
by sampnum;
model followup*status(0)=sex tment age log2bili alb/eventcode=2;
baseline out=cuminc2 covariates=cov cif=cif2;
run;

/* Find F2(2 yr) og beregn bootstrap mean og SD */

data c22; set cuminc2;
by sampnum;
retain told cold;
if followup>=2 and told<= 2 then do;
output;
end;
told=followup; cold=cif2;
run;

data c22rens; set c22;
drop tment sex alb age log2bili told cif2 followup;
run;

data betafin2; set beta2;
betatreat=tment; betasex=sex; betaage=age; betaalb=alb; betalog2bili=log2bili;
keep sampnum betatreat betasex betaage betaalb betalog2bili;
run;

data gformel2; merge betafin2 c22rens; by sampnum; run;

data gfinal2; merge gformel2 bootpbc; by sampnum; run;

data gfinal2; set gfinal2;
pred0=1-(1-cold)**exp(betatreat*0+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
pred1=1-(1-cold)**exp(betatreat*1+betasex*sex+betaage*age+betaalb*alb+betalog2bili*log2bili);
diff=pred1-pred0;
run;

proc means data=gfinal2;
by sampnum;
var pred0;
output out=gresults20 mean=xbar0;
run;

proc means data=gresults20 mean stddev;
var xbar0;
run;

/* Her er ingen prediktioner=1 */

proc means data=gfinal2;
by sampnum;
var pred1;
output out=gresults21 mean=xbar1;
run;

proc means data=gresults21 mean stddev;
var xbar1;
run;


proc means data=gresults1 mean stddev;
where xbar1<1;
var xbar1;
run;



proc means data=gfinal2;
by sampnum;
var diff;
output out=gresults2d mean=xbard;
run;

proc means data=gresults2d mean stddev;
var xbard;
run;

/* Og nu Fig 6.1 - koden burde virke, hvis ptno 415 og  458 erstattes af id 305 og 325.
   variablen fail skal lige tilføjes: */

data pbc3; set pbc3; fail=(status>0); run; 

data minus1; set pbc3; /* drop failure with t~1 yr */
if ptno=415 then delete;
run;

proc lifetest data=pbc3;
time days*fail(0);
survival out=sall;
run;

data sall; set sall; survall=survival; run;

proc lifetest data=minus1;
time days*fail(0);
survival out=sminus1;
run;

data ekstra;
input days survival;
datalines;
366 0.9225440684
;
run;

data sminus1ny; set sminus1 ekstra; 
survny=survival;
proc sort; by days;
run;

data final1; merge sall sminus1ny; by days; 
followup=days/365;
pseudo1=349*survall-348*survny;
run;

proc gplot data=final1;
plot pseudo1*followup;
run;

data minus2; set pbc3; /* drop cens with t~1 yr */
if ptno=458 then delete;
run;

data ekstra;
input days survival;
datalines;
365 0.9225440684
;
run;

proc lifetest data=minus2;
time days*fail(0);
survival out=sminus2;
run;


data sminus2ny; set sminus2 ekstra; 
survny=survival;
proc sort; by days;
run;

data final2; merge sall sminus2ny; by days; 
followup=days/365;
pseudo2=349*survall-348*survny;
run;

proc gplot data=final2;
plot pseudo2*followup;
run;

data final; merge final1 final2; by days; 
en=1; nul=0;run;

proc gplot data=final;
plot pseudo1*followup pseudo2*followup en*followup nul*followup
      /overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 6 by 1 minor=none label=('Years');
axis2 order=-0.1 to 1.1 by 0.1 minor=none label=(a=90 'Pseudo-values');
symbol1  v=none i=join c=blue;
symbol2  v=none i=join c=red;
symbol3 v=none i=join c=black;
symbol4 v=none i=join c=black;
run;
