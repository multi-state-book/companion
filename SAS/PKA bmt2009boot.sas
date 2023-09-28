options nocenter ls=80 ps=3000;

data one; infile 'bmt-lida-data.txt';
 input id team disease age dnrage sex dnrsex karnofpr 
       gsource regimprg prevgvh1 indxtx indxcr1 incr1tx 
       anc500 dytxanc5 agvhd dytxagvh cgvhd intxcgvh
       trm rel lfs intxrel dead intxsurv;

data one; set one;
label disease='Disease type';
label age='Patient Age';
label dnrage='Donor Age';
label sex='Patient Gender';
label dnrsex='Donor Gender';
label karnofpr='Karnofsky Score at Pre-TX';
label gsource='Graftype';
label regimprg='Conditioning Regimen';
label prevgvh1='GVHD Prophylaxis';
label indxtx='Interval DX-TX in Mons';
label indxcr1='Interval DX-CF1 in Mons';
label incr1tx='Interval CR1-TX in mons';
label anc500='Achieve ANC 500';
label dytxanc5='Interval TX-ANC500 in Days';
label agvhd='Develop AGVHD';
label dytxagvh='Interval TX-AGVHD in Days';
label cgvhd='Develop CGVHD';
label intxcgvh='Interval TX-CGVHD in Mons';
label trm='Indicator of TRM';
label rel='Indicator of Relapse';
label lfs='Indicator of LFS (Leukemia-Free-Survival)';
label intxrel='Interval TX-Relapse in Mons';
label dead='Indicator of Death';
label intxsurv='Interval TX-Death in Mons';

proc format;
  value disease
    10= '10 AML'
    20= '20 ALL';
  value sex
    1 = '1 Male'
    2 = '2 Female';
  value dnrsex
    1 = '1 Male'
    2 = '2 Female';
  value karnofpr
   -9 = '-9 Unknown';
  value gsource
    1='1 BM' 
    2='2 PB/PB+BM';
  value regimprg
    1 = '1 CY + TBI +/- oth'
    2 = '2 TBI + other'
    3 = '3 Busulf+CY+/-oth'
    9 = '9 Other/Unknown';
  value prevgvhf
    1 ='1 mtx +/- other'
    2 ='2 csa +/- other'
    3 ='3 mtx + csa +/- other'
    4 ='4 tdep +/- other'
    9 ='9 Other/Unknown';
  value anc500f
    1 = '1 Present'
    0 = '0 Absent'
   -9 = '-9 Unknown';
  value agvhd
    1 = '1 Moderate to Severe'
    0 = '0 None/Mild'
   -9 = '-9 Unknown';
  value cgvhd
    1 = '1 Mild to Severe'
    0 = '0 None'
   -9 = '-9 Unknown';
  value event
    1='1 Yes'
    0='0 No';
run;

/* Tab 1.4 (mm) */

data one; set one;
state0=rel+2*trm;
run;


/* State probabilities */

proc phreg data=one; /* Relapse-free surv */
model intxrel*state0(0)=;
baseline out=surv survival=km;
run;

proc phreg data=one; /* Relapse */
model intxrel*state0(0)=/eventcode=1;
baseline out=cif1 cif=cif1;
run;

proc phreg data=one; /* Death in remission */
model intxrel*state0(0)=/eventcode=2;
baseline out=cif2 cif=cif2;
run;

proc phreg data=one; /* Overall surv. */
model intxsurv*dead(0)=/eventcode=1;
baseline out=dead cif=cif23;
run;

/* We need the same time variable for all prob's */

data dead; set dead; time=intxsurv; run;
data surv; set surv; time=intxrel; run;
data cif1; set cif1; time=intxrel; run;
data cif2; set cif2; time=intxrel; run;
data all; merge surv cif1 cif2 dead; by time; run;

data allrev; set all;
by time;
retain last1 last2 last3 last4;
if km=. then rfs=last1; if km ne . then rfs=km; 
if cif1=. then c1=last2; if cif1 ne . then c1=cif1;
if cif2=. then c2=last3; if cif2 ne . then c2=cif2;
if cif23=. then c23=last4; if cif23 ne . then c23=cif23;
output;
last1=rfs; last2=c1; last3=c2; last4=c23;
run;

data allrev; set allrev;
q0=rfs; q2=c2; q3=c23-c2; q1=c1-q3; sum=q0+q1+q2+q3; prev=q1/(q0+q1); tment=0;
run;

/* Fig 4.13 */

proc gplot data=allrev;
plot prev*time q1*time/overlay haxis=axis1 vaxis=axis2;
axis1 order=0 to 120 by 10 minor=none label=('Months');
axis2 order=0 to 0.05 by 0.01 minor=none label=(a=90 'Relapse prev. and prob.');
symbol1  v=none i=stepjl c=blue;
symbol2  v=none i=stepjl c=red;
run;

/* Expected time spent in states, Sect. 4.1 */

%areastep(allrev,tment,0,time,q1,120);
%areastep(allrev,tment,0,time,q0,120);
%areastep(allrev,tment,0,time,q2,120);
%areastep(allrev,tment,0,time,q3,120);
%areastep(allrev,tment,0,time,c23,120);

/* Bootstrap */
data bootbmt;
do sampnum = 1 to 1000; /* nboot=1000*/
do i = 1 to 2009; /*nobs=2009*/
x=round(ranuni(0)*2009); /*nobs=2009*/
set one
point=x;
output;
end;
end;
stop;
run;


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


proc phreg data=bootbmt; /* Relapse-free surv */
by sampnum;
model intxrel*state0(0)=;
baseline out=surv survival=km;
run;

proc phreg data=bootbmt; /* Relapse */
by sampnum;
model intxrel*state0(0)=/eventcode=1;
baseline out=cif1 cif=cif1;
run;

proc phreg data=bootbmt; /* Death in remission */
by sampnum;
model intxrel*state0(0)=/eventcode=2;
baseline out=cif2 cif=cif2;
run;

proc phreg data=bootbmt; /* Overall surv. */
by sampnum;
model intxsurv*dead(0)=/eventcode=1;
baseline out=dead cif=cif23;
run;

data dead; set dead; time=intxsurv; drop intxsurv;  run;  
data surv; set surv; time=intxrel; drop intxrel; run;
data cif1; set cif1; time=intxrel; drop intxrel; run;
data cif2; set cif2; time=intxrel; drop intxrel; run;
data all; merge surv cif1 cif2 dead ; by sampnum time; run;

data allrev; set all;
by sampnum time; 
retain last1 last2 last3 last4;
if km=. then rfs=last1; if km ne . then rfs=km; 
if cif1=. then c1=last2; if cif1 ne . then c1=cif1;
if cif2=. then c2=last3; if cif2 ne . then c2=cif2;
if cif23=. then c23=last4; if cif23 ne . then c23=cif23;
output;
last1=rfs; last2=c1; last3=c2; last4=c23;
run;

data allrev; set allrev;
q0=rfs; q2=c2; q3=c23-c2; q1=c1-q3; sum=q0+q1+q2+q3; prev=q1/(q0+q1); tment=0;
run;

%areastepby(allrev,sampnum,tment,0,time,q0,120);

proc means data=sidste stddev mean;
var mu;
run;

/* macro skal ændres for cuminc (skal starte i 0) */

%macro areastepby0(data,byvar,beh,grp,tid,y,tau);

data select; set &data;
where &beh=&grp;
run;

data select; set select;

by &byvar;
retain mu oldt oldy;

if first.&byvar then do oldt=0; oldy=0; mu=0;  end;

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

%areastepby0(allrev,sampnum,tment,0,time,q1,120);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby0(allrev,sampnum,tment,0,time,q2,120);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby0(allrev,sampnum,tment,0,time,q3,120);

proc means data=sidste stddev mean;
var mu;
run;

%areastepby0(allrev,sampnum,tment,0,time,c23,120);

proc means data=sidste stddev mean;
var mu;
run;
