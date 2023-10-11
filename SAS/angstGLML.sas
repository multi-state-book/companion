data angst;
infile 'angstGLML.txt' firstobs=2;
input id prev stop status bip;
run;

proc phreg data=angst;
class bip(ref='0');
model stop*status(0)=bip/entry=prev eventcode=1 rl;
run;

data angstML; set angst;
if status ne 2 then output;
if status = 2 then do; status=1; output; 
                       prev=stop; stop=stop+0.5; status=2; output;
end; 
run;

proc phreg data=angstML;
class bip(ref='0');
model stop*status(0)=bip/entry=prev eventcode=1 rl;
run;
