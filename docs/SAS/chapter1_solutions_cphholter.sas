*----------------------------------------------------------------------------------------------------------------------------------;
*----------------------------- CHAPTER 1 - EXERCISE SOLUTIONS - COPENHAGEN HOLTER STUDY -------------------------------------------;
*----------------------------------------------------------------------------------------------------------------------------------;

* The data in wide format is created below;

data wide_data;
input id timeafib timestroke timedeath lastobs;
datalines;
1 . . . 100
2 10 . . 90
3 . 20 . 80
4 15 30 . 85
5 . . 70 70
6 30 . 75 75 
7 . 35 95 95
8 25 50 65 65
;

* To convert the data into long format we must first add indicators of whether the subjects transited to AF, stroke and death 
  during the study. Furthermore, we have to fill in timeafib, timestroke and timedeath for all subjects, such that they correspond to
  the last time where transitions to AF, stroke and death are possible.

* We can now create a long format version of the data, i.e. for each subject all relevant transitions h->j correspond to one
  row in the data set and where
  	trans: transition h -> j
  	Tstart: Time of entry into h,
  	Tstop: Time last seen in h,
  	Status: Transition to j or not at time Tstop;

* The seven possible transitions are:
	trans 1: 0 -> AF
	trans 2: 0 -> Stroke
	trans 3: 0 -> Dead
	trans 4: AF -> Stroke
	trans 5: AF -> Dead
	trans 6: Stroke -> Dead
    trans 7: Stroke -> AF;

data wide_data;
	set wide_data;
	* Introducing indicators for transitions and last time at risk in each state.;
	if timedeath = . then do; death = 0; end; else do death = 1; end;
    if death = 0 then timedeath = lastobs;
    if timestroke = . then do; stroke = 0; timestroke = timedeath; end; else do; stroke = 1; end;
	if timeafib = . then do; afib = 0; timeafib = timedeath; end; else do; afib = 1; end;
run;

proc print data = wide_data; run;


data long_data;
	set wide_data;
	* Path: 0;
	if ((afib = 0) and (stroke = 0) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF;
	if ((afib = 1) and (stroke = 0) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke;
	if ((afib = 0) and (stroke = 1) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF -> Stroke;
	if ((afib = 1) and (stroke = 1) and (death = 0) and (timestroke >= timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timestroke; status = 1; output;
	trans = 5; Tstart = timeafib; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke -> AF;
    if ((afib = 1) and (stroke = 1) and (death = 0) and (timestroke < timeafib)) then do;
    trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timeafib; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timeafib; status = 1; output; end;
	* Path: 0 -> Dead;
	if ((afib = 0) and (stroke = 0) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timedeath; status = 1; output; end;
	* Path: 0 -> AF -> Dead;
	if ((afib = 1) and (stroke = 0) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 1; output; end;
	* Path: 0 -> Stroke -> Death;
	if ((afib = 0) and (stroke = 1) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 1; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF -> Stroke -> Dead;
	if ((afib = 1) and (stroke = 1) and (death = 1) and (timestroke >= timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timestroke; status = 1; output;
	trans = 5; Tstart = timeafib; Tstop = timestroke; status = 0; output; 
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 1; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke -> AF -> Dead;
	if ((afib = 1) and (stroke = 1) and (death = 1) and (timestroke < timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 1; output; 
	trans = 6; Tstart = timestroke; Tstop = timeafib; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timeafib; status = 1; output; end;
	keep id trans Tstart Tstop status;
run;

proc print data = long_data; run;


proc import out = chs_data
    datafile = 'C:/HRfinal/holter/cphholter.csv'
	dbms= csv replace;
	getnames=yes;
run;

* We can now repeat the procedure from above to transform the entire data set into long format;

data chs_long;
	set chs_data;
	if timestroke = . then timestroke = timedeath;
	if timeafib = . then timeafib = timedeath;
	* Path: 0;
	if ((afib = 0) and (stroke = 0) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF;
	if ((afib = 1) and (stroke = 0) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke;
	if ((afib = 0) and (stroke = 1) and (death = 0)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF -> Stroke;
	if ((afib = 1) and (stroke = 1) and (death = 0) and (timestroke >= timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timestroke; status = 1; output;
	trans = 5; Tstart = timeafib; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke -> AF;
    if ((afib = 1) and (stroke = 1) and (death = 0) and (timestroke < timeafib)) then do;
    trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timeafib; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timeafib; status = 1; output; end;
	* Path: 0 -> Dead;
	if ((afib = 0) and (stroke = 0) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timedeath; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timedeath; status = 1; output; end;
	* Path: 0 -> AF -> Dead;
	if ((afib = 1) and (stroke = 0) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 1; output; end;
	* Path: 0 -> Stroke -> Death;
	if ((afib = 0) and (stroke = 1) and (death = 1)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 1; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> AF -> Stroke -> Dead;
	if ((afib = 1) and (stroke = 1) and (death = 1) and (timestroke >= timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timeafib; status = 1; output;
	trans = 2; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 3; Tstart = 0; Tstop = timeafib; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timestroke; status = 1; output;
	trans = 5; Tstart = timeafib; Tstop = timestroke; status = 0; output; 
	trans = 6; Tstart = timestroke; Tstop = timedeath; status = 1; output; 
    trans = 7; Tstart = timestroke; Tstop = timedeath; status = 0; output; end;
	* Path: 0 -> Stroke -> AF -> Dead;
	if ((afib = 1) and (stroke = 1) and (death = 1) and (timestroke < timeafib)) then do;
	trans = 1; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 2; Tstart = 0; Tstop = timestroke; status = 1; output;
	trans = 3; Tstart = 0; Tstop = timestroke; status = 0; output;
	trans = 4; Tstart = timeafib; Tstop = timedeath; status = 0; output;
	trans = 5; Tstart = timeafib; Tstop = timedeath; status = 1; output; 
	trans = 6; Tstart = timestroke; Tstop = timeafib; status = 0; output; 
    trans = 7; Tstart = timestroke; Tstop = timeafib; status = 1; output; end;
	keep id trans Tstart Tstop status;
run;

proc print data = chs_long; run;
