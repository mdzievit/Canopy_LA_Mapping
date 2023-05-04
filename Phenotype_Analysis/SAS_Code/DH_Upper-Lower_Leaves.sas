proc import 
	out = DHLines
	datafile = "Upper_Lower_LA_2016_2017_2018.txt"
	replace
	DBMS = TAB;
run;
proc sort data = DHLines;
	by Trait Pop Loc Range Pass;
run;

data DHLines_v4;
set DHLines;
	if Gen < 4 then Grp = Gen;
	else Grp = 999;
	if Gen < 4 then Switch = 0;
	else Switch = 1;
run;

/*Use this to obtain variance components for estimating heritability */
proc mixed data = DHLines_v4 method=reml covtest;
	by Trait Pop;
	class Year Loc Range Gen Grp;
	model Value = Grp / DDFM=SATTERTHWAITE;
	random Year Switch*Gen*Grp*Gen Switch*Grp*Gen*Year*Loc Year*Loc Year*Loc*Range Year*Loc*Grp;
	ods output COVPARMS = Var;
run;

ODS HTML CLOSE;
ODS HTML;

/* Use this to obtain BLUEs and BLUPs */
proc mixed data = DHLines_v4 method=reml covtest;
	by Trait Pop;
	class Year Loc Range Gen Grp;
	model Value = Grp Switch*Gen*Grp Year Year*Loc Grp*Year*Loc Switch*Gen*Grp*Year*Loc / Solution DDFM=SATTERTHWAITE;
	random Year*Loc*Range / Solution;
	ods output SolutionR = S_R SolutionF = S_F;
run;

/* Use this to obtain ANOVA table */
proc glm data = DHLines_v4;
By Trait Pop;
class Loc Range Gen Grp;
model Value = Grp Switch*Gen*Grp Loc Grp*Loc Switch*Gen*Loc*Grp Range*Loc;
random Range*Loc;
test h = Switch*Gen*Grp e = Switch*Gen*Loc*Grp;
run;
