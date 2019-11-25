PROC IMPORT OUT= WORK.zhanmin 
            DATAFILE= "C:\Users\838035\Documents\ErasmusMC_datasets\CPO\
Zhanmin Lin\Zhanmin.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data zhanmin;
set zhanmin;
ID = cat("Molecule", Index);
drop T_No_ C_No_ I_No_ Index;
run;

/* simple statistics*/
proc means data =  zhanmin; 
run;

/*correlations*/
proc corr data=zhanmin out=corr_zhanmin;
run;

/*Standardize the data */
proc standard data = zhanmin mean=0 std=1 out=stan_zhanmin;
run;

/* Make principal components */
proc princomp data=stan_zhanmin std out=prin_zhanmin;
run;

/* We will choose 4 principal components. ALmost 73% variance */
data prin_zhanmin;
set prin_zhanmin;
drop Prin5-Prin41;
run;

proc corr data=prin_zhanmin out=loadings_princomp;
run;

data loadings_princomp;
  set loadings_princomp(where=(_type_='CORR'));
  format _numeric_ 4.2;
  keep _NAME_ Prin1-Prin4;
run;

proc sgplot data=prin_zhanmin;
	title 'Prin1 vs Prin2';
    scatter x=Prin1 y=Prin2;
run;

proc sgplot data=prin_zhanmin;
	title 'Prin1 vs Prin3';
    scatter x=Prin1 y=Prin3;
run;

proc sgplot data=prin_zhanmin;
	title 'Prin2 vs Prin3';
    scatter x=Prin2 y=Prin3;
run;

/* Cluster analysis */
ods pdf file="C:\Users\838035\Dropbox\PhD\src\cpo\Zhanmin Lin\cluster.pdf";
ods graphics on;
proc cluster data=prin_zhanmin method=ward ccc pseudo print=25 outtree=Tree RSQUARE RMSSTD;
	var Prin1-Prin4;
	id ID;
run;
ods graphics off;
ods pdf close;

proc tree data=Tree out=clus_zhamin_ward nclusters=15;
id ID;
copy Prin1-Prin4;
run;

proc sgplot data=clus_zhamin_ward;
	title 'Prin1 vs Prin2';
    scatter x=Prin1 y=Prin2 /group=cluster;
run;

proc sgplot data=clus_zhamin_ward;
	title 'Prin1 vs Prin3';
    scatter x=Prin1 y=Prin3/group=cluster;
run;

proc sgplot data=clus_zhamin_ward;
	title 'Prin2 vs Prin3';
    scatter x=Prin2 y=Prin3/group=cluster;
run;

/* K-Means */
proc fastclus data=prin_zhanmin maxc=11 maxiter=100 distance out=clus_zhanmin_kmeans;
   var Prin1-Prin4;
run;

ods pdf file="C:\Users\838035\Dropbox\PhD\src\cpo\Zhanmin Lin\11 clusters.pdf";
proc sgplot data=clus_zhanmin_kmeans;
	title 'Micro-Level Movement (Prin1) vs Macro-Level Movement (Prin2)';
    scatter x=Prin1 y=Prin2 /group=cluster;
	xaxis label='Micro-Level Movement (Prin1)';
	yaxis label='Macro-Level Movement (Prin2)';
run;

proc sgplot data=clus_zhanmin_kmeans;
	title 'Micro-Level Movement (Prin1) vs Movement ellipsoid features (Prin3)';
    scatter x=Prin1 y=Prin3/group=cluster;
	xaxis label='Micro-Level Movement (Prin1)';
	yaxis label='Movement ellipsoid features (Prin3)';
run;

proc sgplot data=clus_zhanmin_kmeans;
	title 'Macro-Level Movement (Prin1) vs Movement ellipsoid features (Prin3)';
    scatter x=Prin2 y=Prin3/group=cluster;
	xaxis label='Macro-Level Movement (Prin2)';
	yaxis label='Movement ellipsoid features (Prin3)';
run;
ods pdf close;

/********************/
data zhanmin_log10msd;
set zhanmin;
log10msd = log10(abs(MSD));
run;

data zhanmin_log10msd_filtered;
set zhanmin_log10msd;
where log10msd>0.0062;
drop MSD;
run;

/*Standardize the data */
proc standard data = zhanmin_log10msd_filtered mean=0 std=1 out=stan_zhanmin_filtered;
run;

/* Make principal components */
proc princomp data=stan_zhanmin_filtered std out=prin_zhanmin_filtered;
run;

proc cluster data=prin_zhanmin_filtered method=ward ccc pseudo print=25 outtree=Tree_filtered RSQUARE RMSSTD;
	var Prin1-Prin4;
	id ID;
run;

proc tree data=Tree_filtered out=clus_zhamin_ward_filtered nclusters=2;
id ID;
copy Prin1-Prin4;
run;

proc fastclus data=prin_zhanmin_filtered maxc=2 maxiter=100 distance out=clus_zhanmin_kmeans_filtered;
   var Prin1-Prin4;
run;

proc sgplot data=clus_zhanmin_kmeans_filtered;
	title 'Micro-Level Movement (Prin1) vs Macro-Level Movement (Prin2)';
    scatter x=Prin1 y=Prin2 /group=cluster;
	xaxis label='Micro-Level Movement (Prin1)';
	yaxis label='Macro-Level Movement (Prin2)';
run;

proc sgplot data=clus_zhanmin_kmeans_filtered;
	title 'Micro-Level Movement (Prin1) vs Movement ellipsoid features (Prin3)';
    scatter x=Prin1 y=Prin3/group=cluster;
	xaxis label='Micro-Level Movement (Prin1)';
	yaxis label='Movement ellipsoid features (Prin3)';
run;

proc sgplot data=clus_zhanmin_kmeans_filtered;
	title 'Macro-Level Movement (Prin1) vs Movement ellipsoid features (Prin3)';
    scatter x=Prin2 y=Prin3/group=cluster;
	xaxis label='Macro-Level Movement (Prin2)';
	yaxis label='Movement ellipsoid features (Prin3)';
run;
