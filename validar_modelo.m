%%VALIDACION DEL MODELO NP 15.10.19 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

%aca se puede cambiar el solver
%changeCobraSolver('glpk','all');
%changeCobraSolver('gurobi','all');

%nombre del archivo de salida con los flujos
archivo='res.xlsx';

%modelo de entrada
load('modelo_1481c.mat') %este nombre debe coincidir con el modelo usado en genexp.m

%%%%%%%%%%%%%%%%%%%%%%%%%%[REACCIONES NUEVAS

%sucrose alpha glucohidrolase
model = addReaction(model, 'SAGH','metaboliteList',{'sucr[c]','h2o[c]','f1p[c]','glc_D[c]'},'stoichCoeffList',[-1,-1,1,1])
model.lb(findRxnIDs(model,'SAGH'))=-1000;
model.ub(findRxnIDs(model,'SAGH'))=1000; 

model = addReaction(model, 'SUCRtex','metaboliteList',{'sucr[e]','sucr[p]'},'stoichCoeffList',[-1,1])
model.lb(findRxnIDs(model,'SUCRtex'))=-1000;
model.ub(findRxnIDs(model,'SUCRtex'))=1000; 

model = addReaction(model, 'SUCRtpp','metaboliteList',{'sucr[p]','sucr[c]'},'stoichCoeffList',[-1,1])
model.lb(findRxnIDs(model,'SUCRtpp'))=-1000;
model.ub(findRxnIDs(model,'SUCRtpp'))=1000;

%model = addReaction(model, 'F1PP','metaboliteList',{'h2o[c]','f1p[c]','pi[c]','fru[c]'},'stoichCoeffList',[-1,-1,1,1])
%model.lb(findRxnIDs(model,'F1PP'))=-1000;
%model.ub(findRxnIDs(model,'F1PP'))=1000;

%model = addReaction(model, 'EX_fru(e)','metaboliteList',{'fru[e]'},'stoichCoeffList',[-1])
%model.lb(findRxnIDs(model,'EX_fru(e)'))=0;
%model.ub(findRxnIDs(model,'EX_fru(e)'))=1000;

%printRxnFormula(model, 'rxnAbbrList', 'SAGH')
%printRxnFormula(model, 'rxnAbbrList', 'SUCRtex')
%printRxnFormula(model, 'rxnAbbrList', 'SUCRtpp')

%aca se pueden quitar reacciones
%model = removeRxns(model, {'ACONT'});
%model = removeRxns(model, {'FBA2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE BIOMASA (VER EXCEL)
%biomasa original con hierro
model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','hemeO[c]','sheme[c]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0005,-0.0005,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											
printRxnFormula(model, 'BIOMASS_F1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE GENERALES
model.lb(findRxnIDs(model,'EX_co2(e)'))=-3.1;
model.ub(findRxnIDs(model,'EX_co2(e)'))=-3.1;

model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=0;
model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=1000; 

model.lb(findRxnIDs(model,'EX_nh4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_nh4(e)'))=0; 

model.lb(findRxnIDs(model,'EX_o2(e)'))=-30;
model.ub(findRxnIDs(model,'EX_o2(e)'))=0;

model.lb(findRxnIDs(model,'EX_so4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_so4(e)'))=0;

model.lb(findRxnIDs(model,'EX_mal_L(e)'))=0;
model.ub(findRxnIDs(model,'EX_mal_L(e)'))=0;

%%%
%sacarosa qs=4.01  u=0.36 ??no esta el exchange
%citrate qs=13.45 u=0.5 ->fba: 0.326007841	integr:0.309665385
%gluconate qs=7.67 u=0.46 ->fba: 0.507745773 integr:0.480737413
%glicerol qs=4.39 u=0.18  ->fba:0.031699118	integr:0,029894802
%succinate qs=14.07 u=0.38 ->fba: 0.479345816 integr:0.466889542
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso='VALIDACION SACAROSA (bm emp 0.36)';

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0;

model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

model.lb(findRxnIDs(model,'EX_sucr(e)'))=-4.01; 
model.ub(findRxnIDs(model,'EX_sucr(e)'))=-4.01;

model.lb(findRxnIDs(model,'EX_cit(e)'))=0;
model.ub(findRxnIDs(model,'EX_cit(e)'))=0;

model.lb(findRxnIDs(model,'EX_glyc(e)'))=0;
model.ub(findRxnIDs(model,'EX_glyc(e)'))=0;

model.lb(findRxnIDs(model,'EX_glcn(e)'))=0;
model.ub(findRxnIDs(model,'EX_glcn(e)'))=0;

model.lb(findRxnIDs(model,'EX_succ(e)'))=0;
model.ub(findRxnIDs(model,'EX_succ(e)'))=0;

solution_FBA = optimizeCbModel(model);
fprintf('BIOMASA  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EX_CO2  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso='VALIDACION CITRATO: (bm emp 0.5)';

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0;


model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

model.lb(findRxnIDs(model,'EX_sucr(e)'))=0; 
model.ub(findRxnIDs(model,'EX_sucr(e)'))=0;

model.lb(findRxnIDs(model,'EX_cit(e)'))=-13.45;
model.ub(findRxnIDs(model,'EX_cit(e)'))=-13.45;

model.lb(findRxnIDs(model,'EX_glyc(e)'))=0;
model.ub(findRxnIDs(model,'EX_glyc(e)'))=0;

model.lb(findRxnIDs(model,'EX_glcn(e)'))=0;
model.ub(findRxnIDs(model,'EX_glcn(e)'))=0;

model.lb(findRxnIDs(model,'EX_succ(e)'))=0;
model.ub(findRxnIDs(model,'EX_succ(e)'))=0;

solution_FBA = optimizeCbModel(model);
fprintf('BIOMASA  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EX_CO2  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso='VALIDACION GLUCONATO (bm emp 0.46)';

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0;


model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

model.lb(findRxnIDs(model,'EX_sucr(e)'))=0; 
model.ub(findRxnIDs(model,'EX_sucr(e)'))=0;

model.lb(findRxnIDs(model,'EX_cit(e)'))=0;
model.ub(findRxnIDs(model,'EX_cit(e)'))=0;

model.lb(findRxnIDs(model,'EX_glyc(e)'))=0;
model.ub(findRxnIDs(model,'EX_glyc(e)'))=0;

model.lb(findRxnIDs(model,'EX_glcn(e)'))=-7.67;
model.ub(findRxnIDs(model,'EX_glcn(e)'))=-7.67;

model.lb(findRxnIDs(model,'EX_succ(e)'))=0;
model.ub(findRxnIDs(model,'EX_succ(e)'))=0;

solution_FBA = optimizeCbModel(model);
fprintf('BIOMASA  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EX_CO2  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso='VALIDACION GLICEROL (bm emp 0.18)';

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0;

model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

model.lb(findRxnIDs(model,'EX_sucr(e)'))=0; 
model.ub(findRxnIDs(model,'EX_sucr(e)'))=0;

model.lb(findRxnIDs(model,'EX_cit(e)'))=0;
model.ub(findRxnIDs(model,'EX_cit(e)'))=0;

model.lb(findRxnIDs(model,'EX_glyc(e)'))=-4.39;
model.ub(findRxnIDs(model,'EX_glyc(e)'))=-4.39;

model.lb(findRxnIDs(model,'EX_glcn(e)'))=0;
model.ub(findRxnIDs(model,'EX_glcn(e)'))=0;

model.lb(findRxnIDs(model,'EX_succ(e)'))=0;
model.ub(findRxnIDs(model,'EX_succ(e)'))=0;

solution_FBA = optimizeCbModel(model);
fprintf('BIOMASA  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EX_CO2  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso='VALIDACION SUCCINATO (bm emp 0.38)';

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0;

model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

model.lb(findRxnIDs(model,'EX_sucr(e)'))=0; 
model.ub(findRxnIDs(model,'EX_sucr(e)'))=0;

model.lb(findRxnIDs(model,'EX_cit(e)'))=0;
model.ub(findRxnIDs(model,'EX_cit(e)'))=0;

model.lb(findRxnIDs(model,'EX_glyc(e)'))=0;
model.ub(findRxnIDs(model,'EX_glyc(e)'))=0;

model.lb(findRxnIDs(model,'EX_glcn(e)'))=0;
model.ub(findRxnIDs(model,'EX_glcn(e)'))=0;

model.lb(findRxnIDs(model,'EX_succ(e)'))=-14.07;
model.ub(findRxnIDs(model,'EX_succ(e)'))=-14.07;

solution_FBA = optimizeCbModel(model);
fprintf('BIOMASA  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EX_CO2  %s:%f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));
