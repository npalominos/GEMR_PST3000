%%Analisis para distintas fuentes de carbon version 20.05.19 NP2019

%%Este codigo compara la biomasa resultante segun las siguientes  
%%fuentes de carbon:

%%VALIDACION DEL SUCCINATO

%%SUCCINATO  BIOMASA
%%-1.46     1.387
%%-1.67     0.1586
%%-1.42     0.1349
%%-1,25     0.1187

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

%changeCobraSolver('glpk','all');
%changeCobraSolver('gurobi','all');

load('modelo_1408.mat')

caso='f1'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE BIOMASA (VER EXCEL)
%%biomasa modificada sin hierro
%model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45,-38.5657,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.0001,-0.2323,-0.2288,-0.0001,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.0001,-0.0833,-0.0786,-0.0737,-0.051,-0.0001,-0.0492,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-2.2674,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,45,45,45,1.1372]);																		

%biomasa original sin hierro
model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE GENERALES


model.lb(findRxnIDs(model,'EX_nh4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_nh4(e)'))=0; 

model.lb(findRxnIDs(model,'EX_o2(e)'))=-30;
model.ub(findRxnIDs(model,'EX_o2(e)'))=0;

model.lb(findRxnIDs(model,'EX_so4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_so4(e)'))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE POR CASO

if caso=='f1'
  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;
  
  %model.lb(findRxnIDs(model,'EX_co2(e)'))=3.06;
  %model.ub(findRxnIDs(model,'EX_co2(e)'))=3.06;
  
  model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
  model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;   
  
  model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=-1000;
  model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0; 

elseif caso=='f2'
  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-1.64;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-1.64; 
  
  %model.lb(findRxnIDs(model,'EX_co2(e)'))=7.36;
  %model.ub(findRxnIDs(model,'EX_co2(e)'))=7.36;
  
  model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;   
  
  model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0;  

elseif caso=='f4'
  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.66;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.66; 
  
  %model.lb(findRxnIDs(model,'EX_co2(e)'))=0.65;
  %model.ub(findRxnIDs(model,'EX_co2(e)'))=0.65;
  
  model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;
  
  model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=-1000;
  model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE LA FUENTE DE CARBONO
model.lb(findRxnIDs(model,'EX_glc_D(e)'))=0;
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0; 
  
model.lb(findRxnIDs(model,'EX_succ(e)'))=-1.46;
model.ub(findRxnIDs(model,'EX_succ(e)'))=-1.46; 

model.lb(findRxnIDs(model,'EX_mal_L(e)'))=0;
model.ub(findRxnIDs(model,'EX_mal_L(e)'))=0; 

model.lb(findRxnIDs(model,'EX_co2(e)'))=0;
model.ub(findRxnIDs(model,'EX_co2(e)'))=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solution_FBA = optimizeCbModel(model);

%%printFluxVector(model, solution_FBA.x)

fprintf('->BIOMASA     %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')));
fprintf('EDA         %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EDA')));
fprintf('EDD         %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EDD')));

fprintf('EX_fe2(e)   %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EX_fe2(e)')));
fprintf('EX_co2(e)   %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')));
fprintf('EX_glc_D(e) %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EX_glc_D(e)')));
fprintf('->EX_succ(e)   %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EX_succ(e)')));
fprintf('EX_mal_L(e)   %s:%f\n',caso,solution_FBA.x(findRxnIDs(model,'EX_mal_L(e)')));
