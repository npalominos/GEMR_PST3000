%%Robustness Analysis + proteomica NP2019

%%Este codigo obtiene un grafico del cambio de flujo en una reaccion objetivo
%%a partir del cambio de flujo de una reaccion de control manteniendo 
%%el resto del sistema constante, ademas de trabajar en el modelo con proteomica integrada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

%aca se puede cambiar el solver
%changeCobraSolver('glpk','all');
%changeCobraSolver('gurobi','all');

%reacciones de control y objetivo
controlRxn='EX_fe2(e)';
objRxn='EX_fe3(e)';

%controlRxn='EX_glc_D(e)';
%objRxn='BIOMASS_F1';

archivo='out.xlsx';
hoja='sensibilidad';

xlswrite(archivo,[{'controlRxn_FBA','objRxn_FBA','controlRxn_GIMME','objRxn_GIMME'}],hoja,'A1');

%modelo de entrada
load('modelo_1481c.mat') %este nombre debe coincidir con el modelo usado en genexp.m

%calcula vector triestado de expresion aplicado al modelo automaticamente
genexp 

threshold_lb=-0.5;
threshold_ub=0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%[OPCIONES DE PREPROCESADO DEL MODELO
%aca se pueden quitar reacciones
model = removeRxns(model, {'SUCDi'});
%model = removeRxns(model, {'FBA2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE BIOMASA (VER EXCEL)
%biomasa original con hierro
model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','hemeO[c]','sheme[c]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0005,-0.0005,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											
printRxnFormula(model, 'BIOMASS_F1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE GENERALES
model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=-1000;
model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0; 

model.lb(findRxnIDs(model,'EX_nh4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_nh4(e)'))=0; 

model.lb(findRxnIDs(model,'EX_o2(e)'))=-30;
model.ub(findRxnIDs(model,'EX_o2(e)'))=0;

model.lb(findRxnIDs(model,'EX_so4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_so4(e)'))=0;

parsedGPR=GPRparser(model);
expData.gene=model.genes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE POR CASO

%%%CASO HIERRO LIMITADO (F4)  
  
model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.66;
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.66; 
  
model.lb(findRxnIDs(model,'EX_co2(e)'))=0.65;
model.ub(findRxnIDs(model,'EX_co2(e)'))=0.65;
  
model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[INTEGRACION PROTEOMICA
expData.value=gen_expr.f4;
  
[gene_names, gene_expr] = findUsedGenesLevels(model, expData);
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expData)

rxn_exp.f4=expressionRxns;
rxn_exp.f4=selectGeneFromGPR(model, gene_names, gene_expr, parsedGPR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[DEFINICION DE REACCION DE CONTROL Y REACCION OBJETIVO
n=50;

x = zeros(n+1, 1);
y = zeros(n+1, 1);
w = zeros(n+1, 1);
z = zeros(n+1, 1);

k=1;

xmin=-1.0;
xmax=1.0;
paso=(xmax-xmin)/n;

for i = xmin:paso:xmax
  %fprintf('%d:%f',k,i)
  model.lb(findRxnIDs(model,controlRxn))=i; 
  model.ub(findRxnIDs(model,controlRxn))=i;

  try
   modelo_GIMME=GIMME(model, rxn_exp.f4,0.5);
   solution_FBA = optimizeCbModel(model);
   solution_GIMME = optimizeCbModel(modelo_GIMME);

   x(k)=solution_GIMME.full(findRxnIDs(modelo_GIMME,controlRxn));
   y(k)=solution_GIMME.full(findRxnIDs(modelo_GIMME,objRxn));

   w(k)=solution_FBA.full(findRxnIDs(model,controlRxn));
   z(k)=solution_FBA.full(findRxnIDs(model,objRxn));
  catch
   warning('excepcion!');
  end
  k=k+1;
   
end

xlswrite(archivo,x,hoja,'A2');
xlswrite(archivo,y,hoja,'B2');
xlswrite(archivo,w,hoja,'C2');
xlswrite(archivo,z,hoja,'D2');
  
plot (xmin:paso:xmax, y);
plot (xmin:paso:xmax, z);
xlabel(controlRxn);
ylabel(objRxn);


