%%FBA,GIMME version 11.08.19 NP2019

%%Este codigo crea una tabla comparativa de flujos para 
%%FBA y GIMME en formato Excel

%%Tiene como entrada un modelo de reconstruccion metabolica en cobra
%%(formato mat) y un vector de proteomica para cada caso a analizar (f1,f2
%%y f4)

%%se recomienda usar gurobi como solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

%%comandos para importar y exportar un modelo en excel
%model = xls2model('ModeloNP.xls')
%writeCbModel(model,'xls','out')

%aca se puede cambiar el solver
%changeCobraSolver('glpk','all');
%changeCobraSolver('gurobi','all');

%nombre del archivo de salida con los flujos
archivo='res2.xlsx';

%modelo de entrada
load('modelo_1481c.mat') %este nombre debe coincidir con el modelo usado en genexp.m

%calcula vector triestado de expresion aplicado al modelo automaticamente
genexp 

threshold_lb=-0.5;
threshold_ub=0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%[OPCIONES DE PREPROCESADO DEL MODELO

%%aca se pueden restringir reacciones
%model.lb(findRxnIDs(model,'FCLT'))=1;
%model.ub(findRxnIDs(model,'FCLT'))=1; 

%model.lb(findRxnIDs(model,'FERO'))=0.1;
%model.ub(findRxnIDs(model,'FERO'))=0.1;

%aca se pueden agregar reacciones
%model = addReaction(model, 'ACONTex2','metaboliteList',{'acon_C[p]','acon_C[c]'},'stoichCoeffList',[-1,1])
%model.lb(findRxnIDs(model,'ACONTex2'))=-1000;
%model.ub(findRxnIDs(model,'ACONTex2'))=1000; 

%sucrose alpha glucohidrolase
%model = addReaction(model, 'SAGH','metaboliteList',{'sucr[c]','h2o[c]','f1p[c]','glc_D[c]'},'stoichCoeffList',[-1,-1,1,1])
%model.lb(findRxnIDs(model,'SAGH'))=-1000;
%model.ub(findRxnIDs(model,'SAGH'))=1000; 

%model = addReaction(model, 'SUCRtex','metaboliteList',{'sucr[e]','sucr[p]'},'stoichCoeffList',[-1,1])
%model.lb(findRxnIDs(model,'SUCRtex'))=-1000;
%model.ub(findRxnIDs(model,'SUCRtex'))=1000; 

%model = addReaction(model, 'SUCRtpp','metaboliteList',{'sucr[p]','sucr[c]'},'stoichCoeffList',[-1,1])
%model.lb(findRxnIDs(model,'SUCRtpp'))=-1000;
%model.ub(findRxnIDs(model,'SUCRtpp'))=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%aca se pueden modificar gprs
%model = changeGeneAssociation(model, 'ACONT', 'PSPTO_2016');

%aca se pueden quitar reacciones
model = removeRxns(model, {'ACONT'}); %SUCDi,ACONT,FBA2
model = removeRxns(model, {'FBA2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE BIOMASA (VER EXCEL)
%%biomasa modificada sin hierro
%model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45,-38.5657,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.0001,-0.2323,-0.2288,-0.0001,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.0001,-0.0833,-0.0786,-0.0737,-0.051,-0.0001,-0.0492,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-2.2674,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,45,45,45,1.1372]);																		

%biomasa original sin hierro
%model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											

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


hoja='DC3000';
xlswrite(archivo,[{'Reaction','FBA_F1','GIMME_F1','FBA_F2','GIMME_F2','FBA_F4','GIMME_F4'}],hoja,'A1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE POR CASO

%%%CASO CARBONO LIMITADO (F1)

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.96; 
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;

model.lb(findRxnIDs(model,'EX_co2(e)'))=3.06;
model.ub(findRxnIDs(model,'EX_co2(e)'))=3.06;

model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

expData.value=gen_expr.f1;
  
[gene_names, gene_expr] = findUsedGenesLevels(model, expData);
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expData);

rxn_exp.f1=expressionRxns;
rxn_exp.f1=selectGeneFromGPR(model, gene_names, gene_expr, parsedGPR);
modelo_GIMME=GIMME(model, rxn_exp.f1,0.5);
  
solution_FBA = optimizeCbModel(model);
solution_GIMME = optimizeCbModel(modelo_GIMME);

res.FBA_rxn=model.rxns;
res.FBA_sol=solution_FBA.full;

res.GIMME_rxn=modelo_GIMME.rxns;
res.GIMME_sol=solution_GIMME.full;


xlswrite(archivo,res.FBA_rxn,hoja,'A2');
xlswrite(archivo,res.FBA_sol,hoja,'B2');

%%FBA vs GIMME
tmp=zeros(numel(res.FBA_rxn),1);
for i=1:numel(res.FBA_rxn)
  for j=1:numel(res.GIMME_rxn)
    if strcmpi(res.FBA_rxn{i,1},res.GIMME_rxn{j,1})==1
       %fprintf('%s:%i %i %s\n',res.FBA_rxn{i,1},i,j,res.GIMME_sol(j));
       tmp(i)=res.GIMME_sol(j);
    end
  end  
end  
xlswrite(archivo,tmp,hoja,'C2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k=0
%for i=1:numel(res.FBA_rxn)
%    if res.GIMME_sol(i)>2*res.FBA_sol(i) && res.GIMME_sol(i)>0.001
%        fprintf('%d %s\t%i\t%i\tSOBREEXP\n',k,char(res.FBA_rxn(i)),res.FBA_sol(i),res.GIMME_sol(i));
%        k=k+1;
%    end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caso='f1';
fprintf('BIOMASA F1   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'BIOMASS_F1')));
fprintf('EX_glc_D(e) %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_glc_D(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_glc_D(e)')));
fprintf('EX_co2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_co2(e)')));
fprintf('EX_fe2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_fe2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_fe2(e)')));
fprintf('EX_cobalt2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_cobalt2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_cobalt2(e)')));


%%%CASO ELEMENTOS TRAZA LIMITADO (F2)

model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-1.64;
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-1.64; 
  
model.lb(findRxnIDs(model,'EX_co2(e)'))=7.36;
model.ub(findRxnIDs(model,'EX_co2(e)'))=7.36;
  
model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0; 
  
model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=0; 
model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0;  
  
expData.value=gen_expr.f2;
  
[gene_names, gene_expr] = findUsedGenesLevels(model, expData);
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expData)

rxn_exp.f2=expressionRxns;
rxn_exp.f2=selectGeneFromGPR(model, gene_names, gene_expr, parsedGPR);
modelo_GIMME=GIMME(model, rxn_exp.f2,0.5);
  
solution_FBA = optimizeCbModel(model);
solution_GIMME = optimizeCbModel(modelo_GIMME);

res.FBA_rxn=model.rxns;
res.FBA_sol=solution_FBA.full;

res.GIMME_rxn=modelo_GIMME.rxns;
res.GIMME_sol=solution_GIMME.full;

xlswrite(archivo,res.FBA_sol,hoja,'D2');

%%FBA vs GIMME
tmp=zeros(numel(res.FBA_rxn),1);
for i=1:numel(res.FBA_rxn)
  for j=1:numel(res.GIMME_rxn)
    if strcmpi(res.FBA_rxn{i,1},res.GIMME_rxn{j,1})==1
       %fprintf('%s:%i %i %s\n',res.FBA_rxn{i,1},i,j,res.GIMME_sol(j));
       tmp(i)=res.GIMME_sol(j);
    end
  end  
end  
xlswrite(archivo,tmp,hoja,'E2');
  
caso='f2';
fprintf('BIOMASA F2  %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'BIOMASS_F1')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'BIOMASS_F1')));
fprintf('EX_glc_D(e) %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_glc_D(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_glc_D(e)')));
fprintf('EX_co2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_co2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_co2(e)')));
fprintf('EX_fe2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_fe2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_fe2(e)')));
fprintf('EX_cobalt2(e)   %s:%f %f \n',caso,solution_FBA.x(findRxnIDs(model,'EX_cobalt2(e)')),solution_GIMME.full(findRxnIDs(modelo_GIMME,'EX_cobalt2(e)')));

%%%CASO HIERRO LIMITADO (F4)  
  
model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.66;
model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.66; 
  
model.lb(findRxnIDs(model,'EX_co2(e)'))=0.65;
model.ub(findRxnIDs(model,'EX_co2(e)'))=0.65;
  
model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

expData.value=gen_expr.f4;
  
[gene_names, gene_expr] = findUsedGenesLevels(model, expData);
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expData)

rxn_exp.f4=expressionRxns;
rxn_exp.f4=selectGeneFromGPR(model, gene_names, gene_expr, parsedGPR);
modelo_GIMME=GIMME(model, rxn_exp.f4,0.5);
  
solution_FBA = optimizeCbModel(model);
solution_GIMME = optimizeCbModel(modelo_GIMME);

res.FBA_rxn=model.rxns;
res.FBA_sol=solution_FBA.full;

res.GIMME_rxn=modelo_GIMME.rxns;
res.GIMME_sol=solution_GIMME.full;

xlswrite(archivo,res.FBA_sol,hoja,'F2');

%%FBA vs GIMME
tmp=zeros(numel(res.FBA_rxn),1);
for i=1:numel(res.FBA_rxn)
  for j=1:numel(res.GIMME_rxn)
    if strcmpi(res.FBA_rxn{i,1},res.GIMME_rxn{j,1})==1
       %fprintf('%s:%i %i %s\n',res.FBA_rxn{i,1},i,j,res.GIMME_sol(j));
       tmp(i)=res.GIMME_sol(j);
    end
  end  
end  
xlswrite(archivo,tmp,hoja,'G2');

%%%%%%[MOSTRAR REACCIONES
%for i=1:numel(model.rxns)
%    fprintf('%i: ',i)
%    printRxnFormula(model, 'rxnAbbrList', model.rxns{i});
%end


%%printFluxVector(model, solution_FBA.x)

%al finalizar mostrar un resumen de estas reaccioens
rxn1=[{'BIOMASS_F1','EX_glc_D(e)','EX_co2(e)','EX_fe2(e)','EX_cobalt2(e)'}];
rxn2=[{'EDA','EDD','CS','ACONT','ADK1','HEX1','6PG'}];

for i=1:numel(rxn1)
  fprintf('%s:\t%f %f \n',char(rxn1(i)),solution_FBA.x(findRxnIDs(model,rxn1(i))),solution_GIMME.full(findRxnIDs(modelo_GIMME,rxn1(i))));
end


