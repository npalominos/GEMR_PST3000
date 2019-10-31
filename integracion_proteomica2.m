%%FBA,iNEP version 11.08.19 NP2019

%%Este codigo crea una tabla comparativa de flujos para 
%%FBA e iNEP en formato Excel

%%el algoritmo iNEP, obtiene aquellas reacciones asociadas a un unico gen
%%en su GPR y luego limita su flujo a 0, para posteriormente realizar un
%%FBA. 

%%Tiene como entrada un modelo de reconstruccion metabolica en cobra
%%(formato mat) y un vector de proteomica para cada caso a analizar (f1,f2
%%y f4)

%%se recomienda usar gurobi como solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

%changeCobraSolver('glpk','all');
%changeCobraSolver('gurobi','all');

archivo='res.xlsx';
 
load('proteomica_794.mat');
load('gen.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=[{'f1','f2','f4'}];

for z=1:3
    
caso=c{1,z};
load('modelo_1481c.mat');
model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											
%model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45,-38.5657,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.0001,-0.2323,-0.2288,-0.0001,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.0001,-0.0833,-0.0786,-0.0737,-0.051,-0.0001,-0.0492,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-2.2674,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,45,45,45,1.1372]);																		


%%%%%%%%%%%%%%%%[MODIFICACIONES AL MODELO
%%GK1 REVERSIBLE

%model = addReaction(model, 'GK1','metaboliteList',{'atp[c]','gmp[c]','adp[c]','gdp[c]'},'stoichCoeffList',[-1,-1,1,1])											
printRxnFormula(model, 'rxnAbbrList', 'GK1');

%model.grRules(findRxnIDs(model,'GK1'))='';      %PSPTO_0075 (gmk)

%model = changeGeneAssociation(model, 'ICS', '');
%model = changeGeneAssociation(model, 'GK1', '');     %PSPTO_0075 (gmk)
%model = changeGeneAssociation(model, 'QULNS', '');   %PSPTO_3959 (nadA)
%model = changeGeneAssociation(model, 'TMDS', '');    %PSPTO_5282 (thyA)
%model = changeGeneAssociation(model, 'UAMAGS', '');  %PSPTO_4410 (murD)
%model = changeGeneAssociation(model, 'UAMAS', '');   %PSPTO_4407 (murC)
%model = changeGeneAssociation(model, 'AICART', '');  %PSPTO_4866 (purH)
%model = changeGeneAssociation(model, 'APSR', '');    %PSPTO_2280 (cysH)

%model = changeGeneAssociation(model, 'ACONTa', 'PSPTO_2016');    
%model = changeGeneAssociation(model, 'ACONTb', 'PSPTO_2016');  
%model = changeGeneAssociation(model, 'PPC', 'PSPTO_1508');

%model.lb(findRxnIDs(model,'ACONTa'))=-100; 
%model.ub(findRxnIDs(model,'ACONTa'))=100;

%model.lb(findRxnIDs(model,'ACONTb'))=-100; 
%model.ub(findRxnIDs(model,'ACONTb'))=100;

%model.lb(findRxnIDs(model,'AKGDH'))=-100; 
%model.ub(findRxnIDs(model,'AKGDH'))=100;

model.lb(findRxnIDs(model,'PC'))=0; 
model.ub(findRxnIDs(model,'PC'))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE GENERALES
model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=-1000;
model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0; 

model.lb(findRxnIDs(model,'EX_nh4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_nh4(e)'))=0; 

model.lb(findRxnIDs(model,'EX_o2(e)'))=-30;
model.ub(findRxnIDs(model,'EX_o2(e)'))=0;

model.lb(findRxnIDs(model,'EX_so4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_so4(e)'))=0;


if caso=='f1'
  %%%CASO CARBONO LIMITADO (F1)
  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;

  model.lb(findRxnIDs(model,'EX_co2(e)'))=3.06;
  model.ub(findRxnIDs(model,'EX_co2(e)'))=3.06;
end

if caso=='f2'
  %%%CASO ELEMENTOS TRAZA LIMITADO (F2)

  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-1.64;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-1.64; 
  
  model.lb(findRxnIDs(model,'EX_co2(e)'))=7.36;
  model.ub(findRxnIDs(model,'EX_co2(e)'))=7.36;
  
  model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_fe2(e)'))=0; 
  
  model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0;
end

if caso=='f4'
 %%%CASO HIERRO LIMITADO (F4)  
  
  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.66;
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.66; 
  
  model.lb(findRxnIDs(model,'EX_co2(e)'))=0.65;
  model.ub(findRxnIDs(model,'EX_co2(e)'))=0.65;
  
  model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
  model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SE BUSCA LA EXPRESION ASOCIADA A CADA GEN DEL MODELO SEGUN LOS DATOS
%%EXPERIMENTALES DE PROTEOMICA.

%%orden de las columnas de gen_expr2
%%pspto; name; f1; f2; f4

n=numel(model.genes);      %genes del modelo
m=numel(proteomica.gene);  %genes de proteomica

gen_expr2=cell(n,5);

for i=1:1:n
  gen_expr2{i,1}=model.genes{i};
  
  encontrado=0;
  for j=1:1:m
      if encontrado==1
          encontrado=0;
          break;
      end
          
      if strcmp(model.genes{i,1},proteomica.gene{j,1})
          encontrado=1;
          %fprintf('%s\n', char(model.genes(i)));
          gen_expr2{i,3}=proteomica.f1{j,1};
          gen_expr2{i,4}=proteomica.f2{j,1};
          gen_expr2{i,5}=proteomica.f4{j,1};
      end
          
    end
end
  
%%busca el nombre del gen segun su locus tag 'pspto'
for i=1:numel(model.genes)
  gen_expr2{i,2}={'nn'};
  for k=1:numel(gen.pspto)
             
    if strcmp(model.genes{i,1},gen.pspto{k,1})       
      %fprintf('**%d %s %s %s\n', i,char(model.genes{i,1}),char(gen.pspto(k)), char(gen.name(k)));
      gen_expr2{i,2}=gen.name(k);
      break;
    end
    
  end  %%fin del for j
end    %%fin del for i

%%en este punto ya esta formada la matriz de expresiones.
%%ahora se debe buscar aquellas reacciones asociadas a un unico gen 
%%y limitar su flujo como restriccion anexa al fba

%%calculo la cantidad de genes asociados a la reaccion en su gpr
%donde se extran las reacciones asociadas a un unico gen
rules.id=zeros(n,1);
rules.num=zeros(n,1);

r=0  %contador de reacciones asociadas a un unico gen

%%aca se extraen las reacciones con gpr no vacio y se buscan las de gpr
%%unitario
rxn_unico_gen=cell(1,2);
for i=1:numel(model.rxns)
  x=isempty(char(model.rules(i)));
  y=length(strfind(char(model.rules(i)),'x'));
  rules.id(i)=x;
  rules.num(i)=y;
  
  if y==1
      %%agregar a listado de reacciones asociadas a un unico gen
      fprintf('%d %s %s\n',r,char(model.rxns(i)),char(model.grRules(i)));
      rxn_unico_gen{r+1,1}=char(model.rxns(i));
      rxn_unico_gen{r+1,2}=char(model.grRules(i));
      r=r+1;

  end

end

modelo_orig=model;

%%fba sin modificaciones
fprintf('FBA MODELO ORIGINAL');
solution_FBA = optimizeCbModel(modelo_orig);

%%en este punto ya se tiene el listado de reacciones asociadas a un unico
%%gen en su gpr
%%ahora se extraera la expresion asociada a estos genes y se agregara a la
%%tabla rxn_unico_gen
rxn_letales=[];
p=0;
for i=1:size(rxn_unico_gen,1)
  for j=1:size(gen_expr2,1)
    if strcmp(gen_expr2{j,1},rxn_unico_gen{i,2})  
      nam=char(gen_expr2{j,2});
      ge1=gen_expr2{j,3};
      ge2=gen_expr2{j,4};
      ge4=gen_expr2{j,5};
      
      rxn_unico_gen{i,3}=ge1;
      rxn_unico_gen{i,4}=ge2;
      rxn_unico_gen{i,5}=ge4;
      
       %%para cada reaccion asociada a un gen, la reprimo limitando el flujo
      
       %%model.lb(findRxnIDs(model,rxn_unico_gen{i,1}))=0;
       %%model.ub(findRxnIDs(model,rxn_unico_gen{i,1}))=0;
       
       model2=model;
       %if strcmp(rxn_unico_gen{i,1},'GADT')==0 && strcmp(rxn_unico_gen{i,1},'VALTA')==0  && strcmp(rxn_unico_gen{i,1},'UAMAGS')==0  && strcmp(rxn_unico_gen{i,1},'UAPGR')==0  && strcmp(rxn_unico_gen{i,1},'UAMAS')==0
       
       if caso=='f1' 
         if ge1==0
          model = removeRxns(model, {rxn_unico_gen{i,1}});
          fprintf('F1 se elimino: %s',rxn_unico_gen{i,1});
         end
       end

       if caso=='f2'
         if ge2==0
          model = removeRxns(model, {rxn_unico_gen{i,1}});
          fprintf('F2 se elimino: %s',rxn_unico_gen{i,1});
         end
       end
       
       if caso=='f4' 
         if ge4==0
          model = removeRxns(model, {rxn_unico_gen{i,1}});
          fprintf('F4 se elimino: %s',rxn_unico_gen{i,1});
         end
       end
   
       solution_iNEP = optimizeCbModel(model);
       
       
       if solution_iNEP.obj<0.03
           model=model2;
           fprintf('-------> %s',char(rxn_unico_gen{i,1}));

           rxn_letales{p+1,1}=char(rxn_unico_gen{i,1});
           fprintf('killed\n');
           p=p+1;
       end
       fprintf('**%d %s %s %s f1:%d f2:%d f4:%d fba:%f\n', i,char(rxn_unico_gen{i,1}),char(rxn_unico_gen{i,2}),nam,ge1,ge2,ge4,solution_iNEP.obj);

      break;
    end
  end
end

fprintf('REACCIONES ASOCIADAS A UN UNICO GEN: %d\n',r);

%%en este punto, el modelo tiene restringidos los flujos para als
%%reacciones asociadas a un unico gen

%se realizara entonces un nuevo FBA
solution_iNEP = optimizeCbModel(model)

archivo='excel.xls';
hoja='INTEGRACION PROTEOMICA'; %caso;

%%EXPORTANDO FLUJOS A EXCEL
xlswrite(archivo,[{'Reaction','FBA_F1','iNEP_F1','FBA_F2','iNEP_F2','FBA_F4','iNEP_F4'}],hoja,'A1');

xlswrite(archivo,modelo_orig.rxns,hoja,'A2');

if caso=='f1' 
 xlswrite(archivo,solution_FBA.full,hoja,'B2');
end

if caso=='f2' 
 xlswrite(archivo,solution_FBA.full,hoja,'D2');
end

if caso=='f4' 
 xlswrite(archivo,solution_FBA.full,hoja,'F2');
end

tmp=zeros(numel(solution_FBA.full),1);
for i=1:numel(modelo_orig.rxns)
  for j=1:numel(model.rxns)

    if strcmpi(modelo_orig.rxns(i),model.rxns(j))==1
       tmp(i)=solution_iNEP.full(j);
       break;
    end
  end  
end  

if caso=='f1' 
 xlswrite(archivo,tmp,hoja,'C2');
end

if caso=='f2' 
 xlswrite(archivo,tmp,hoja,'E2');
end

if caso=='f4' 
 xlswrite(archivo,tmp,hoja,'G2');
end
%al finalizar mostrar un resumen de estas reaccioens
rxn1=[{'BIOMASS_F1','EX_glc_D(e)','EX_co2(e)','EX_fe2(e)','EX_cobalt2(e)'}];
rxn2=[{'GK1','QULNS','TMDS','UAMAGS','UAMAS','AICART','EDA','EDD','CS','ACONT','ADK1','HEX1','6PG'}];

for i=1:numel(rxn1)
  fprintf('%s:\t%f %f \n',char(rxn1(i)),solution_FBA.x(findRxnIDs(modelo_orig,rxn1(i))),solution_iNEP.full(findRxnIDs(model,rxn1(i))));
end
%for i=1:numel(rxn2)
 
 %fprintf('%s:\t%f \n',char(rxn2(i)),solution_FBA.x(findRxnIDs(modelo_orig,rxn2(i))));

 % try
 %  fprintf('%s:\t%f\n',char(rxn2(i)),solution_iNEP.full(findRxnIDs(model,rxn2(i))));
%  catch
 %  warning('excepcion!');
 % end
%end

end %casos