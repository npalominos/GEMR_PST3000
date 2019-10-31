%%SISTEMA DE APAGADO DE REACCIONES
%%NESTOR PALOMINOS 2018

echo off;

%elimina de a una reaccion desde la 1040 en adelante,puesto a que el modelo
%original es de 1039 reacciones

for i = 1822:-1:1039

  load('modelo_1822.mat')
  preprocesamiento
  
  rxn=model.rxns{i};
  
  try
   model = removeRxns(model, {rxn});

   changeObjective(model,'BIOMASS_F1', 1);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE POR CASO

   %%%CASO CARBONO LIMITADO (F1)

   model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;
   model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.96;

   model.lb(findRxnIDs(model,'EX_co2(e)'))=3.06;
   model.ub(findRxnIDs(model,'EX_co2(e)'))=3.06;

   model.lb(findRxnIDs(model,'EX_fe2(e)'))=-30; 
   model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;   

   res.FBA_rxn=model.rxns;  

   solution_FBA_F1 = optimizeCbModel(model);
   res.FBA_sol=solution_FBA_F1.full;

   %%%CASO ELEMENTOS TRAZA LIMITADO (F2)

   model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-1.64;
   model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-1.64; 
  
   model.lb(findRxnIDs(model,'EX_co2(e)'))=7.36;
   model.ub(findRxnIDs(model,'EX_co2(e)'))=7.36;
  
   model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
   model.ub(findRxnIDs(model,'EX_fe2(e)'))=0; 
  
   model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=0; 
   model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0;  
  
   solution_FBA_F2 = optimizeCbModel(model);
   res.FBA_sol=solution_FBA_F2.full;

   %%%CASO HIERRO LIMITADO (F4)  
  
   model.lb(findRxnIDs(model,'EX_glc_D(e)'))=-0.66;
   model.ub(findRxnIDs(model,'EX_glc_D(e)'))=-0.66; 
  
   model.lb(findRxnIDs(model,'EX_co2(e)'))=0.65;
   model.ub(findRxnIDs(model,'EX_co2(e)'))=0.65;
  
   model.lb(findRxnIDs(model,'EX_fe2(e)'))=0; 
   model.ub(findRxnIDs(model,'EX_fe2(e)'))=0;

   solution_FBA_F4 = optimizeCbModel(model);
   res.FBA_sol=solution_FBA_F4.full;

   
   %al finalizar mostrar un resumen de estas reacciones en un excel
   rxn1=[{'BIOMASS_F1','HEX1','GNK','G6PDH2','EDA','EDD','CS','ACONT','ADK1','EX_glc_D(e)','EX_co2(e)','EX_fe2(e)','EX_cobalt2(e)'}];

   file = fopen('apagado.csv','a+t');
   fprintf('\n');
   n=numel(rxn1);
   fprintf(file,'\n%d;%s;',i,rxn);

   for i=1:numel(rxn1)
    x1=solution_FBA_F1.x(findRxnIDs(model,rxn1(i)));
    x2=solution_FBA_F2.x(findRxnIDs(model,rxn1(i)));
    x3=solution_FBA_F4.x(findRxnIDs(model,rxn1(i)));
    fprintf('%s:\t%f %f %f\n',char(rxn1(i)),x1,x2,x3);
    fprintf(file,'%s:\t%f %f %f;',char(rxn1(i)),x1,x2,x3);
   end

   fclose(file);
  catch
   warning('excepcion!');
  end
  
  
end


  