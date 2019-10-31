%%ROBUSTNESS ANALYSIS
%%NP2019

%%En el modelo utilizado, se agrego el prefijo 'EX_' a todas las reacciones
%%Exchange. A dicho modelo, se le nombro 'model_small2.mat'


%changeCobraSolver ('gurobi', 'all', 1);
%changeCobraSolver ('glpk', 'all', 1);

load('model_small2.mat');

%changeObjective(model,'C10-PHA', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[RESTRICCIONES GENERALES

%%De acuerdo a cobratoolbox, los exchange reversibles
%%deben ser transformados a reversibles, por lo que el lower bound 
%%de estos debe ser fijado en cero

model.lb(findRxnIDs(model,'EX_NH3(c)'))=-100;
model.ub(findRxnIDs(model,'EX_NH3(c)'))=0;

model.lb(findRxnIDs(model,'EX_SO4(c)'))=-100;
model.ub(findRxnIDs(model,'EX_SO4(c)'))=0;

model.lb(findRxnIDs(model,'EX_O2(c)'))=-100;
model.ub(findRxnIDs(model,'EX_O2(c)'))=0;

model.lb(findRxnIDs(model,'EX_Rxn5Biomass'))=-100;
model.ub(findRxnIDs(model,'EX_Rxn5Biomass'))=0;

model.lb(findRxnIDs(model,'EX_CAT_ex'))=0.02;
model.ub(findRxnIDs(model,'EX_CAT_ex'))=0.02;

model.lb(findRxnIDs(model,'EX_BEN(e)'))=-2.26;
model.ub(findRxnIDs(model,'EX_BEN(e)'))=-2.26;

%%%

model.lb(findRxnIDs(model,'EX_CO2(c)'))=0;
model.ub(findRxnIDs(model,'EX_CO2(c)'))=100;

model.lb(findRxnIDs(model,'EX_MUC_ex(e)'))=0.688;
model.ub(findRxnIDs(model,'EX_MUC_ex(e)'))=0.688;

model.lb(findRxnIDs(model,'catA'))=-5;
model.ub(findRxnIDs(model,'catA'))=5;

model.lb(findRxnIDs(model,'C10-PHA'))=0.26;
model.ub(findRxnIDs(model,'C10-PHA'))=0.26;
%%%%%%%%%

for i=1:length(model.rxns)
  if strncmp(model.rxns{i},'EX_',3)
    fprintf('%s\n',model.rxns{i});
    model.subSystems{i}='Exchange/demand reaction';
  else
    model.subSystems{i}=''
  end
end

idx=strmatch('Exchange/demand reaction', model.subSystems);

c=0;

for i=1:length(idx)
 if model.lb(idx(i))~=0
   c=c+1;
   uptakes{c}=model.rxns{idx(i)};
 end
end

modelalter = model;
%modelalter = changeRxnBounds(modelalter, uptakes, 0, 'l');
modelrobust = modelalter;

%%ACA SE PUEDE CAMBIAR LA VARIABLE DE CONTROL Y EL OBJETIVO
controlRxn='catA';
objRxn='C10-PHA';

x = zeros(401, 1);
k=0;
for i = -1.0:0.01:3.0
  %fprintf('%d:%f',k,i)
  
  modelrobust = changeRxnBounds(modelrobust, controlRxn, -i, 'b');
  modelrobust = changeObjective(modelrobust, objRxn);
  FBArobust = optimizeCbModel(modelrobust, 'max');
  x(k+1) = FBArobust.f;
  k=k+1;
end

plot (-1.0:0.01:3.0, x);
xlabel(controlRxn);
ylabel(objRxn);

%npoints=20;
%plotResFlag=true
%objType='max'
%[controlFlux, objFlux] = robustnessAnalysis(model, controlRxn, npoints,plotResFlag, objRxn,objType)

