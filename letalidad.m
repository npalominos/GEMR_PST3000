%%ANALISIS DE LETALIDAD

load('nestor101018.mat')
solution = optimizeCbModel(model)
k=1;

A=[]
echo off;
for i=-3.0:0.02:0.0

  model.lb(findRxnIDs(model,'EX_glc_D(e)'))=i %%-1.04
  model.ub(findRxnIDs(model,'EX_glc_D(e)'))=0
  solution = optimizeCbModel(model);

  disp('***');
  disp(solution.obj);
  disp('***');

  A=[A;i,solution.obj]

  k=k+1;
end
echo on;

for i=1:1:size(A);
  disp(A')
end;

f=figure;
p = plot(A(:,1),A(:,2));
p(1).LineWidth = 2;
p(2).Marker = '.';
ax = gca
%ax.XLim = [0 0.55];
%ax.YLim = [0 0.2];