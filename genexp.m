%%PREPROCESAMIENTO DE PROTEOMICA Y GENERACION DE VECTOR TRIESTADO
%%NESTOR PALOMINOS 2019

%%ENTRADAS:
%%MODELO DE RECONSTRUCCION METABOLICA (rxns,mets,genes)
%%PROTEOMICA SEGUN CONTEXTOS (f1,f2,f4,gene)

%%PARA CADA GEN PERTENECIENTE AL MODELO SE BUSCA EN LA TABLA DE EXPRESIONES
%%SI EXISTE, SE CLASIFICA EN -1,0,1 SEGUN LA SIGUIENTE REGLA
%%-1<E-k*STD<0<E+k*STD<1 (-1:REP 0:MODERADO 1:SOBREEXP)

%%CTE k DEFINIDA POR EL USUARIO (COMUNMENTE USADO ENTRE 0.3 A 0.5)

%%SI NO EXISTE, SE CONSIDERA 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0.5  %cte desv standard
ne=0   %valor no exp

load('modelo_1481c.mat') 
load('proteomica_794.mat')

n=numel(model.genes)      %genes del modelo
m=numel(proteomica.gene)  %genes de proteomica

gen_expr.f1=zeros(n,1);
gen_expr.f2=zeros(n,1);
gen_expr.f4=zeros(n,1);


for i=1:1:n
  encontrado=0;
  for j=1:1:m
      
      if encontrado==1
          encontrado=0;
          break;
      end
          
      if strcmp(model.genes{i,1},proteomica.gene{j,1})
          encontrado=1;
          fprintf('%s\n', char(model.genes(i)));
          
          f1=proteomica.f1{j,1};
          f2=proteomica.f2{j,1};
          f4=proteomica.f4{j,1};
          F=[f1;f2;f4];
          
          E=mean(F);
          S=std(F);
          
          u1=E-k*S;
          u2=E+k*S;
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(f1<u1)
              gen_exp.f1(i)=-1;
          elseif(f1>u2)
              gen_exp.f1(i)=1;
          else
              gen_exp.f1(i)=0;
          end
          
          if(f2<u1)
              gen_exp.f2(i)=-1;
          elseif(f2>u2)
              gen_exp.f2(i)=1;
          else
              gen_exp.f2(i)=0;
          end
          
          if(f4<u1)
              gen_exp.f4(i)=-1;
          elseif(f4>u2)
              gen_exp.f4(i)=1;
          else
              gen_exp.f4(i)=0;
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %else
       %   gen_exp.f1=[gen_exp.f1,ne];
       %   gen_exp.f2=[gen_exp.f2,ne];
       %   gen_exp.f4=[gen_exp.f4,ne]; 
      end
      
      
  end
end
    
    
    






