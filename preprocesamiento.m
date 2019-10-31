%%PREPROCESAMIENTO DE BIOMASA Y EXCHANGE GENERALES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[MODIFICACION DE BIOMASA (VER EXCEL)
%%biomasa modificada sin hierro
%model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45,-38.5657,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.0001,-0.2323,-0.2288,-0.0001,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.0001,-0.0833,-0.0786,-0.0737,-0.051,-0.0001,-0.0492,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-2.2674,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,45,45,45,1.1372]);																		

%biomasa original sin hierro
model = addReaction(model, 'BIOMASS_F1','metaboliteList',{'atp[c]','h2o[c]','glu_L[c]','ala_L[c]','gly[c]','asn_L[c]','leu_L[c]','val_L[c]','ile_L[c]','gtp[c]','lys_L[c]','gln_L[c]','ctp[c]','arg_L[c]','asp_L[c]','thr_L[c]','ser_L[c]','phe_L[c]','pro_L[c]','utp[c]','tyr_L[c]','met_L[c]','his_L[c]','trp_L[c]','5mthf[c]','cys_L[c]','pe161[p]','ptrc[c]','peptido_kt[c]','dgtp[c]','pg161[c]','dctp[c]','pe181[p]','dttp[c]','datp[c]','pg181[c]','pe160[p]','pg160[c]','clpn161[p]','udpg[c]','pe180[p]','nad[c]','pg180[c]','clpn181[p]','amp[c]','clpn160[p]','nadph[c]','clpn180[p]','nadp[c]','accoa[c]','nadh[c]','fad[c]','coa[c]','succoa[c]','pi[c]','adp[c]','h[c]','ppi[c]'},'stoichCoeffList',[-45.7318,-45.5608,-0.5516,-0.4595,-0.4178,-0.3827,-0.3524,-0.2997,-0.2526,-0.235,-0.2323,-0.2288,-0.214,-0.2117,-0.2096,-0.189,-0.1792,-0.176,-0.1566,-0.136,-0.0833,-0.0786,-0.0737,-0.051,-0.05,-0.0492,-0.0419,-0.035,-0.028,-0.0232,-0.0231,-0.0231,-0.019,-0.0148,-0.0143,-0.0105,-0.0086,-0.0046,-0.0031,-0.003,-0.0028,-0.0022,-0.0015,-0.0014,-0.001,-0.0006,-0.0004,-0.0002,-0.0001,-0.0001,-0.0001,0,0,0,45.5628,45.5608,45.5604,0.7302])											

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[LIMITACIONES EXCHANGE GENERALES
model.lb(findRxnIDs(model,'EX_cobalt2(e)'))=-1000;
model.ub(findRxnIDs(model,'EX_cobalt2(e)'))=0; 

model.lb(findRxnIDs(model,'EX_nh4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_nh4(e)'))=0; 

model.lb(findRxnIDs(model,'EX_o2(e)'))=-30;
model.ub(findRxnIDs(model,'EX_o2(e)'))=0;

model.lb(findRxnIDs(model,'EX_so4(e)'))=-30;
model.ub(findRxnIDs(model,'EX_so4(e)'))=0;
