function [biomass]= campaniensis_after_gapfill()
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs\gapfill_results');

model = readCbModel([baseFolder filesep 'campaniensis' filesep 'modelGapfilled_1_v3.mat']);

% revision del gapfilling de boliviensis dice que hay que eliminar rxn
% PPGO3
model = removeRxns(model,{'PPGO3'});

% vemos si genera biomasa
sol = optimizeCbModel(model);
% valor de objetivo
biomass = sol.f;

% obtenemos composicion del medio
media = getMediaCompsFromModel(model);

% defined media descrito por Salvador et al. (2018) Insights into metabolic
% osmoadaptation of C. salexigens
% % condicion high salinity
model = changeRxnBounds(model,'EX_glc-D_e',-2.1,'l');
model = changeRxnBounds(model,'EX_pyr_e',-0.3,'l');
model = changeRxnBounds(model,'EX_ac_e',-0.02,'l');
model = changeRxnBounds(model,'EX_nh4_e',-2.48,'l');
model = changeRxnBounds(model,'EX_glcn_e',0.17,'l');

sol1 = optimizeCbModel(model);
if sol1.f > 0
    fprintf('El modelo produce %f de biomasa',sol1.f)
end

% cambiamos nombre a ec biomasa
model.rxns(getPosOfElementsInArray({'BIOMASS_high salinity '},model.rxns)) = {'BIOMASS_high_salinity'};
% guardamos modelo
writeCbModel(model,'C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs\after_gapfill\campaniensis\campaniensis_gapfill_checked_v2.mat')

end