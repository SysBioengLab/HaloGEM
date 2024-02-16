function add_unique_Rxns_campaniensis
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs\after_gapfill\campaniensis');
% leemos modelo
model = readCbModel([baseFolder filesep 'campaniensis_gapfill_checked_v2.mat']);
% leemos rxns unicas
% usamos las rxns formato kegg pq coinciden mas id con el modelo
% las otras se revisan a mano posteriormente
% rxns1 = readcell([baseFolder filesep 'unique_rxns_camp_kegg.xls']);
unique_rxns = load([baseFolder filesep 'unique_rxns_camp_kegg_translated_full_v5.mat']);
rxnsToAdd = unique_rxns.unique_rxns;

% rxns3 = readcell([baseFolder filesep 'unique_rxns_camp_metacyc.xls']); % no se usa pq id no coincide, ver a mano

% consideramos mets no nulos
R = rxnsToAdd.unique_mets;
rxnsToAdd.unique_mets = R(~cellfun('isempty',R)); 

rep = intersect(rxnsToAdd.rxns,model.rxns);% hay una rxn repetida
% eliminamos la repetida
pos = getPosOfElementsInArray(rep,rxnsToAdd.rxns);
rxnsToAdd.rxns(pos) = [];
rxnsToAdd.rxn_formula(pos) = [];
rxnsToAdd.rxnNames(pos) = [];
rxnsToAdd.mets4rxn(pos,:) = [];
rxnsToAdd.stoich_coef(pos) = [];
rxnsToAdd.lb(pos) = [];
rxnsToAdd.ub(pos) = [];

% buscamos metabolitos no existentes en el modelo 
r = load([baseFolder filesep 'unique_mets_camp_kegg_translated_full_v2.mat']);
unique_mets = r.unique_mets;
new_mets = setdiff(unique_mets.mets,model.mets);
new_metFormulas = unique_mets.metFormulas(getPosOfElementsInArray(new_mets,unique_mets.mets));
new_metNames = unique_mets.metNames(getPosOfElementsInArray(new_mets,unique_mets.mets)); 
% agregamos mets nuevos
% % To add metabolites, with charges, formulas and KEGG ids:
% model = addMultipleMetabolites(model,{'A','b','c'},'metCharges',...
% [ -1 1 0], 'metFormulas', {'C','CO2','H2OKOPF'}, 'metKEGGID',{'C000012','C000023','C000055'})
model_new = addMultipleMetabolites(model,new_mets,'metNames',new_metNames,'metFormulas',new_metFormulas);
% agregamos rxns nuevas
% 'newRxn1','metaboliteList',{'A','B','C'},'stoichCoeffList',[-1 1 2], 'reversible',false);
%%
% agregamos todas las rxns
for i = 1:length(rxnsToAdd.rxns)
    model_new = addReaction(model_new,rxnsToAdd.rxns{i},'reactionName',rxnsToAdd.rxnNames{i},...
        'reactionFormula',rxnsToAdd.rxn_formula{i},...
        'lowerBound',rxnsToAdd.lb(i),...
        'upperBound',rxnsToAdd.ub(i));
    model_new.subSystems(length(model.rxns)+i) = unique_rxns.unique_rxns.subSystems(i);
end
%writeCbModel(model_new);% guardamos como excel
writeCbModel(model_new,'fileName', [baseFolder filesep 'campaniensis_w_unique_rxns_v5-1.mat']);
end