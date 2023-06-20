function  getUniqueRxnsHcamp()
% determinar rxns unicas de campaniensis
% setdiff(A,B) gives the data of A that is not in B
%% %% kegg id
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\');
% hcs = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_kegg.mat']);
% csa = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_kegg.mat']);
% model = hcs.model;
% % difference
% shouldRemove1 = setdiff(model.rxns,csa.model.rxns);
% 
% % save to mat
% rxns1 = getRxn_cobraFormat(model,getPosOfElementsInArray(shouldRemove1, model.rxns));
% rxnNames = getRxn_cobraFormat(model,getPosOfElementsInArray(shouldRemove1,model.rxns),1);
% mets = getMetsFromRxns(model,1:length(shouldRemove1));
% info_full = [shouldRemove1 model.rxnNames(getPosOfElementsInArray(shouldRemove1,model.rxns))...
%     rxns1 rxnNames mets.mets mets.posMets ]; 
% xlswrite([baseFolder filesep 'outputs' filesep 'after_gapfill' filesep 'unique_rxns_camp_kegg_v2'],info_full);
Folder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs\after_gapfill\campaniensis');
% leemos modelo
model1 = readCbModel([Folder filesep 'campaniensis_gapfill_checked_v2.mat']);

%% kegg translated to bigg full
% load both models
hcs = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'campaniensis' filesep 'hcs_ncbi_raven_kegg_translated_full.mat']);
csa = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'campaniensis' filesep 'csa_ncbi_raven_kegg_translated_full.mat']);
model = hcs.model_translated_to_bigg_full; % modelo generado por RAVEN
% difference
shouldadd = setdiff(model.rxns,csa.model_translated_to_bigg_full.rxns);
shouldAddMets = cell(length(shouldadd),1); % hay mas rxns que metabolitos
helper = (setdiff(model.mets,csa.model_translated_to_bigg_full.mets));
shouldAddMets(1:length(helper)) = (setdiff(model.mets,csa.model_translated_to_bigg_full.mets));
setdiff(model.mets,model1.mets)
% calculate values
unique_rxns.rxns = shouldadd;
unique_rxns.rxn_formula = getRxn_cobraFormat(model,getPosOfElementsInArray(shouldadd,model.rxns),0);
unique_rxns.mets4rxn = getMetsFromRxns(model,getPosOfElementsInArray(shouldadd,model.rxns)).mets;
unique_rxns.stoich_coef = getRxnsStoichCoef(model,getPosOfElementsInArray(shouldadd,model.rxns));
unique_rxns.rxnNames = model.rxnNames(getPosOfElementsInArray(shouldadd,model.rxns));
unique_rxns.lb = model.lb(getPosOfElementsInArray(shouldadd,model.rxns));
unique_rxns.ub = model.ub(getPosOfElementsInArray(shouldadd,model.rxns));
unique_rxns.subSystems = model.subSystems(getPosOfElementsInArray(shouldadd,model.rxns));
unique_rxns.unique_mets = shouldAddMets;

unique_mets.mets = setdiff(model.mets,csa.model_translated_to_bigg_full.mets);
unique_mets.metNames = model.metNames(getPosOfElementsInArray(unique_mets.mets,model.mets));
unique_mets.metFormulas = model.metFormulas(getPosOfElementsInArray(unique_mets.mets,model.mets));
%unique_mets.KEGGID = model.
%data = [unique_mets.mets unique_mets.metNames unique_mets.metFormulas];
% xlswrite([baseFolder filesep 'outputs' filesep 'after_gapfill' filesep 'campaniensis' filesep...
%    'unique_mets_camp_kegg_translated_full_v1'],data);
mets = [baseFolder filesep 'outputs' filesep 'after_gapfill' filesep 'campaniensis' filesep 'unique_mets_camp_kegg_translated_full_v2'];
save(mets,'unique_mets');

%%
%info_full = [shouldadd rxn_formula rxns1 rxnNames (mets4rxn.mets) mets4rxn.posMets shouldAddMets stoich_coef]; 
%xlswrite([baseFolder filesep 'outputs' filesep 'after_gapfill' filesep ...
%    'unique_rxns_camp_kegg_translated_full'],info_full);
% save to mat
scalar = [baseFolder filesep 'outputs' filesep 'after_gapfill' filesep 'campaniensis' filesep 'unique_rxns_camp_kegg_translated_full_v5'];
save(scalar,'unique_rxns');
%% %% metacyc
% hcs = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_metacyc.mat']);
% csa = load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_metacyc.mat']);
% model = hcs.model;
% % difference
% shouldRemove3 = setdiff(model.rxns,csa.model.rxns);
% 
% % save to xls
% rxns1 = getRxn_cobraFormat(model,getPosOfElementsInArray(shouldRemove3, model.rxns));
% rxnNames = getRxn_cobraFormat(model,getPosOfElementsInArray(shouldRemove3,model.rxns),1);
% %mets = getMetsFromRxns(model,1:length(shouldRemove3));
% info_full = [shouldRemove3 model.rxnNames(getPosOfElementsInArray(shouldRemove3,model.rxns))...
%     rxns1 rxnNames]; 
% xlswrite([baseFolder filesep 'outputs' filesep 'after_gapfill' filesep 'unique_rxns_camp_metacyc_v2'],info_full);

end
