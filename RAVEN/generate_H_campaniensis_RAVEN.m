function generate_H_campaniensis_RAVEN 
%% csa; Chromohalobacter salexigens DSM 3043
% misma que la cepa usada para la generación del template
baseFolder = ('G:\My Drive\Magister\Metadraft');
%regexprep(which('init_h_campaniensis_toolbox.m'), ['\' filesep 'init_h_campaniensis_toolbox.m'],'');
% % fasta proteina del template
% inputFile = [baseFolder filesep 'Templates' filesep 'v2' filesep 'GCF_000055785.1_ASM5578v1_protein.faa'];
% 
% % codigo de KEGG del template, get kegg model
% model=getKEGGModelForOrganism('csa',inputFile,'prok90_kegg100','csa');
% save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_kegg.mat'],'model')
% 
% % get metacyc model
% model = getMetaCycModelForOrganism( 'csa_raven_MetaCyc',inputFile,1);
% save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_metacyc.mat'],'model')

% loading data
load([baseFolder filesep filesep 'Templates' filesep 'v2' filesep 'model_c_sal.mat'])
rootFolder = fileparts(which('initSystemsBioinformaticsToolbox.m'));
load([rootFolder filesep 'MNX' filesep 'metOtherIDs.mat'])
load([rootFolder filesep 'BIGG' filesep 'bigg_85_more_fixediNF517.mat'])
bigg = bigg_more_fixediNF517;
load([rootFolder filesep 'MNX' filesep 'metMNXIDs.mat'])

load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_kegg.mat'])
% corri hasta aca bien
model = refineRAVENModel(model,'kegg',bigg,metOtherIDs,metMNXIDs);
model = removeDuplicatedMetabolitesFromModel(model);
model = removeEmptyRxns(model);
model_translated_to_bigg_just_mets = translateModelToTargetLanguage(model, 'kegg', 'bigg', bigg, metOtherIDs, metMNXIDs, 1);
model_translated_to_bigg_full =      translateModelToTargetLanguage(model, 'kegg', 'bigg', bigg, metOtherIDs, metMNXIDs, 0);

save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_kegg_translated_mets.mat'],'model_translated_to_bigg_just_mets')
save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_kegg_translated_full.mat'],'model_translated_to_bigg_full')

load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_metacyc.mat'])
model = refineRAVENModel(model, 'metacyc', bigg, metOtherIDs, metMNXIDs);
model = removeDuplicatedMetabolitesFromModel(model);
model = removeEmptyRxns(model);
model_translated_to_bigg_just_mets = translateModelToTargetLanguage(model, 'metacyc', 'bigg', bigg, metOtherIDs, metMNXIDs, 1);
model_translated_to_bigg_full =      translateModelToTargetLanguage(model, 'metacyc', 'bigg', bigg, metOtherIDs, metMNXIDs, 0);

save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_metacyc_translated_mets.mat'],'model_translated_to_bigg_just_mets')
save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'csa_ncbi_raven_metacyc_translated_full.mat'],'model_translated_to_bigg_full')


%% H campaniensis
inputFile = [baseFolder filesep 'genome' filesep 'H_campaniensis_5AG.faa'];

% model=getKEGGModelForOrganism('hcs',inputFile,'prok90_kegg100','hcs');
% save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_kegg.mat'],'model')
% % 
% model = getMetaCycModelForOrganism( 'hcs_raven_MetaCyc',inputFile,1);
% save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_metacyc.mat'],'model')

load([rootFolder filesep 'BIGG' filesep 'bigg_85_more_fixediNF517.mat'])
rootFolder = fileparts(which('initSystemsBioinformaticsToolbox.m'));
load([rootFolder filesep 'MNX' filesep 'metOtherIDs.mat'])
load([rootFolder filesep 'MNX' filesep 'metMNXIDs.mat'])


%% KEGG
load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_kegg.mat'])
model = refineRAVENModel(model, 'kegg', bigg, metOtherIDs, metMNXIDs);
model = removeDuplicatedMetabolitesFromModel(model);
model = removeEmptyRxns(model);
model_translated_to_bigg_just_mets = translateModelToTargetLanguage(model, 'kegg', 'bigg', bigg, metOtherIDs, metMNXIDs, 1);
model_translated_to_bigg_full =      translateModelToTargetLanguage(model, 'kegg', 'bigg', bigg, metOtherIDs, metMNXIDs, 0);

save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_kegg_translated_mets.mat'],'model_translated_to_bigg_just_mets')
save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_kegg_translated_full.mat'],'model_translated_to_bigg_full')

%% metacyc

load([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_metacyc.mat'])
model = refineRAVENModel(model, 'metacyc', bigg, metOtherIDs, metMNXIDs);
model = removeDuplicatedMetabolitesFromModel(model);
model = removeEmptyRxns(model);
model_translated_to_bigg_just_mets = translateModelToTargetLanguage(model, 'metacyc', 'bigg', bigg, metOtherIDs, metMNXIDs, 1);
model_translated_to_bigg_full =      translateModelToTargetLanguage(model, 'metacyc', 'bigg', bigg, metOtherIDs, metMNXIDs, 0);

save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_metacyc_translated_mets.mat'],'model_translated_to_bigg_just_mets')
save([baseFolder filesep 'outputs' filesep 'RAVEN' filesep 'hcs_ncbi_raven_metacyc_translated_full.mat'],'model_translated_to_bigg_full')



end