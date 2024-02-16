function preProcessingHalomonas_v1 %(database)

rootFolder = fileparts(which('initSystemsBioinformaticsToolbox.m'));
load([rootFolder filesep 'BIGG' filesep 'bigg_85_more_fixediNF517.mat'])
database = bigg_more_fixediNF517;

tic
changeCobraSolver('glpk')
biomassEquation = 'BIOMASS_high_salinity';
%load('Streptococcus_model_NL');
% old = model;
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs');
load([baseFolder filesep 'manual_adds' filesep 'campaniensis_manual_fix_v6-1.mat'])
% model_cbmpy = readCbModel('campaniensis_w_phb_v3.xml');
%load('G:\My Drive\Magister\Metadraft\manual_adds\database.mat') % ubicacion de database

% model.grRules = cell(size(model.rxns));
% for i = 1:length(model.rxns) 
%     pos = find(strcmp(old.rxns, model.rxns(i)));
%     if ~isempty(pos) && length(pos)==1
%         model.grRules{i} = old.grRules{pos};
%     end
% end
% model = createGenesFromGrRules(model);
media = getMediaCompsFromModel(model);
model = changeRxnBounds(model, media.reactions, -10,'l');
model = changeObjective(model, biomassEquation);
model = changeRxnBounds(model, biomassEquation, 1000,'u'); 
fba1 = optimizeCbModel(model);

speciesName = 'Halomonas_campaniensis_5AG';
sufix = 'v6-1';
modelName = 'iHCS5AG_v6'; % todos parten con i, identificador kegg, cepa
baseFolder = pwd;
model = creategrRulesField(model); % genera campo GPR del modelo
model = fixMetIDsAndCompartments(model); % revisa los compartimentos de los mets

if ~isdir([baseFolder filesep 'Reports'])
    mkdir([baseFolder filesep 'Reports'])
end
bigg = database;
model.mets = regexprep(model.mets,'-','__');
cd([baseFolder filesep 'Reports'])
model= removeDuplicatedRxnsFromModel(model, 1, ['removedRxns1_' sufix]);
model = createRulesFromgrRules(model);
metrics = assessDatabaseCompatibility(model, bigg, ['reportBeforeCuration_' sufix]);
model = changeRxnsToMatchDatabase(model, metrics, 1, ['automaticReport1_' sufix]);
metrics2 = assessDatabaseCompatibility(model, bigg, ['reportAfterAutomaticCuration_' sufix]);
model= removeDuplicatedRxnsFromModel(model, 1, ['removedRxns2_' sufix]);
metrics3 = assessDatabaseCompatibility(model, bigg, ['reportAfterAutomaticDuplicateRxnsDeletion_' sufix]);
model = performGeneralIdentifierCuration(model, bigg,1, ['ManualCurationReport1_' sufix]);
metrics4 = assessDatabaseCompatibility(model, bigg, ['reportAfterManualCuration_' sufix]);
model = changeRxnsToMatchDatabase(model, metrics4, 1, ['automaticReport2_' sufix]);
metrics5 = assessDatabaseCompatibility(model, bigg, ['reportAfterAutomaticCuration2_' sufix]);
model= removeDuplicatedRxnsFromModel(model, 1, ['removedRxns3_' sufix]);
metrics6 = assessDatabaseCompatibility(model, bigg, ['reportAfterAutomaticDuplicateRxnsDeletion_' sufix]);
model = changeRxnsToMatchDatabase(model, metrics6, 1, ['automaticReport3_' sufix]);
metrics7 = assessDatabaseCompatibility(model, bigg, ['reportAfterAutomaticCuration3_' sufix]);

model = addChargesFromBigg(model, bigg);
model = creategrRulesField(model);
%model = create_rxnGeneMat(model);

if ~isdir([baseFolder filesep 'updatedModels'])
    mkdir([baseFolder filesep 'updatedModels'])
end
cd([baseFolder filesep 'updatedModels'])

format long
dt = datestr(datetime('now'));
dateToday = regexprep(dt, {':','-',' '},{'_','_','_'});
save([speciesName '_' sufix '_' dateToday], 'model');
save(modelName,'model');
model.mets = regexprep(model.mets,'(.*)_([ce])$' ,'$1\[$2\]');
save([modelName '_cobra_format'],'model');
writeCbModel(model, 'fileName', [modelName '_cobra_format'], 'format', 'sbml');
removeSquareBracketsFromBIGGModel([modelName '_cobra_format.xml'], [modelName '.xml']);
save([speciesName '_' sufix '_cobra_format_' dateToday],'model');
writeCbModel(model, 'fileName', [speciesName '_' sufix '_cobra_format_' dateToday], 'format', 'sbml')
removeSquareBracketsFromBIGGModel([speciesName '_' sufix '_cobra_format_' dateToday '.xml'], [speciesName '_' sufix '_' dateToday '.xml']);

model = changeObjective(model, biomassEquation);
model = changeRxnBounds(model, biomassEquation, 1000,'u'); 
fba2 = optimizeCbModel(model);
toc
end