% Script pruebas produccion modelo
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs');
% leemos modelo
model = readCbModel([baseFolder filesep 'manual_adds' filesep 'campaniensis_w_phb_v5-1.mat']);
% length(model.mets)
% length(model.rxns)

%% add rxns based on pangenomic info
[data1] = get_pangenomic_info_from_drafts();
[data2] = get_pangenomic_info_from_gapfill();
dif1 = setdiff(data1.Code,model.rxns);
dif2= setdiff(data2.Code,model.rxns);

rxnsToAdd.rxns = data1.Code;
rxnsToAdd.rxnNames = data1.Description;
rxnsToAdd.rxnformula = data1.Reaction;
rxnsToAdd.subsystem = data1.Subsystem;
rxnsToAdd.gpr = data1.GPR;
rxnsToAdd.lb = data1.lb;
rxnsToAdd.ub = data1.ub;
% agregamos todas las rxns
for i = 1:length(rxnsToAdd.rxns)
    model = addReaction(model,rxnsToAdd.rxns{i},'reactionName',rxnsToAdd.rxnNames{i},...
        'reactionFormula',rxnsToAdd.rxnformula{i},'subSystem',rxnsToAdd.subsystem{i},...
        'lowerBound',rxnsToAdd.lb(i),'upperBound',rxnsToAdd.ub(i));
end

if isempty(dif2)
    disp("All pangenomic reactions from the gapfill are already present in the model")
else
    disp("Added pangenomic reactions to the model")
end
%% convert mets to bigg codes
dataFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\manual_adds');
data = readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],...
    'Sheet','met_equivalencies', 'Range','A1:E140');

for i = 1:length(table2cell(data(:,1)))
model = changeMetIdentifier(model,data.Actual_code{i},data.New_code{i});
model.metKEGGID(getPosOfElementsInArray(data.New_code(i),model.mets)) = data.kegg_code(i);
end

length(model.mets)
length(model.rxns)

%% fixing metFormulas and KeggID
for i = 1:length(table2cell(data(:,1)))
model = changeMetFormula(model,data.New_code{i},data.met_formula{i});
end
model.metFormulas = erase(model.metFormulas,'.');
%% add new reactions from kegg
folder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\manual_adds');
% new mets
data1 = readtable([folder filesep 'new_rxns_kegg_campaniensis.xlsx'],'Sheet','new_mets');
model = addMultipleMetabolites(model,data1.Code,'metNames',data1.Description,'metCharges',...
data1.Charge, 'metFormulas', data1.Charged_formula, 'metKEGGID',data1.kegg_code);

% new rxns
data2 = readtable([folder filesep 'new_rxns_kegg_campaniensis.xlsx'],'Sheet','new_rxns');
rxnsToAdd.rxns = data2.Code;
rxnsToAdd.rxnNames = data2.Description;
rxnsToAdd.rxnformula = data2.Reaction;
rxnsToAdd.subsystem = data2.Subsystem;
rxnsToAdd.gpr = data2.GPR;
% agregamos todas las rxns
for i = 1:length(rxnsToAdd.rxns)
    model = addReaction(model,rxnsToAdd.rxns{i},'reactionName',rxnsToAdd.rxnNames{i},...
        'reactionFormula',rxnsToAdd.rxnformula{i},'subSystem',rxnsToAdd.subsystem{i});%,...
        %'geneRule', rxnsToAdd.gpr{i});
end
% homologamos genes y proteinas
%model.geneNames{length(model.genes),1}=[];
% cuando agregamos los genes nos genera problemas para el modelo, resolver
% dps
%% guardamos el modelo
writeCbModel(model,'fileName', [baseFolder filesep 'manual_adds' filesep 'campaniensis_manual_fix_v6-1.mat'])
%writeCbModel(model,'fileName', [baseFolder filesep 'manual_adds' filesep 'campaniensis_manual_fix_v6-1.xml'])
%% vemos si produce phb desde distintas fuentes C que vimos experimentalmente
if testPathway(model, 'glc-D_e','phb_e') > 0 
    disp('Produce con glucosa')
end
if testPathway(model, 'malt_e','phb_e') > 0 
    disp('Produce con maltosa')
end
if testPathway(model, 'fru_e','phb_e') > 0
    disp('Produce con fructosa')
end
if testPathway(model, 'glyc_e','phb_e') > 0
    disp('Produce con glicerol')
end
if testPathway(model, 'sucr_e','phb_e') > 0
    disp('Produce con sacarosa')
end
%% 
% exc_rxn = model.rxns(383:526);
% prueba = changeRxnBounds(model,exc_rxn,-10,'l');
%gapFind(prueba)
%% generamos un nuevo medio
media = getMediaCompsFromModel(model); % bloquear compuestos nitrogenados
% dejamos todo el medio en -10
model = changeRxnBounds(model, media.reactions, -10,'l');
nitro = {'EX_ura_e','EX_pyr_e','EX_sucr_e'};
model = changeRxnBounds(model,nitro,0,'l');

%% acoplamiento N y C
 model.lb(getPosOfElementsInArray({'EX_nh4_e'},model.rxns)) = -1; % nh4
 model.lb(getPosOfElementsInArray({'EX_glc-D_e'},model.rxns)) = -2; % glc-D
% % optimizamos biomasa
% model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 0;
sol2 = optimizeCbModel(model);
% valor de objetivo
sol2.x(getPosOfElementsInArray({'EX_glc-D_e'},model.rxns)) % glc-D
sol2.x(getPosOfElementsInArray({'PHB-exchange'},model.rxns)) % phb
sol2.x(getPosOfElementsInArray({'EX_nh4_e'},model.rxns)) % nh4
sol2.x(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) % biomass

% N y C no estan acoplados

%%
% limitamos nitrogeno y optimizamos phb
%model = changeRxnBounds(model,'EX_nh4_e',-5,'l');
model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 1;
sol3 = optimizeCbModel(model);
sol3.x(getPosOfElementsInArray({'EX_glc-D_e'},model.rxns)) % glc-D
sol3.x(getPosOfElementsInArray({'PHB-exchange'},model.rxns)) % phb
sol3.x(getPosOfElementsInArray({'EX_nh4_e'},model.rxns)) % nh4
sol3.x(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) % biomass


% revisamos costos reducidos
src = sol3.w.*sol3.x; % scaled reduced cost (multiplicacion de costos reducidos * flujos)
pos_scr = find(sol3.w.*sol3.x);
rc = sol3.w;
rc3 = [model.rxns(pos_scr) getRxn_cobraFormat(model,pos_scr) num2cell(rc(pos_scr)) num2cell(src(pos_scr)) num2cell(sol3.x(pos_scr))];

%% revisamos el template
Folder2 = ('G:\My Drive\Magister\Metadraft\Templates\v2');
% leemos modelo
model2 = readCbModel([Folder2 filesep 'model_c_sal.mat']);
model2 = changeRxnBounds(model2,'EX_nh4_e',-5,'l');
media = getMediaCompsFromModel(model2);
sol4 = optimizeCbModel(model2);
sol4.x(getPosOfElementsInArray({'EX_glc-D_e'},model.rxns)) % glc-D
sol4.x(getPosOfElementsInArray({'EX_nh4_e'},model.rxns)) % nh4
sol4.x(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) % biomass
src = sol4.w.*sol4.x; % scaled reduced cost (multiplicacion de costos reducidos * flujos)
pos_scr = find(sol4.w.*sol4.x);
rc = sol4.w;
rc3 = [model2.rxns(pos_scr) getRxn_cobraFormat(model2,pos_scr) num2cell(rc(pos_scr)) num2cell(src(pos_scr)) num2cell(sol3.x(pos_scr))];

% template no tiene problemas de limitacion N