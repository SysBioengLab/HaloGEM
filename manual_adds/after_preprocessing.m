%% checking preprocessing
% version 6 adds genes
baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\manual_adds');
load([baseFolder filesep 'updatedModels' filesep 'iHCS5AG_v6.mat'])
model.rxns(getPosOfElementsInArray({'PHB-exchange'},model.rxns)) = {'EX_phb_e'};
%changeCobraSolver('glpk')
%% fixing nh3 to nh4
media = getMediaCompsFromModel(model);
model = changeRxnBounds(model, media.reactions, -10,'l');
model = changeMetIdentifier(model,'nh3_c','nh4_c');
model = changeMetFormula(model,'nh4_c','NH4');
model = removeMetabolites(model,'nh3_c');
model = changeRxnBounds(model,{'EX_nh4_e'},1000,'u');
%% acoplamiento C N

lb = {'EX_nh4_e','EX_glc__D_e','EX_o2_e'};
lower = model.lb(getPosOfElementsInArray(lb,model.rxns))';  
% optimizamos biomasa
model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 0;
model.c(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) = 1;
sol2 = optimizeCbModel(model);
% valor de objetivo
rxns = {'EX_nh4_e','EX_glc__D_e','EX_o2_e','BIOMASS_high_salinity','EX_phb_e'};
after_preprocess = [lower sol2.x(getPosOfElementsInArray(rxns,model.rxns))']';

% C y N no están acoplados, genera biomasa
model = removeDuplicatedMetabolitesFromModel(model);
%% generamos un nuevo medio
dataFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\manual_adds');
media2 = getMediaCompsFromModel(model);
model = changeRxnBounds(model, media2.reactions, 0,'l');
MMmedia = readtable([dataFolder filesep 'yeast_extract_composition.xlsx'],'Sheet','model_bounds');
% dejamos todo el medio en -10
model = changeRxnBounds(model, MMmedia.rxn, -MMmedia.mM,'l');

%% fix rxns based on escher map construction
% new mets
data1 = readtable([dataFolder filesep 'add_rxns_with_escher.xlsx'],'Sheet','new_mets');
model = addMultipleMetabolites(model,data1.Code,'metNames',data1.Description,'metCharges',...
data1.Charge, 'metFormulas', data1.Charged_formula, 'metKEGGID',data1.kegg_code);

% add exchange rxns
data1 = readtable([dataFolder filesep 'add_rxns_with_escher.xlsx'],'Sheet','exchanges');
model = addExchangeRxn(model,data1.Met,data1.lb,data1.ub);


% new rxns
data2 = readtable([dataFolder filesep 'add_rxns_with_escher.xlsx'],'Sheet','new_rxns');
rxnsToAdd.rxns = data2.Code;
rxnsToAdd.rxnNames = data2.Description;
rxnsToAdd.rxnformula = data2.Reaction;
rxnsToAdd.subsystem = data2.Subsystem;
rxnsToAdd.gpr = data2.GPR;
rxnsToAdd.lb = data2.lb;
rxnsToAdd.ub = data2.ub;
rxnsToAdd.ConfidenceScore = data2.Confidence_score;
rxnsToAdd.geneRule = data2.GPR;
prev_rxns = length(model.rxns);
% agregamos todas las rxns
for i = 1:length(rxnsToAdd.rxns)
    model = addReaction(model,rxnsToAdd.rxns{i},'reactionName',rxnsToAdd.rxnNames{i},...
        'reactionFormula',rxnsToAdd.rxnformula{i},'subSystem',rxnsToAdd.subsystem{i},...
        'lowerBound',rxnsToAdd.lb(i),'upperBound',rxnsToAdd.ub(i),'geneRule',...
        rxnsToAdd.geneRule{i});
    model.rxnConfidenceScores(prev_rxns+i) = rxnsToAdd.ConfidenceScore(i);
end
%% adding extra genes
genes_data= readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','genes');
for i= 1:length(genes_data.gene)
    rxn_to_link = model.rxns(getPosOfElementsInArray(genes_data.rxn_code(i),model.rxns));
    model = changeGeneAssociation(model, rxn_to_link, genes_data.gene{i});
end
%model.geneNames(765:(length(model.genes(765:end))+764)) = model.genes(765:end);
%model.proteins(765:(length(model.genes(765:end))+764)) = model.genes(765:end);
if isfield(model,'grRules')
    model = generateRules(model);
end

%% biomass medium salinity
[coef,compound] = biomass_medium_salinity(56.81,1.03); % 56.81 se obtiene interpolando coefs
% de low y high salinity
model = addReaction(model,'BIOMASS_medium_salinity','reactionName','Biomass medium salinity',...
    'metaboliteList',compound,'stoichCoeffList',coef,'lowerBound',0,...
    'upperBound',1000);

%% fixing rxns with bad stoi
data4 = readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','fix_rxns');
for i = 1:length(data4.Rxn)
    model = changeRxnStoi(model,data4.Rxn(i),data4.new_rxn{i});
end
%% add missing or fix met formulas
data3 = readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','met_formulas');

for i = 1:length(data3.Actual_code)
    model.metFormulas(getPosOfElementsInArray(data3.Actual_code(i),model.mets)) = data3.Formula(i);
    model.metCharges(getPosOfElementsInArray(data3.Actual_code(i),model.mets)) = data3.Charge(i);
    if isempty(model.metKEGGID{i})
        model.metKEGGID(getPosOfElementsInArray(data3.Actual_code(i),model.mets)) = data3.kegg_code(i);
    end
    if ~isempty(data3.Met_name{i})
        model.metNames(getPosOfElementsInArray(data3.Actual_code(i),model.mets)) = data3.Met_name(i);
    end
end
%% fix ATP reaction bounds
data5 = readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','bounds');
for i = 1:length(data5.Rxn)
    model.lb(getPosOfElementsInArray(data5.Rxn(i),model.rxns)) = data5.lb(i);
    model.ub(getPosOfElementsInArray(data5.Rxn(i),model.rxns)) = data5.ub(i);
end
%% remove rxns
data1 = readtable([baseFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','rxn_equivalencies',...
'Range','A1:B228');
rxns_to_remove = [];
for i= 1:length(data1.action)
    if (contains(data1.action(i),'Remove')==1) 
        rxns_to_remove = [rxns_to_remove data1.rxn_code(i)];
    end
end
model = removeRxns(model,rxns_to_remove); % removes metabolites not present in any reaction

%% fix extracellular pH to 9
data4 = readtable([dataFolder filesep 'manual_check_unique_campaniensis.xlsx'],'Sheet','fix_pH');

for i = 1:length(data4.Abbreviation)
    model.metFormulas(getPosOfElementsInArray(data4.Abbreviation(i),model.mets)) = data4.Charged_formula(i);
    model.metCharges(getPosOfElementsInArray(data4.Abbreviation(i),model.mets)) = data4.Charge(i);
end

%% fix biomass equations
data1 = readtable('biomass_equation.xlsx','Sheet','adjust');
high_coef = zeros(length(data1.compound),1);
low_coef = zeros(length(data1.compound),1);
for i=1:length(data1.compound)
    high_coef(i) = -data1.high(i); 
    low_coef(i) = -data1.low(i);
end
compound = data1.compound;
model = addReaction(model,'BIOMASS_high_salinity','reactionName','Biomass high salinity',...
    'metaboliteList',compound,'stoichCoeffList',high_coef);
model = addReaction(model,'BIOMASS_low_salinity','reactionName','Biomass low salinity',...
    'metaboliteList',compound,'stoichCoeffList',low_coef);

%% acoplamiento N y C
% %aa = {'EX_ala__D_e','EX_arg__L_e','EX_asp__L_e','EX_cys__L_e','EX_glu__L_e',...
%     'EX_gly_e','EX_his__L_e','EX_ile__L_e','EX_leu__L_e','EX_lys__L_e',...
%     'EX_met__L_e','EX_phe__L_e','EX_pro__L_e','EX_ser__L_e','EX_trp__L_e',...
%     'EX_tyr__L_e','EX_val__L_e'};
% model= changeRxnBounds(model,aa,0,'l')
% model.lb(getPosOfElementsInArray({'EX_nh4_e'},model.rxns))= 0;
% model.lb(getPosOfElementsInArray({'EX_glc__D_e'},model.rxns))= -10;
% model.lb(getPosOfElementsInArray({'EX_o2_e'},model.rxns))=-5;
lower = model.lb(getPosOfElementsInArray(lb,model.rxns))'; 

% optimizamos biomasa
model = changeObjective(model,{'BIOMASS_high_salinity'},1);
rxns = {'EX_nh4_e','EX_glc__D_e','EX_o2_e','BIOMASS_high_salinity','EX_phb_e'};
%media = getMediaCompsFromModel(model);
%table(media.reactions,media.lb);
%model.c(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) = 1;
sol2 = optimizeCbModel(model);
manual_adds = [lower sol2.x(getPosOfElementsInArray(rxns,model.rxns))']';

%
% model = changeRxnBounds(model,{'BIOMASS_high_salinity'},0,'l');
% limitamos nitrogeno y optimizamos phb
lower = model.lb(getPosOfElementsInArray(lb,model.rxns))'; 

%model = changeRxnBounds(model,'EX_nh4_e',-1,'l');
model = changeObjective(model,{'EX_phb_e'},1);
sol3 = optimizeCbModel(model);
%exportBIGGModelToSBML(model,'model_phb',1);

optPHB = [lower sol3.x(getPosOfElementsInArray(rxns,model.rxns))']';

% revisamos costos reducidos
src = sol2.w.*sol2.x; % scaled reduced cost (multiplicacion de costos reducidos * flujos)
pos_scr = find(abs(src)> 1e-6);
rc = sol2.w;
rc3 = [model.rxns(pos_scr) getRxn_cobraFormat(model,pos_scr) num2cell(rc(pos_scr)) num2cell(src(pos_scr)) num2cell(sol3.x(pos_scr))];

labels = {'nh4 lb','glc lb','o2 lb','nh4 used','glc used','o2 used','biom','phb'};
table(labels',after_preprocess,manual_adds,optPHB)


%% revisando
% estequiometria
%checkMassChargeBalance(model,-1)
%writeCbModel(model);