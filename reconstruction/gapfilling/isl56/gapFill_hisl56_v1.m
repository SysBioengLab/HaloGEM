function gapFill_hisl56_v1
%load('G:\My Drive\Magister\Metadraft\Templates\v2\model_Caro.xml');
model = readCbModel('G:\My Drive\Magister\Metadraft\Templates\v2\model_Caro.xml');
database = model;
% los mets y lo cambia por _e
%database.mets = regexprep(database.mets,{'\[e\]$','\[p\]$','\[c\]$'},{'_e','_p','_c'});
database = createGrRulesFromRules(database);
database = transformModelToCBMPYFormat(database);
media = getMediaFromModel(database);
database = changeRxnBounds(database, media.reactions, -10,'l');
database = changeRxnBounds(database,'BIOMASS_high salinity ',0,'l');
database = changeRxnBounds(database,'BIOMASS_high salinity ',1000,'u');
database.c = zeros(size(database.c));
database.c(getPosOfElementsInArray({'BIOMASS_high salinity '},database.rxns)) =1;
database = removeRxns(database,{'BIOMASS_low salinity '});
database.ub(intersect(find(database.lb<0), find(database.ub<0))) = 0;
pos = intersect(find(database.lb>0), find(database.ub>0));
database.lb(pos) = 0;
database.ub(pos) = 1000;


fba_database = optimizeCbModel(database);

draft = readCbModel('G:\My Drive\Magister\Metadraft\drafts\v2\h_asl10.xml','fileType','SBML');% los mets y lo cambia por _e
%draft = h_bol_lc1;
draft = transformModelToCBMPYFormat(draft);
% draft = addReactionFromModelRef(model, {'biomass_LPL60'},database);
% draft.rxns{getPosOfElementsInArray({'lpl_biomass'},draft.rxns)} = 'biomass_LPL60';
draft.c(getPosOfElementsInArray({'BIOMASS_high salinity '},draft.rxns)) =1;
draft = removeRxns(draft,{'BIOMASS_low salinity '});
pos_growth = find(database.c);
growth = database.rxns{pos_growth};
fbaWT = optimizeCbModel(database);


media2 = getMediaFromModel(database);

%%
constrOpt = struct('rxnList', {{growth}},'values', 0.01*fbaWT.f, 'sense', 'G');
constraints = constrOpt;

[modelGapfilled, solutions_rxns, solutions, lengthSolutions, translatedSolutions] = gapFill(draft, database, constraints,0);


save('modelGapfilled_001_wp_2','modelGapfilled')

% model_gapfilled = removeRxns(database, translatedSolutions{1});
fba_check1 = optimizeCbModel(modelGapfilled);
added = setdiff(modelGapfilled.rxns,draft.rxns);

for i = 1:length(added); 
    fba_i = optimizeCbModel(removeRxns(modelGapfilled, added(i))); disp(fba_i.f); 
    if fba_i.f>fba_check1.f*0.9; 
        disp(added(i)); 
    end; 
end


info_added = [added getRxn_cobraFormat(database,added)];
info_full = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled)];
info_solution = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled) num2cell(fba_check1.x)];
gene_associatios = [added database.grRules(getPosOfElementsInArray(added, database.rxns))];
info_media = [media2.reactions, num2cell(media2.lb), num2cell(media2. ub)];


xlswrite('gap_filling_v2_001_wp_isl56',info_added,'added')
xlswrite('gap_filling_v2_001_wp_isl56',info_full,'full_model')
xlswrite('gap_filling_v2_001_wp_isl56',info_solution,'solution')
xlswrite('gap_filling_v2_001_wp_isl56',gene_associatios,'gene_asso')
xlswrite('gap_filling_v2_001_wp_isl56',info_media,'media')


%%
constrOpt = struct('rxnList', {{growth}},'values', 0.8*fbaWT.f, 'sense', 'G');
constraints = constrOpt;

[modelGapfilled, solutions_rxns, solutions, lengthSolutions, translatedSolutions] = gapFill(draft, database, constraints);

save('modelGapfilled_01_2','modelGapfilled')
% model_gapfilled = removeRxns(database, translatedSolutions{1});
fba_check2 = optimizeCbModel(modelGapfilled);
added = setdiff(modelGapfilled.rxns,draft.rxns);

info_added = [added getRxn_cobraFormat(database,added)];
info_full = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled)];
info_solution = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled) num2cell(fba_check2.x)];
gene_associatios = [added database.grRules(getPosOfElementsInArray(added, database.rxns))];
info_media = [media2.reactions, num2cell(media2.lb), num2cell(media2. ub)];


xlswrite('gap_filling_v2_08_isl56',info_added,'added')
xlswrite('gap_filling_v2_08_isl56',info_full,'full_model')
xlswrite('gap_filling_v2_08_isl56',info_solution,'solution')
xlswrite('gap_filling_v2_08_isl56',gene_associatios,'gene_asso')
xlswrite('gap_filling_v2_08_isl56',info_media,'media')

%%
constrOpt = struct('rxnList', {{growth}},'values', 1*fbaWT.f, 'sense', 'G');
constraints = constrOpt;

[modelGapfilled, solutions_rxns, solutions, lengthSolutions, translatedSolutions] = gapFill(draft, database, constraints);


save('modelGapfilled_1_2','modelGapfilled')

% model_gapfilled = removeRxns(database, translatedSolutions{1});
fba_check3 = optimizeCbModel(modelGapfilled);
added = setdiff(modelGapfilled.rxns,draft.rxns);

info_added = [added getRxn_cobraFormat(database,added)];
info_full = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled)];
info_solution = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled) num2cell(fba_check3.x)];
gene_associatios = [added database.grRules(getPosOfElementsInArray(added, database.rxns))];
info_media = [media2.reactions, num2cell(media2.lb), num2cell(media2. ub)];


xlswrite('gap_filling_v1_2_isl56',info_added,'added')
xlswrite('gap_filling_v1_2_isl56',info_full,'full_model')
xlswrite('gap_filling_v1_2_isl56',info_solution,'solution')
xlswrite('gap_filling_v1_2_isl56',gene_associatios,'gene_asso')
xlswrite('gap_filling_v1_2_isl56',info_media,'media')

%%
constrOpt = struct('rxnList', {{growth}},'values', 0.01*fbaWT.f, 'sense', 'G');
constraints = constrOpt;

[modelGapfilled, solutions_rxns, solutions, lengthSolutions, translatedSolutions] = gapFill(draft, database, constraints);


save('modelGapfilled_001_2','modelGapfilled')

% model_gapfilled = removeRxns(database, translatedSolutions{1});
fba_check4 = optimizeCbModel(modelGapfilled);
added = setdiff(modelGapfilled.rxns,draft.rxns);

info_added = [added getRxn_cobraFormat(database,added)];
info_full = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled)];
info_solution = [modelGapfilled.rxns getRxn_cobraFormat(modelGapfilled) num2cell(fba_check4.x)];
gene_associatios = [added database.grRules(getPosOfElementsInArray(added, database.rxns))];
info_media = [media2.reactions, num2cell(media2.lb), num2cell(media2. ub)];


xlswrite('gap_filling_v2_001_isl56',info_added,'added')
xlswrite('gap_filling_v2_001_isl56',info_full,'full_model')
xlswrite('gap_filling_v2_001_isl56',info_solution,'solution')
xlswrite('gap_filling_v2_001_isl56',gene_associatios,'gene_asso')
xlswrite('gap_filling_v2_001_isl56',info_media,'media')



end