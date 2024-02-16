% Script adición via PHB

baseFolder = ('C:\Users\carod\OneDrive - Universidad Católica de Chile\Magister\Metadraft\outputs');
% leemos modelo
model = readCbModel([baseFolder filesep 'after_gapfill' filesep 'campaniensis' filesep 'campaniensis_w_unique_rxns_v5-1.mat']);
model.rxns = regexprep(model.rxns,{'\[e\]$','\[p\]$','\[c\]$'},{'_e','_p','_c'});

%% optimizacion inicial
media = getMediaCompsFromModel(model);
model = changeRxnBounds(model,'EX_glc-D_e',-10,'l');
sol = optimizeCbModel(model);
% valor de objetivo
biomass = sol.f;

%% agregamos via phb

%%%(S)-3-Hydroxybutanoyl-CoA a (R)-3-Hydroxybutanoyl-CoA 5.1.2.3
% catalizado por (S)-3-Hydroxybutanoyl-CoA 3-epimerase
% NCBI-ProteinID: 	AIA75630
% UniProt: 	A0A060BEB0
% https://www.genome.jp/entry/K01782+K01825+5.1.2.3+R03276
% agregamos metabolito
% http://bigg.ucsd.edu/models/universal/metabolites/3hbcoa__R
model = addMetabolite(model,'3hbcoa__R_c','metName','(R)-3-Hydroxybutanoyl-CoA',...
    'metFormula','C25H42N7O18P3S','KEGGId','C03561');
% agregamos reaccion
model = addReaction(model,'3HBC3E','reactionName','(S)-3-Hydroxybutanoyl-CoA 3-epimerase',...
    'metaboliteList',{'3hbcoa_c','3hbcoa__R_c'},'stoichCoeffList',[-1 1], 'reversible',true);

%%%(R)-3-Hydroxybutanoyl-CoA a PHB 2.3.1.304
% catalizado por poly(3-hydroxyalkanoate) synthetase
% NCBI-ProteinID: 	AIA75487
% UniProt: 	A0A060B971
% https://www.genome.jp/entry/hcs:FF32_11725
% agregamos PHB como metabolito
% estamos asumiendo largo PHB = 4 monomeros, ver word supuestos
model = addMetabolite(model,'phb_c','metName','Poly-beta-hydroxybutyrate-PHB',...
    'metFormula','(C4H6O2)n','KEGGId','C06143');
% coenzima A esta como CoA en el modelo
% agregamos reaccion
model = addReaction(model,'PHB_syn_1','reactionName','poly(3-hydroxyalkanoate) synthetase',...
    'metaboliteList',{'3hbcoa__R_c','phb_c','coa_c'},'stoichCoeffList',[-4 1 4],...
    'reversible',true,'lowerBound',0);

%%%PHB Exchange
% no hay exchange pq lo acumula, lo usamos como artefacto
model = addMetabolite(model,'phb_e','metName','Poly-beta-hydroxybutyrate-PHB-exchange',...
     'metFormula','(C4H6O2)n','KEGGId','C06143');
model = addReaction(model,'PHB-trasport','metaboliteList',{'phb_c',...
     'phb_e'},'stoichCoeffList',[-1 1], 'reversible',false);
model = addReaction(model,'PHB-exchange','metaboliteList',{'phb_e'},'stoichCoeffList',[-1], 'reversible',false);


%%%Acetoacetyl-CoA a (R)-3-Hydroxybutanoyl-CoA
% catalizado por (R)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase
% NCBI-ProteinID: 	AIA74865
% UniProt: 	A0A060BC75
% https://www.genome.jp/entry/R01779/R01977
% http://bigg.ucsd.edu/universal/reactions/AACOAR_syn
model = addReaction(model,'AACOAR_syn','reactionName','(R)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase','metaboliteList',...
{'3hbcoa__R_c','nadp_c','aacoa_c','nadph_c','h_c'},'stoichCoeffList',...
[-1 -1 1 1 1], 'reversible',true);
     
%%%(S)-3-Hydroxybutanoyl-CoA a Acetoacetyl-CoA
% catalizado por (S)-3-Hydroxybutanoyl-CoA:NAD+ oxidoreductase
% RHEA: 	30802
% https://www.genome.jp/entry/R01778/R01975/R04203/R04737/R04739/R04741/R04743/R04745/R04748/R06941/R07890/R07894/R07898/R05066/R05305/R05575/R08094
% Ya está en el modelo, 3HYDEHY1

%%%(S)-3-Hydroxybutanoyl-CoA a Acetoacetyl-CoA
% catalizado por (S)-3-Hydroxybutanoyl-CoA:NADP+ oxidoreductase
% NCBI-ProteinID: 	AIA74137
% UniProt: 	A0A060B5G7
model = addReaction(model,'3BTCOAD','reactionName','(S)-3-Hydroxybutanoyl-CoA:NAD+ oxidoreductase','metaboliteList',...
{'3hbcoa_c','nadp_c','aacoa_c','nadph_c','h_c'},'stoichCoeffList',...
[-1 -1 1 1 1], 'reversible',true);

%%% acetoacetyl-CoA a acetyl-CoA
% ya está en el modelo, 3KETOTHIO9
%%% crotonoyl-CoA a (S)-3-Hydroxybutanoyl-CoA
% ya está en el modelo, 3HYDEHYDRA1
%%% butanoyl-CoA a crotonoyl-CoA 
% ya está en el modelo, ACYLDEHY1
% revisamos que modelo funcione bien
%% por datos experimentales sabemos crece con sacarosa, agregamos la vía
%%% Exchange Sucrose
model = addMetabolite(model,'sucr_c','metName','Sucrose','metFormula','C12H22O11',...
    'KEGGID','C00089');
model = addMetabolite(model,'sucr_e','metName','Sucrose','metFormula','C12H22O11',...
    'KEGGID','C00089');
model = addReaction(model,'SUCRtpp','reactionName','Sucrose exchange','metaboliteList',...
    {'sucr_e','atp_c','sucr_c','adp_c','pi_c'},'stoichCoeffList',[-1 -1 1 1 1],'reversible',true);
model = addExchangeRxn(model,'sucr_e');

%%% Sucrose a D-glucose
% catalizado por sucrose glucohydrolase
% NCBI-ProteinID: 	AIA75699
% UniProt: 	A0A060B964
% https://www.genome.jp/entry/hcs:FF32_12860+hcs:FF32_12890
model = addReaction(model,'SGHB','reactionName','sucrose glucohydrolase','metaboliteList',...
    {'sucr_c','h2o_c','fru_c','glc-D_c'},'stoichCoeffList',[-1 -1 1 1],'reversible',true);


%% Revisamos modelo
model.lb(491) = -1;
model.lb(440) = -20;
sol2 = optimizeCbModel(model);
% valor de objetivo
biomass2 = sol2.f;
sol2.x(440)
sol2.x(1504)
sol2.x(491)

%return

% vemos produccion phb
sol2.x(1502)
% si maximizamos produccion phb manteniendo biomasa
model.c(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) = 1;
model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 1;
sol3 = optimizeCbModel(model);
% valor de objetivo
biomass3 = sol3.f;

% vemos produccion phb
sol3.x(1502) % ! produce PHB ahora hay que revisar consistencia de la via

% si maximizamos produccion phb sin forzar la biomasa
model.c(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) = 0;
model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 1;
sol4 = optimizeCbModel(model);
% valor de objetivo
biomass4 = sol4.f;

% vemos produccion phb
sol4.x(1502)
% NO hay producción de PHB

% reseteamos a biomasa
model.c(getPosOfElementsInArray({'BIOMASS_high_salinity'},model.rxns)) = 1;
model.c(getPosOfElementsInArray({'poly(3-hydroxyalkanoate) synthetase'},model.rxnNames)) = 0;

writeCbModel(model,'fileName', [baseFolder filesep 'manual_adds' filesep 'campaniensis_w_phb_v5-1.mat'])
%writeCbModel(model,'fileName', [baseFolder filesep 'manual_adds' filesep 'campaniensis_w_phb_v5.xml'])
