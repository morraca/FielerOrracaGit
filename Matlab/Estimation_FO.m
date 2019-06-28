% May 23rd 2018
% Author: Maria Jose Orraca

% Estimation of parameters in model

%%
clear 
clc
%cd('E:\maria_orraca\Documents\Personales\PaperCecilia')
cd('~/Documents/Fieler Orraca/201906 matlab')
%cd('C:\Users\maria\Dropbox\Fieler Orraca\201906 Matlab')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Step 1: Fix Parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
sigma=5;
rho=(sigma-1)/sigma;
Pareto_scale=1.5;
wu=1;
ws=2.92; % Skill premium considering hours worked
N=100000;
tau=1.2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Step 2: Fix alpha_kl*R_k + alpha_kh*R_k
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% I define sector 1 is export-oriented (those with total foreign sales/Total sales greater than 10%)
% and domestic-oriented, with total foreign sales/Total sales less than 10%. 
% From Economic Census Data, the total domestic sales in the sample for each 
% of those is given by

TotDomSalesData=[1869690342470,961377058980];
TotImportsData=[1077156383030,3450629447030];


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Step 3: Fix verctor of random variables u
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
rng(3333)
u=rand(N,2);
u=sort(u);


% FIRST COLUMN IS EXPORT-ORIENTED
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Step 4: Split firms according to productivity
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
FracProd=[0.9,0.9]; % I set this arbitrarily
auxExport=[0.0058899,0.0169077]; % Data: Fraction of exporters in each sector
FracExport=auxExport.*FracProd;
auxHighTask=[0.5363637,0.5463974]; % Data: Fraction of exporters that invest in the creation of new products.
FracHighTask=auxHighTask.*FracExport;

LowCutoff=round(N*(ones(1,2)-FracProd));
ExportCutoff=round(N*(ones(1,2)-FracExport));
HighCutoff=round(N*(ones(1,2)-FracHighTask));

AuxLowTask=zeros(N,2);
AuxExporter=zeros(N,2);
AuxHighTask=zeros(N,2);

AuxLowTask(LowCutoff(1):end,1)=1;
AuxLowTask(LowCutoff(2):end,2)=1;
AuxHighTask(HighCutoff(1):end,1)=1;
AuxHighTask(HighCutoff(2):end,2)=1;
AuxExporter(ExportCutoff(1):end,1)=1;
AuxExporter(ExportCutoff(2):end,2)=1;

ActiveDom_nonexp_orient=[AuxLowTask(:,1),AuxHighTask(:,1)];
ActiveFor_nonexp_orient=[AuxExporter(:,1),AuxHighTask(:,1).*AuxExporter(:,1)];

ActiveDom_exp_orient=[AuxLowTask(:,2),AuxHighTask(:,2)];
ActiveFor_exp_orient=[AuxExporter(:,2),AuxHighTask(:,2).*AuxExporter(:,2)];

ActiveDom=reshape([ActiveDom_nonexp_orient,ActiveDom_exp_orient],[N,2,2]);
ActiveFor=reshape([ActiveFor_nonexp_orient,ActiveFor_exp_orient],[N,2,2]);
NonExporters=ActiveDom.*(ActiveFor==0);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Step 5: Guess vector of seven parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
varsigma_init=2;
beta_l_init=0.15;
beta_h_init=0.66;
a_init=7.5;
alpha_h_init=0.5;
Qstar_l_init=0.4;
Qstar_h_init=0.2;

Params_init=[varsigma_init,beta_l_init,beta_h_init,a_init,alpha_h_init];



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Steps 6-
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
numberFirms_data=[186761,108353];


% Moments in data

AvgSkintNonExp=[0.0616771,0.0968911];  %s/u
delta1dat=[0.0698,0.055];  %s/u with ln sales
meanSkintExporter=[0.6563934,0.3728212];  %s/u
gamma1dat=[-0.42,-0.338];  % lndepvar s/u
RelativeDomSalesNonexpExpdat=[0.005855117,0.012085438];
AvgXShare=[0.2804691,0.3070189];
mean_sd_domsalesdat=[0.02557199,0.047840274];
NumFirms=[186761,108353];


%%
mom_delta=NaN(1,2,8);
mom_avgSkintNonExp=NaN(1,2,8);
mom_avgSkillExp=NaN(1,2,8);
mom_RelDomSalesNonExpExp=NaN(1,2,8);
mom_mean_sd_sales=NaN(1,2,8);
mom_avgXshare=NaN(1,2,8);
mom_gamma=NaN(1,2,8);

diff=NaN(1,2,8);
QuantFstar=NaN(2,2,8);
UnskTot=NaN(1,2,8);
SkTot=NaN(1,2,8);
fll=NaN(2,2,8);
fxx=NaN(2,2,8);
fhh=NaN(1,2,8);
meanNetProfTot=NaN(1,2,8);
posfll=NaN(2,2,8);
posfxx=NaN(2,2,8);
Pstar_dom=NaN(2,2,8);
RevFor_mat=NaN(2,2,8);
XShareByFirm_mat=NaN(N,2,8);
SkillIntensityByFirm_mat=NaN(N,2,8);


penaliz=NaN(1,2,8);

ParamsEstim=NaN(5,2,8);



%% Estimations

eta_vec=[1.6,2,2.5,4];
wht_momskintsales_vec=[1,10];


for ii=1:1:length(eta_vec)
    eta=eta_vec(ii);
    fix_par=[sigma,rho,eta];

    for jj=1:1:2
        weight_momskintsales=wht_momskintsales_vec(jj);
        
        modd=ii*2+jj-2;
        
       Equilibrium_kFctParams=@(Pparams) Equilibrium(u(:,1),Pparams,fix_par,tau,Pareto_scale,wu,ws,N,TotDomSalesData(1),...
            ActiveDom(:,:,1),ActiveFor(:,:,1),NonExporters(:,:,1),...
            AvgSkintNonExp(1), delta1dat(1), meanSkintExporter(1),...
            RelativeDomSalesNonexpExpdat(1), AvgXShare(1), gamma1dat(1), mean_sd_domsalesdat(1),NumFirms(1),weight_momskintsales);
        
        ParamsEstim(:,1,modd)=fminsearch(Equilibrium_kFctParams,Params_init);
        [diff(1,1,modd),QuantFstar(:,1,modd),UnskTot(1,1,modd),SkTot(1,1,modd),fll(:,1,modd),fxx(:,1,modd),fhh(1,1,modd),...
            meanNetProfTot(1,1,modd),posfll(:,1,modd),posfxx(:,1,modd),mom_delta(1,1,modd),...
            mom_avgSkintNonExp(1,1,modd),mom_avgSkillExp(1,1,modd),mom_RelDomSalesNonExpExp(1,1,modd),...
            mom_mean_sd_sales(1,1,modd),mom_avgXshare(1,1,modd), mom_gamma(1,1,modd),Pstar_dom(:,1,modd),...
            RevFor_mat(:,1,modd),XShareByFirm_mat(:,1,modd),SkillIntensityByFirm_mat(:,1,modd),penaliz(1,1,modd)]=...
            Equilibrium_kFctParams(ParamsEstim(:,1,modd));
        
        Equilibrium_kFctParams=@(Pparams) Equilibrium(u(:,2),Pparams,fix_par,tau,Pareto_scale,wu,ws,N,TotDomSalesData(2),...
            ActiveDom(:,:,2),ActiveFor(:,:,2),NonExporters(:,:,2),...
            AvgSkintNonExp(2), delta1dat(2), meanSkintExporter(2),...
            RelativeDomSalesNonexpExpdat(2), AvgXShare(2), gamma1dat(2), mean_sd_domsalesdat(2),NumFirms(2),weight_momskintsales);

        ParamsEstim(:,2,modd)=fminsearch(Equilibrium_kFctParams,Params_init);
        [diff(1,2,modd),QuantFstar(:,2,modd),UnskTot(1,2,modd),SkTot(1,2,modd),fll(:,2,modd),fxx(:,2,modd),fhh(1,2,modd),...
            meanNetProfTot(1,2,modd),posfll(:,2,modd),posfxx(:,2,modd),mom_delta(1,2,modd),...
            mom_avgSkintNonExp(1,2,modd),mom_avgSkillExp(1,2,modd),mom_RelDomSalesNonExpExp(1,2,modd),...
            mom_mean_sd_sales(1,2,modd),mom_avgXshare(1,2,modd), mom_gamma(1,2,modd),Pstar_dom(:,2,modd),...
            RevFor_mat(:,2,modd),XShareByFirm_mat(:,2,modd),SkillIntensityByFirm_mat(:,2,modd),penaliz(1,2,modd)]=...
            Equilibrium_kFctParams(ParamsEstim(:,2,modd));
    end
end

%
% Calculating Pforeignstar
% Calculating price of foreign goods sold in domestic market
Pstar_for=NaN(2,2,length(eta_vec));
diff_pstar=NaN(1,2,length(eta_vec));

for modd=1:1:length(eta_vec)
    for iind= 1:1:2
        funPstar_fun=@(Ppstar_kl) funPstar(Ppstar_kl,TotImportsData(iind),ParamsEstim(5,iind,modd),...
        Pstar_dom(1,iind,modd),Pstar_dom(2,iind,modd),TotDomSalesData(iind),sigma,RevFor_mat(1,iind,modd),RevFor_mat(2,iind,modd));

        Pstar_for(1,iind,modd)=fsolve(funPstar_fun,Pstar_dom(1,iind,modd));

        [diff_pstar(1,iind,modd),Pstar_for(2,iind,modd)]=funPstar_fun(Pstar_for(1,iind,modd));

    end
end
%
% Bar graphs, skilli intensity and export share

for modd=1:1:length(eta_vec)
    Skint_TotFirms=[SkillIntensityByFirm_mat(:,1,modd);SkillIntensityByFirm_mat(:,2,modd)];
    XShare_TotFirms=[XShareByFirm_mat(:,1,modd);XShareByFirm_mat(:,2,modd)];

    Together=[Skint_TotFirms,XShare_TotFirms];
    ExportersInd=(Together(:,2)>0);
    SkintExporters=Skint_TotFirms(ExportersInd==1);
    XShareExporters=XShare_TotFirms(ExportersInd==1);
    TogetherExporter=[SkintExporters,XShareExporters];
    [~,idx] = sort(TogetherExporter(:,1)); % sort just the first column
    sorteTogetherExporter = TogetherExporter(idx,:);

    aux=floor(length(sorteTogetherExporter(:,2))/10);
    decileSkint=NaN(length(sorteTogetherExporter(:,2)),1);
    decileSkint(1:aux)=1;
    decileSkint(aux+1:2*aux)=2;
    decileSkint(2*aux+1:3*aux)=3;
    decileSkint(3*aux+1:4*aux)=4;
    decileSkint(4*aux+1:5*aux)=5;
    decileSkint(5*aux+1:6*aux)=6;
    decileSkint(6*aux+1:7*aux)=7;
    decileSkint(7*aux+1:8*aux)=8;
    decileSkint(8*aux+1:9*aux)=9;
    decileSkint(9*aux+1:end)=10;
   
    meanXShareBySkintdec=grpstats(sorteTogetherExporter(:,2),decileSkint);
    bar(meanXShareBySkintdec)
    saveas(gcf,sprintf('grXSh_bySkint_modd%d.png',modd));   
end

 
SU_economy=NaN(1,1,length(eta_vec));
for modd=1:1:length(eta_vec)
   SU_economy(1,1,modd)=(SkTot(1,1,modd)+SkTot(1,2,modd))./(UnskTot(1,1,modd)+UnskTot(1,2,modd))
end  

%%
FALTA ARRELGAR ESTOO PARA QUE EL 7 Y 8 QUEDEN EN TODOS
TAMBIEN FALTA EXPLICAR LOQ UE DECIDI DE PSTARFOR 
TAMBIEN FALTA HACER EL DOFILE PARA HACER EL BOOTSTRAP SAMPLE Y LLEVARLO AL INEGI!!!

T11=[ParamsEstim(:,1,1),ParamsEstim(:,1,3),ParamsEstim(:,1,4),ParamsEstim(:,1,5),ParamsEstim(:,1,6),...
    ParamsEstim(:,2,1),ParamsEstim(:,2,3),ParamsEstim(:,2,4),ParamsEstim(:,2,5),ParamsEstim(:,2,6),ParamsEstim(:,2,7),ParamsEstim(:,2,8)]
T12=[QuantFstar(:,1,1),QuantFstar(:,1,3),QuantFstar(:,1,4),QuantFstar(:,1,5),QuantFstar(:,1,6),...
    QuantFstar(:,2,1),QuantFstar(:,2,3),QuantFstar(:,2,4),QuantFstar(:,2,5),QuantFstar(:,2,6),QuantFstar(:,2,7),QuantFstar(:,2,8)]
T13=[fll(1,1,1),fll(1,1,3),fll(1,1,4),fll(1,1,5),fll(1,1,6),fll(1,2,1),fll(1,2,3),fll(1,2,4),fll(1,2,5),fll(1,2,6),fll(1,2,7),fll(1,2,8)]
T14=[fxx(1,1,1),fxx(1,1,3),fxx(1,1,4),fxx(1,1,5),fxx(1,1,6),fxx(1,2,1),fxx(1,2,3),fxx(1,2,4),fxx(1,2,5),fxx(1,2,6),fxx(1,2,7),fxx(1,2,8)]
T15=[fhh(1,1,1),fhh(1,1,3),fhh(1,1,4),fhh(1,1,5),fhh(1,1,6),fhh(1,2,1),fhh(1,2,3),fhh(1,2,4),fhh(1,2,5),fhh(1,2,6),fhh(1,2,7),fhh(1,2,8)]
T16=[meanNetProfTot(1,1,1),meanNetProfTot(1,1,3),meanNetProfTot(1,1,4),meanNetProfTot(1,1,5),meanNetProfTot(1,1,6),meanNetProfTot(1,2,1),...
    meanNetProfTot(1,2,3),meanNetProfTot(1,2,4),meanNetProfTot(1,2,5),meanNetProfTot(1,2,6),meanNetProfTot(1,2,7),meanNetProfTot(1,2,8)]


T2=[mom_delta(1,1,1),mom_delta(1,1,3),mom_delta(1,1,4),mom_delta(1,1,5),mom_delta(1,1,6)...
    mom_delta(1,2,1),mom_delta(1,2,3),mom_delta(1,2,4),mom_delta(1,2,5),mom_delta(1,2,6);...
    mom_avgSkintNonExp(1,1,1),mom_avgSkintNonExp(1,1,3),mom_avgSkintNonExp(1,1,4),mom_avgSkintNonExp(1,1,5),mom_avgSkintNonExp(1,1,6)...
    mom_avgSkintNonExp(1,2,1),mom_avgSkintNonExp(1,2,3),mom_avgSkintNonExp(1,2,4),mom_avgSkintNonExp(1,2,5),mom_avgSkintNonExp(1,2,6);...    
    mom_avgSkillExp(1,1,1),mom_avgSkillExp(1,1,3),mom_avgSkillExp(1,1,4),mom_avgSkillExp(1,1,5),mom_avgSkillExp(1,1,6)...
    mom_avgSkillExp(1,2,1),mom_avgSkillExp(1,2,3),mom_avgSkillExp(1,2,4),mom_avgSkillExp(1,2,5),mom_avgSkillExp(1,2,6);...    
    mom_mean_sd_sales(1,1,1),mom_mean_sd_sales(1,1,3),mom_mean_sd_sales(1,1,4),mom_mean_sd_sales(1,1,5),mom_mean_sd_sales(1,1,6)...
    mom_mean_sd_sales(1,2,1),mom_mean_sd_sales(1,2,3),mom_mean_sd_sales(1,2,4),mom_mean_sd_sales(1,2,5),mom_mean_sd_sales(1,2,6);...    
    mom_RelDomSalesNonExpExp(1,1,1),mom_RelDomSalesNonExpExp(1,1,3),mom_RelDomSalesNonExpExp(1,1,4),mom_RelDomSalesNonExpExp(1,1,5),mom_RelDomSalesNonExpExp(1,1,6)...
    mom_RelDomSalesNonExpExp(1,2,1),mom_RelDomSalesNonExpExp(1,2,3),mom_RelDomSalesNonExpExp(1,2,4),mom_RelDomSalesNonExpExp(1,2,5),mom_RelDomSalesNonExpExp(1,2,6);...
    mom_avgXshare(1,1,1),mom_avgXshare(1,1,3),mom_avgXshare(1,1,4),mom_avgXshare(1,1,5),mom_avgXshare(1,1,6)...
    mom_avgXshare(1,2,1),mom_avgXshare(1,2,3),mom_avgXshare(1,2,4),mom_avgXshare(1,2,5),mom_avgXshare(1,2,6);...
    mom_gamma(1,1,1),mom_gamma(1,1,3),mom_gamma(1,1,4),mom_gamma(1,1,5),mom_gamma(1,1,6)...
    mom_gamma(1,2,1),mom_gamma(1,2,3),mom_gamma(1,2,4),mom_gamma(1,2,5),mom_gamma(1,2,6)]

difftot=[diff(1,1,1),diff(1,1,3),diff(1,1,4),diff(1,1,5),diff(1,1,6),...
        diff(1,2,1),diff(1,2,3),diff(1,2,4),diff(1,2,5),diff(1,2,6)]

[SU_economy(1,1,1),SU_economy(1,1,3),SU_economy(1,1,4),SU_economy(1,1,5),SU_economy(1,1,6)] 

save('Estimations')





