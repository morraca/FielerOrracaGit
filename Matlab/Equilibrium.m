%param=[sigma,rho,eta];
% This program does step 7: Calculate costs, skill intensity, price for
% each firm.

function [error,QuantFstar,UnskTot,SkTot,fl,fx,fh,meanNetProfTot,posfl,posfx,...
    delta1,meanSkillIntNonExporters,meanSkillIntExporters,RelativeDomSalesNonexpExp,mean_sd_sales,...
    avgXshare, gamma1,PriceIndx_Dom,RevFor,XShareByFirm,SkillIntensityByFirm,penaliz]=...
    Equilibrium(uu,Pparams,fix_par,tau,Pareto_scale,wwu,wws,NN,TotDomSales,...
    ActiveDom,ActiveFor,NonExporters, AvgSkintNonExpdat, delta1dat,meanSkillIntExportersdat ,...
    RelativeDomSalesNonexpExpdat, AvgXSharedat, gamma1dat, mean_sd_domsalesdat,nnumfirms,wht_momelastsales)
  
    aa=Pparams(4);
    if Pparams(2)>0 && Pparams(2)<1
        beta_l=Pparams(2);
    else
        beta_l=0.01;
    end
    if Pparams(3)>0 && Pparams(3)<1
        beta_h=Pparams(3);
    else 
        beta_h=0.99;
    end
    varsigma=Pparams(1);
    alpha_h=Pparams(5);
    alpha_l=1-alpha_h;
    alphas=[alpha_l,alpha_h];
    QuantsDom=repmat(TotDomSales*alphas,NN,1);
    sigma=fix_par(1);
    
    % Step 6: Transform u into vector of productivities
    Varphis=repmat(repmat(Pareto_scale,NN,1)./((ones(NN,1)-uu).^(1/aa)),1,2);
    
    % Step 7: Calculate costs, skill intensity, price for each firm.
    betas=repmat([beta_l,beta_h],NN,1);
    [~,pricedom,pricefor,skint,auxu]=FirmsChoices(Varphis,betas,varsigma,fix_par,tau,wwu,wws,NN);
    
    % Step 8: Calculate domestic sales:
    % auxu is a firm-level variable that allows to pin down the number of unskilled 
    % workers for a given quantity of production.
    pricedom=pricedom.*ActiveDom;
    pricefor=pricefor.*ActiveFor;
    pricedom(pricedom==0)=NaN;
    pricefor(pricefor==0)=NaN;
    
    PriceIndx_Dom=(nansum(pricedom.^(1-sigma))).^(1/(1-sigma));
    PriceIndx_Dom_mat=repmat(PriceIndx_Dom,NN,1);
    
    quantdom=QuantsDom.*((PriceIndx_Dom_mat).^(sigma-1)).*(pricedom.^(-sigma));
    quantdom(isnan(quantdom)==1)=0;
    revenuedom=pricedom.*quantdom;
    assert(round(sum(nansum(revenuedom)))==TotDomSales);
    revenuedom(isnan(revenuedom)==1)=0;    

    TotDomSales=sum(revenuedom,2);   
    unskilled_dom=ActiveDom.*quantdom./auxu;
    skilled_dom=unskilled_dom.*skint;
   
    % Initializing loop to get Qfstar
    ttol=0.01;
    chs_Xshare=1;
    chs_Xshare_h=0.5;
    diff_XShares=1;
    
    Qqstar_kp1=[1,1];
    it=1;

    while diff_XShares>ttol && it<1000
        QuantFstar=Qqstar_kp1;
        QuantFstar(QuantFstar<=0)=0.001;
        QuantFstarmat=repmat(QuantFstar,NN,1);
        quantfor=QuantFstarmat.*(pricefor).^(-sigma);
        quantfor(isnan(quantfor)==1)=0;
        revenuefor=pricefor.*quantfor;
        revenuefor(isnan(revenuefor)==1)=0; 

        TotForSales=sum(revenuefor,2);
        Xshare=TotForSales./(TotForSales+TotDomSales); % using definition of export share rather than
                                                    % relative foreign sales
        
        RelForSales=TotForSales./TotDomSales;
        
        % 5. Average export share of exporters
        avgXshare=mean(Xshare(ActiveFor(:,1)==1));

        unskilled_for=ActiveFor.*quantfor./auxu;
        skilled_for=unskilled_for.*skint;
        tot_unskilled=unskilled_dom+unskilled_for;    
        tot_skilled=skilled_dom+skilled_for; 

        SkillIntensityByFirm=sum(tot_skilled,2)./sum(tot_unskilled,2);

        SkillIntExporters=SkillIntensityByFirm((ActiveFor(:,1)==1));
        
        % 6. Regression of relative foreign sales on skill intensity
        lnRelForSales=log(RelForSales(ActiveFor(:,1)==1));
        
        XVars_xsh=[ones(length(Xshare(ActiveFor(:,1)==1)),1),log(SkillIntExporters)];
        Betas_Coefs_rfs=XVars_xsh\lnRelForSales;
        gamma1=Betas_Coefs_rfs(2);
        
        auxDiffsXShares_all=(AvgXSharedat-avgXshare).*chs_Xshare;
        DiffsXShares_all=kron(auxDiffsXShares_all,[1,1]);
                   
        auxDiffXShares_h=(gamma1dat-gamma1)*chs_Xshare_h;
        DiffXShares_h=kron(auxDiffXShares_h,[0,1]);
        Qqstar_kp1=QuantFstar.*(ones(1,2)+DiffXShares_h+DiffsXShares_all);
        diff_XShares=abs(gamma1dat-gamma1)+abs(AvgXSharedat-avgXshare);
        it=it+1;

    end
    
    RevFor=sum(revenuefor);
    XShareByFirm=Xshare;
    
    SkillIntNonExporters=SkillIntensityByFirm((NonExporters(:,1)==1));    
        
  % Moment 1. Average skill of non-exporters
    meanSkillIntNonExporters=mean(SkillIntNonExporters);

  % Moment 3. Average skill of exporters  
    meanSkillIntExporters=mean(SkillIntExporters);
    
    TotRevenue=sum(revenuedom+revenuefor,2);
    SizeDomestic=TotRevenue((NonExporters(:,1)==1));
    
    TotDomSalesActive=TotDomSales((ActiveDom(:,1)==1));    
    meanTotDomSales_nonexporter=mean(TotDomSales((NonExporters(:,1)==1)));
    meanTotDomSales_exporter=mean(TotDomSales((ActiveFor(:,1)==1)));
    

    % Calculating rest of the moments in model
    
    % 2. Regression of size and skill intensity 
    XVars_sk=[ones(length(SkillIntNonExporters),1),log(SizeDomestic)];
    Betas_Coefs_sk=XVars_sk\SkillIntNonExporters;
    delta1=Betas_Coefs_sk(2);
   
    % 4. Relative domestic sales of exporters
    RelativeDomSalesNonexpExp=meanTotDomSales_nonexporter/meanTotDomSales_exporter;
        
    % 7. normalized standard deviation of domestic sales
    mean_sd_sales=mean(TotDomSalesActive)/std(TotDomSalesActive);
    
    err1=abs(meanSkillIntNonExporters-AvgSkintNonExpdat);
    err2=wht_momelastsales*abs(delta1-delta1dat);
    err3=abs(meanSkillIntExporters-meanSkillIntExportersdat);
    err4=abs(RelativeDomSalesNonexpExp-RelativeDomSalesNonexpExpdat);
    err5=abs(mean_sd_sales-mean_sd_domsalesdat);
    if it==500
        penaliz=1;
    else
        penaliz=0;
    end
    error=err1+err2+err3+err4+err5+penaliz
   
  % Total moments
  SkTot=sum(sum(tot_skilled,2));
  UnskTot=sum(sum(tot_unskilled,2));
  ProfitDom=revenuedom/sigma;
  ProfitDom(ProfitDom==0)=NaN;
  [fl,posfl]=min(ProfitDom);
  ProfitFor=revenuefor/sigma;
  ProfitFor(ProfitFor==0)=NaN;
  [fx,posfx]=min(ProfitFor);
  fh=fl(2)+fx(2);
  
  omega=nnumfirms/NN;
  NetProfDom=omega*(ProfitDom-repmat(fl,length(ProfitDom),1));
  NetProfFor=omega*(ProfitFor-repmat(fx,length(ProfitFor),1));
  NetProfDom(isnan(NetProfDom)==1)=0;
  NetProfFor(isnan(NetProfFor)==1)=0;
  NetProfTot=sum(NetProfDom,2)+sum(NetProfFor,2);
  meanNetProfTot=mean(NetProfTot);
 
  
end









