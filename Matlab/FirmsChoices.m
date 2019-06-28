%param=[sigma,rho,eta,varsigma];
% This program does step 7: Calculate costs, skill intensity, price for
% each firm.

function [MgCost,pricedom,pricefor,skint,auxu]=...
    FirmsChoices(vvarphis,bbetas,vvarsigma,fix_par,ttau,wwu,wws,NN)

    eeta=fix_par(3);
    rrho=fix_par(2);
    auxMgCost=(bbetas.*((wws*(vvarphis.^(-vvarsigma))).^(1-eeta))+ (ones(NN,2)-bbetas)*(wwu^(1-eeta))).^(1/(1-eeta));
    MgCost=(ones(NN,2)./vvarphis).*auxMgCost;
    pricedom=MgCost./repmat(rrho,NN,2);
    pricefor=ttau*pricedom;
    skint=(bbetas./(1-bbetas))*((wwu/wws)^(eeta)).*(vvarphis.^(vvarsigma*(eeta-1)));
    auxu=vvarphis.*( ( ((ones(NN,2)-bbetas).^(1/eeta)) .* ((bbetas./(1-bbetas)).*(((wwu/wws)^(eeta-1))*(vvarphis.^(vvarsigma*(eeta-1)))) +1)  ).^(eeta/(eeta-1)));
    
end