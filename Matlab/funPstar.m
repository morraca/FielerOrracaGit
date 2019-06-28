% This is an auxiliary function to get the value of Pstar

function [value,Pstar_kh]=funPstar(Pstar_kl,import_k,alpha_kh,Pdom_kl,...
    Pdom_kh,Revenue_k,ssigma,TotExport_kl,TotExport_kh)
    
    alpha_kl=1-alpha_kh;
    Pstar_kh=Pdom_kh/((((TotExport_kl/TotExport_kh)*(alpha_kh/alpha_kl))^(1/(ssigma-1)))*(Pdom_kl/Pstar_kl));
    value=import_k-alpha_kl*Revenue_k*((Pdom_kl/Pstar_kl)^(ssigma-1))-...
        alpha_kh*Revenue_k*(Pdom_kh/Pstar_kh)^(ssigma-1);
    
end