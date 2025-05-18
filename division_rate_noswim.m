function [mu_old] = division_rate_noswim(c_inf,cell_radius,Diat_noDiat)
%This function computes the growth rate.
%inputs are:
%c_inf - ambient nitrogent concentration (in nmole l^-1 = micro mole m-3) 
%cell_radius - cell radius in micron
%Diat_noDiat - equals 1 if cell is a diatom and 0 if it not a diatom.

% mu_old = per day

diam=2*cell_radius;
c_inf2=c_inf*1e-9*14; % c_inf in units of fg N/um^3
cell_volume=pi*4/3*cell_radius^3; % in micron^3


%case of non-diatom
if Diat_noDiat == 0
    C_per_cell=10^(-0.665+0.939*log10(cell_volume)); %from Menden-Deuer & Lessard, 2000 in pgC/cell
    D=0.000015; %diffusion coefficient for N in cm^2/s
    w=(93*diam^0.26)/10000; %swimming rate in cm/s
    Pe=w*cell_radius/D/10000; %Peclet #
    Sh=0.5*(1+(1+2*Pe)^(1/3)); %Sherwood # based on Karp et al., 1996
    Trans_eff=0.9; %Transporter efficiency
    Q=4*pi*D*10^8*cell_radius*c_inf2*Sh*Trans_eff; %max uptake potential fg N/s
    Q_daily=Q*3600*24;
    if diam<=11
        mu_max=-0.05552*diam^2+0.7889*diam+0.64028 ;%divisions per day 
    else
        mu_max=9.2869*diam^(-0.533);
    end
    R_NC_mumax=0.0762*mu_max+0.0389;
    C_percell=C_per_cell*1000;%fg C/cell
    V_max=C_percell*R_NC_mumax ;%fg N /day
    realized_uptake=V_max*Q_daily/(V_max+Q_daily);
    mu=mu_max*realized_uptake/V_max;
    mu_old=0;
    while abs(mu_old/mu-1)>0.0001
        R_NC_mu=0.0762*mu+0.0389;
        V_max=C_percell*R_NC_mu; %fg N /day
        realized_uptake=V_max*Q_daily/(V_max+Q_daily);
        mu=mu_max*realized_uptake/V_max; % this is per day, mu
        mu_old=mu;
    end
    
%case of diatom
elseif Diat_noDiat == 1
    C_per_cell=10^(-0.541+0.811*log10(cell_volume)); %from Menden-Deuer & Lessard, 2000 in pgC/cell
    D=0.000015; %diffusion coefficient for N in cm^2/s
    if cell_radius<=4
        w=0.0; %sinking/swimming rate
        Pe=0.0 ;%Peclet #
        Sh=1; %Sherwood # based on Karp et al., 1996
    else
        
        % COMMENT OR NOT
        w=0.0007*log(cell_radius/10000)+0.0056; % sinking rate
      %            w=(93*diam^0.26)/10000; %swimming rate in cm/s

        Pe=cell_radius*w/D/10000;
Sh=0.5*(1+(1+2*Pe)^(1/3));
    end
    Trans_eff=0.9; %Transporter efficiency
    Q=4*pi*D*10^8*cell_radius*c_inf2*Sh*Trans_eff; %max uptake potential fg N/s
    Q_daily=Q*3600*24;
    if diam<=11
        mu_max=-0.05552*diam^2+0.7889*diam+0.64028; %divisions per day 
    else
        mu_max=9.2869*diam^(-0.533);
    end
    R_NC_mumax=0.0762*mu_max+0.0389;
    C_percell=C_per_cell*1000; %fg C/cell
    V_max=C_percell*R_NC_mumax; %fg N /day
    realized_uptake=V_max*Q_daily/(V_max+Q_daily);
    mu=mu_max*realized_uptake/V_max;
    mu_old=0;
    while abs(mu_old/mu-1)>0.0001
        R_NC_mu=0.0762*mu+0.0389;
        V_max=C_percell*R_NC_mu; %fg N /day
        realized_uptake=V_max*Q_daily/(V_max+Q_daily);
        mu=mu_max*realized_uptake/V_max;
        mu_old=mu;
    end   
else
    disp('problem: Diat_noDiat should be 0 or 1')
    return
end
end

