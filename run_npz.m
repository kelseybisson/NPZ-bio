% NPZ with recycle 
% run_npz.m

% Run simple NPZ model with 25 size classes, non-grazing phyto mortality, and 
% size structured feeding

% author: Kelsey Bisson, OSU, bissonk@oregonstate.edu

%           --> updated May 2025

clear all
close all

%% Case 1. non-diatoms, diatoms (MIXED case)

% 1. establish size classes growing roughly 25% each time
for i = 1:26
    psd(i,1) = 0.6*1.2471^(i-1); % phyto size diameter, um
end

% 2. calculate central bin size
for i = 1:25
    bins(i,1) = geomean([psd(i) psd(i+1)]); % um
    bwidth(i,1) = psd(i+1)-psd(i); %um
end

% 3. Set up model run details
dt=0.01;           % output timestep (days)
t0=0;              % tstep 1
tmax=365*3;        % n output timesteps, 3 years (nb: reaches equilibrium in 3 or less)
tspan=[t0:dt:tmax]; % all output timesteps

% set up nitrate levels with higher sampling resolution at lower levels 
niter = 20; N0iter(1:niter,1)= linspace(1,20,niter); 
N0iter(niter+1:80,1)= linspace(21,500,60); % nanomol/l = micro mole m-3
N0iter(81:100,1)= linspace(600,20000,20); % nanomol/l = micro mole m-3

% initial population numbers 
popinit = 0.18; % phyto mmol N m-3 
popinitz = 0.04; % zoo mmol N m-3 

% 4. loop through different Sinfinity levels
for i = 1:length(N0iter)
    
% params
g1 = 3.24; %  m3 mmol-1 d-1
g2 = 0.5; % unitless
g3 = 0.06;  %d-1
g4 = 1.6; % m3 mmol-1 d-1
fi =bins'; % feeding size range via proportionality factor 0.5
ji =bins'; % feeding size range via proportionality factor 0.5
mort = 0.01; % non grazing mortality, 1%

    % define phytoplankton sizes, ESD
     ESR = bins./2; % radius of cell
     N  = N0iter(i); % S inifinity in nM/L
  
    % calc growth rate (divisions) for diatoms, non-diatoms
    for j = 1:length(ESR)
    mu(1,j) = division_rate_noswim(N,ESR(j),0); % 0 = non diatoms
    mud(1,j) = division_rate_noswim(N,ESR(j),1); % 1 = diatoms
    end

    MU(i,:) = mu; MUD(i,:) = mud; 

    % loop over all time periods 
    for k = 1:length(tspan)
        if k ==1
            Pd = ones(1,25).*popinit; 
            Pnd = ones(1,25).*popinit; 
            P = ones(1,25).*popinit; 
            Z = ones(1,25).*popinitz;

        else      
        Pd = pd(k-1,:); % diatoms
        Pnd = pnd(k-1,:); % non diatoms
        P  = pnd(k-1,:) + pd(k-1,:); % mixed culture of both
        Z  = z(k-1,:); % zooplankton
        end       
      
        
     % calculate rate of change for diatoms (d) non diatoms and zoos
     dp = log(2).*mu.*Pnd -mort*log(2).*mu.*Pnd - fi.*g1.*Pnd.*Z.*log(2); 
     dpd = log(2).*mud.*Pd - mort*log(2).*mud.*Pd - fi.*g1.*Pd.*Z.*log(2); 
     dz = fi.*g1.*g2.*log(2).*(P).*Z-(g3.*Z)-(ji.*g4.*Z.^2);
    
    % calculate new biomass after 
    pnd(k,:) = Pnd + dp.*dt;
    pd(k,:) = Pd + dpd.*dt;
    z(k,:) = Z + dz.*dt;
   
    end
    
    % equlibrium values
    phyto(i,:) = pnd(end,:)+pd(end,:); % mixed assemblage
    dia(i,:) = pd(end,:);
    nondia(i,:) = pnd(end,:);
    zoo(i,:) = z(end,:);
    rec(i,:) = sum(fi.*g1.*(1-g2).*log(2).*(P).*Z + (g3.*Z) + (fi.*g4.*Z.^2))...
        + sum(mort*log(2).*mu.*Pnd + mort*log(2).*mu.*Pd); %recycle term
  
    i
end


% 5 calculate steady state carbon stocks, cell count, psd

phyC = (nondia)*106/16'; % non diatom from mg N to mg C
dC = (dia)*106/16'; % diatom from mg N to mg C

ccel = 10.^(-0.665+0.939*log10(4/3*pi*(bins./2).^3)); % pg/cell, non diatom
cceld = 10.^(-0.541+0.811*log10(4/3*pi*(bins./2).^3)); % pg/cell, diatom

for i = 1:25
cells(:,i) = [1e-9.*phyC(:,i).*1e12./ccel(i)]; % non diatom cells
cellsd(:,i) = [1e-9.*dC(:,i).*1e12./cceld(i)]; % diatom cells
end

% total mixed assemblage cells
Cells = cells+cellsd;

% calculate PSD
for i = 1:length(N0iter)
d= bins;
 
nbl1 = Cells(i,:)./bwidth';
[S,~] = polyfit(log10(d),log10(nbl1),1);
psd(i,1) = S(1); 
end

figure
plot(N0iter,-1*psd,'linewidth',2)
hold on
ylabel('SDS')
xlabel('S_{\infty}, nM/L')
