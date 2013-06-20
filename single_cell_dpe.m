%% AMS 215 Project
% Modeling Cell Division Decisions
% by Raj Maitra
clc
clear all
close all

T = 100;
m_necro = 0.01; %probability of necrosis
mu = 0.5*(1:T)/T; %probability of mutation occurance
rho = 0.4*(1:T)/T;%(T:-1:1)/T; %probability that repair works
nrep = ceil(2*rand(200,1)); %random number of repairs done
nmut = ceil(4*rand(200,1));%0.5*(1:T); %random number of mutations done
p_div = 0.9;
xmax = 30;

%creating decision probabilities 
fapop = 0; %fitness from apoptosing
    amp_div = 5;
    phase = 15;
 
finvest = 0; %fitness from investing
    amp_invest = 15;
    mean = 10;
    width = 5;
    
fdiv = 0; %fitness from dividing
    amp_apop = 60;
    lambda = 1;
    shift = 5;
    
for x = 1:xmax
    fapop(x) = amp_div/(1+exp(-(x-phase)));
    finvest(x) = amp_invest*exp(-(x-mean)^2/(2*width^2));
    fdiv(x) = amp_apop*lambda*exp(-lambda*(x));
end
figure(1)
    subplot(221)
        plot(fapop,'-r')
        hold on 
        plot(finvest,'-b')
        plot(fdiv,'-m')
        title('Fitness Distributions for each Path')
        xlabel('Number of Mutations');
        ylabel('Fitness');
        legend('F_{apop}','F_{rep}','F_{div}')
        
    subplot(222)
        hist(nrep);
        hold on
        hist(nmut);

F = ones(xmax,T);
delta = zeros(xmax,T);

for t = T-1:-1:1
    for x = 1:xmax
        %probability of DNA or soma damage
        nmut = 5;
        nrep = 4;
        
        xdamage = max(x-nmut,1);
        
        [fdec_damage] = [fdiv(xdamage) + F(xdamage,t+1);...%p_div*(fdiv(xdamage) + F(xdamage,t+1)) + (1-p_div)*F(xdamage,t+1);...
                        rho(t)*(finvest(min((xdamage+nrep),xmax)) + F(min((xdamage+nrep),xmax),t+1))+(1-rho(t))*(finvest(xdamage) + F(min((xdamage),xmax),t+1));...
                        fapop(xdamage) + F(xdamage,t+1)];
                                  
        [fdecision] = [fdiv(x) + F(x,t+1);...%p_div*(fdiv(x) + F(x,t+1)) + (1-p_div)*F(x,t+1);...
                       rho(t)*(finvest(min((x+nrep),xmax)) + F(min((x+nrep),xmax),t+1))+(1-rho(t))*(finvest(x) + F(min((x),xmax),t+1));...
                       fapop(x) + F(x,t+1)];
                                  
        [Vmax delta(x,t)] = max(mu(t)*fdec_damage + (1-mu(t))*fdecision);
        F(x,t) = exp(-m_necro)*Vmax;
    end
end
F
    subplot(223);
    imagesc(delta(:,1:end-1))
    xlabel('Time')
    ylabel('Number of Mutations')
    title('Decision Boundaries')
    
    subplot(224);
    plot(F(:,1));
    if 0
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'finalfig2','-dpdf','-r0')
    end

    %% Forward Simulation
    pop = [];
    for j = 1:100
        j
    num_org = 1;
    x = zeros(num_org,T);
    x = x-1;    
    x(1,1) = 10;
    count = 0;
    t = 1;
    num_org = 1;
    
    
    for t = 1:T-1
        for k = 1:num_org
            if x(k,t) > 0
                r = rand; %generate random number (uniform)
                %if necrosis event
                if r <= m_necro
                    x(k,t+1) = -40;
                    continue
                end
                
                %if not a necrosis event
                %if mutations occur
                if rand <= mu(t)
                    nmut = ceil(1*rand);
                    x(k,t) = x(k,t) + nmut;
                end
                
                %if divide
                if delta(x(k,t),t) == 1
                    if count(k) >= 10
                        num_org = num_org + 1;
                        x = [x; zeros(1,T)-1];
                        x(num_org,t+1) = x(k,t);
                        count(k) = 0;
                        count = [count 0];
                        x(k,t+1) = x(k,t);
                        x(k,t) = -10;
                    else
                        count(k) = count(k)+1;
                        x(k,t+1) = x(k,t);
                    end
                    
                %if DNA repair
                elseif  delta(x(k,t),t) == 2
                    %if repair works
                    if rand <= rho(t)
                        nrep = ceil(10*rand);
                        x(k,t+1) = max(x(k,t)-nrep, 1);
                        count(k) = count(k)+1;
                        continue
                    end
                    x(k,t+1) = x(k,t);
                    count(k) = count(k)+1;
                %if apoptosis
                elseif  delta(x(k,t),t) == 3
                    x(k,t+1) = -20;
                    count(k) = count(k)+1;
                end
            end
        end
    end
    
    %exponential plot
    %figure(2)
    %imagesc(x)
    if 0
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'finalfig3','-dpdf','-r0')
    end
    
    %tree plot of offspring from one cell
    arrayi = x == -10;
    nodes = [0];
    for i = 1:T
        nodes = [nodes;find(arrayi(:,i) == 1)];
    end
    %figure(3) 
    %treeplot(nodes')
    if 0
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'finalfig5','-dpdf','-r0')
    end
    
    pop = [pop size(x,1)];
    end
    %%
    
    hold on
    distplot(pop','b')
    
    
    
    