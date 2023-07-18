% b10
% - experimental b9; uses real values
%% setup

traj    = 0;                        % # of trajectories to save
dt      = 6*10^(-4);                % Euler method time step

global  phi_m
global  v0

H0      = 0.7;
v0      = 3*H0^2;

% construction parameters
phi_m   = 1;                        % mesa start
Dphi_e  = 0.0;                      % exit displacement down slope

% NOTE: for eps = 0.01 crossover occurs at 0.047
Dphi_uv = 0.005;                    % allowed backtracking
Dphi_m  = [0.043 + 0.001*(0:12)];

% diffusion regime
Dphi_m  = [0.00005 0.00006 0.00007 0.00008 0.00009 0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001];

% sample size
N       = [1*10^6 * ones(1,length(Dphi_m))];
Pdt     = 0.005;                    % PDF bar width
dphi    = 10^(-5);                  % H block width

% slow-roll parameters
eps     = 0.00;
eta     = 0;

% momentum and slingshot    
M_in    = -sqrt(2*eps*v0/3);        % <!> deprecated by Dphi_0 approx.
Dphi_0  = sqrt(2*eps)/3;

%% data acquisition: recycle, load or simulate as needed

if exist('P','var') == 1 && input('Recycle data? [Y/N]: ', 's') == 'Y'
    % do nothing
elseif input('Load PDF from file? [Y/N]: ', 's') == 'Y'
    
    prefix  = input('File prefix: ','s');
    P       = {};
    i       = 1;
    while isfile([prefix '_' num2str(i) '.dat'])

        P{i}    = readmatrix([prefix '_' num2str(i) '.dat']);
        P{i}    = P{i}/trapz(P{i});     % normalise if needed
        i       = i+1;
    end
    
    disp(P)

else
    
    disp('Running sim...')
    tic
    
    P   = {};     	% PDF tally cell array
    X   = {};       % nested cell arrays for trajectory data
    T   = {};
    
    for m = 1:length(Dphi_m)
        
        % start/exit points
        phi_e   = phi_m - Dphi_e; 
        phi_in  = phi_m + Dphi_m(m);
        phi_uv  = phi_in + Dphi_uv;
        phi_0   = phi_in - Dphi_0;
        
        %% construct H/M arrays

        H       = [H0 + eps*10^(-10)];  % H starting from phi_0
        phi     = phi_0;
        
        while phi < phi_uv              % evaluate up slope -> uv
            
            v_          = v( phi, phi_in, eps, eta );
            % H from HJ eq.
            H_          = H(length(H));
            H           = [H H_ + dphi*sqrt(3/2 * H_^2 - 0.5*v_)];
            phi         = phi + dphi;
        end
        
        if phi_0 > phi_m            % IF non-graceful exit...

            phi         = phi_0;

            while phi > phi_m       % traverse remaining LS; flood with H0

                H       = [H0 H];
                phi     = phi - dphi;
            end
        else                        % IF graceful exit...
            
            % discard H for phi < phi_m
            del         = length(phi_0:dphi:phi_m)-1;
            H(1:del)    = [];
            
            phi         = phi_m;
        end

        M_ln    = [];
        M       = [];               % UNUSED (display only)
        for i = 1:length(H)-1       % get field momentum over all H
            
            M(i)        = -2/dphi * ( H(i+1) - H(i) );      
            M_ln(i)     = -2/dphi * ( log(H(i+1)) - log(H(i)) );
        end
        
        %% H/M verification (debug only)
        %{
        for i = 1:length(H)-1       % get field momentum over all H
            
            dH          = H(i+1) - H(i);
            phi         = phi_e + (i-1)*dphi;
            
            LHS         = (dH/dphi)^2;
            RHS         = 3/2 * H(i)^2 - 1/2 * v(phi, Dphi_m(m));
            err         = LHS - RHS;
            
            disp(err)
            
        end
        %}
        
        %% H/M plot (updated tuesday 05 may)
        %{
        figure
        hold on
        xlabel('\phi')
        
        colororder({'#0072BD','#EDB120'})
        
        Phi     = phi_e:dphi:(phi_uv+dphi);
        Phi_an  = phi_0:dphi:(phi_uv+dphi)
        
        yyaxis left
        plot(Phi, H,'LineWidth',1.2)
        plot(Phi_an, H0*cosh(sqrt(3/2) * 1.02*(Phi_an-phi_0)), '--','LineWidth',1.0)
        ylabel('H(\phi)')
        yticks([H0 max(H)])
        yticklabels({'H_0' 'H_{max}'})
        
        yyaxis right
        plot(phi_e:dphi:phi_uv, M,'LineWidth',1.2)
        ylabel('\Pi(\phi)')
        yticks([M_in 0])
        yticklabels({'\Pi_{in}' '0'})
        
        xline(phi_e, 'k:', 'LineWidth', 1.0);
        xline(phi_in, 'k:', 'LineWidth', 1.0);
        xline(phi_0, 'r:', 'LineWidth', 1.2);
        
        text(phi_e,0.01,'\leftarrow\phi_e')
        text(phi_in,0.01,'\leftarrow\phi_{in}'); return
        %}
        
        
        %% do random walk
        
        P{m}    = [];                   % PDF tally
        NH      = length(H) - 1;        % cache block count
        S       = normrnd(0,1,1,10^5);  % preallocate Gaussian noise

        n       = 0;                    % counter
        X{m}    = [];                   % trajectory data
        T{m}    = [];
        
        %%% <!> DEBUG
        disp(['======== ' num2str(m) ' ========'])
        disp(['D        = ' num2str(Dphi_m(m))])
        disp(['Start    = ' num2str(phi_in)])
        disp(['phi_0    = ' num2str(phi_0)])
        disp(['End      = ' num2str(phi_e)])
        if phi_0 > phi_e
            disp(['Non-graceful exit. Fraction spent in diffusion: ' num2str((phi_0-phi_e)/(phi_in-phi_e))])
        else
            disp(['Graceful exit! Excess fraction: ' num2str((phi_in-phi_0)/(phi_in-phi_e))])
        end
        disp('.')
        %%%

        while n < N(m)

            n   = n + 1;                % move to next sample

            phi     = phi_in;           % inject phi
            r       = [phi_in];
            steps   = 0;
            % do random walk
            while phi > phi_e
                
                k           = S(ceil(rand*10^5));                               % fetch stochastic term
                i           = 1 + floor( NH * (phi-phi_e)/(phi_uv-phi_e) );     % get occupied block
                jump        = M_ln(i)*dt + H(i)/(2*pi) *k*dt;                   % phi jump (field momentum + QD)
                phi         = phi + jump;                                       % move phi
                
                %%% DEBUG
                % disp(['drift/diff: ' num2str(abs((M_ln(i))/(k*H(i)/(2*pi))))])
                
                if phi > phi_uv         % no backtracking
                    phi = phi - jump;
                end

                steps   = steps + 1;    % count step
                
                if n <= traj            % update trajectory
                    r   = [r phi];
                end
                
            end
            
            % sim done -> find right PDF block for exit time
            i       = ceil(steps*dt/Pdt);
            if i > length(P{m})
                P{m}(i) = 1;                % extend tally if it's too short
            else
                P{m}(i) = P{m}(i) + 1;      % add one to tally
            end
            
            % counter
            if mod(n,10^4) == 0
                disp([num2str(Dphi_m(m)) ': ' num2str(n)])
            end
            
            if n <= traj                    % append trajectory
                X{m}{n} = r;
                T{m}{n} = dt*(0:steps);
            end
            
        end
        
        %P{m} = P{m}/trapz(P{m});        % normalise PDF
    end
    
    toc
    disp('Sim complete!')
    
    if input('Save data? [Y/N]: ','s') == 'Y'
        
        prefix = input('File prefix: ','s');
        for i = 1:length(P)
            
            writematrix( P{i}, [prefix '_' num2str(i) '.dat'] )
        end
    end
    
end

%% trajectory plot
%{
if traj > 0
    
    figure
    hold on
    xlabel('Time')
    ylabel('\phi')
    
    for n = 1:traj
        plot(T{end}{n}, X{end}{n}, 'HandleVisibility','off')
    end
    
    % add lines
    yline(phi_uv,'k:');
    yline(phi_e,'k:');
    yline(phi_0,'r:', 'LineWidth', 1.2);
    
    if phi_0 < phi_e
        yticks([phi_0 phi_e phi_in phi_uv])
        yticklabels({'\phi_0' '\phi_e' '\phi_{in}' '\phi_{uv}'})
    else
        yticks([phi_e phi_0 phi_in phi_uv])
        yticklabels({'\phi_e' '\phi_0' '\phi_{in}' '\phi_{uv}'})
    end
end
%}

%% PDF

figure
hold on
xlabel('$\mathcal{N}_e$','interpreter','latex')
ylabel('$P(\mathcal{N}_e)$','interpreter','latex')

% custom colormap (plugin)
map = customcolormap([0 1], {'#e02232', '#226be0'});
colorbar
colormap( map )

for m = 1:length(Dphi_m)
    
    t_e             = Pdt * (1:length(P{m}));           % exit time array
    P{m}            = P{m}/trapz(t_e,P{m});             % normalise

    m_i             = ceil(256 * m/length(Dphi_m));     % colormap index
    pm              = plot(t_e, P{m}, 'Color', map(m_i,:));
    pm.Color(4)     = 0.2 + 0.8*(m/length(Dphi_m));
end

caxis([Dphi_m(1) Dphi_m(length(Dphi_m))])
ylabel(colorbar,'\Delta\phi_m')


%% cumulants/normality #1

K           = {};
for m = 1:length(Dphi_m)
    
    mom     = [0 0 0 0];
    %{
    for i = 1:length(P{m})          % evaluate moments
        for j = 1:4
            mom(j)   = mom(j) + (Pdt*(i-1))^j * P{m}(i);
            
        end
    end
    %}
    
    t_e             = Pdt * (1:length(P{m}));
    for j = 1:4
        mom(j)      = trapz(t_e, t_e.^j .* P{m});
    end
    
    K{1}(m) = mom(1);               % evaluate cumulants
    K{2}(m) = mom(2) - mom(1)^2;
    K{3}(m) = mom(3) - 3*mom(2)*mom(1) + 2*mom(1)^3;
    K{4}(m) = mom(4) - 4*mom(3)*mom(1) - 3*mom(2)^2 + 12*mom(2)*mom(1)^2 - 6*mom(1)^4;
end

%%% OVERWRITE
%K{1}(13)    = 8.86;
%K{2}(13)    = 45.33;
%%%

% cumulant plot
figure
hold on
grid on
plot(Dphi_m, K{1}, 'o-', 'LineWidth', 1.1)
plot(Dphi_m, K{2}, 'o-', 'LineWidth', 1.1)
axis manual
plot(Dphi_m, K{3}, 'o-', 'LineWidth', 1.1)
plot(Dphi_m, K{4}, 'o-', 'LineWidth', 1.1)

xline(Dphi_0,'r:','LineWidth',1.2);

legend('\mu','\sigma^2','\kappa_3','\kappa_4','Location','northwest')
xlabel('\Delta\phi_s')
ylabel('\kappa_i')


%% normality #2: SW
%{
SW          = [];
SW_row      = [];
for k = 1:1000
    for m = 1:length(Dphi_m)
        PX          = floor(N(m)*P{m});
        X           = [];

        % take a small sample
        for i = 1:40

            r       = floor(N(m)*rand());
            tot     = 0;
            j       = 0;
            while tot < r
                
                j   = j + 1;
                tot = tot + PX(j);
            end
            X       = [X j*Pdt];
        end
        
        % perform SW test
        [H, pValue, W]  = swtest(X);
        SW_row          = [SW_row H];
    end
    
    SW      = [SW; SW_row];
    SW_row  = [];
end

% pass fraction plot
pratio      = [];
var         = [];

for m = 1:length(Dphi_m)
    
    X       = transpose( SW(:,m) );     % take data column and transpose
    passes  = length(X(X==0));          % number of passes
    pratio  = [pratio passes/1000];     % append pass ratio
    
    mean    = sum(X)/1000;              % mean/binary variance (http://capone.mtsu.edu/dwalsh/VBOUND2.pdf)
    var     = [var passes*(1000-passes)/1000^2];
    
end

%%
figure
errorbar(Dphi_m, pratio, var, '-s', 'LineWidth', 1.1)
xlabel('\Delta\phi_s')
ylabel('Shapiro-Wilk pass fraction')
axis([-inf inf 0 1])

xline(Dphi_0,'r:','LineWidth',1.2);
%}
%% PBH mass fraction

beta = [];

figure  % temp
hold on
xlabel('\zeta_{cg}')
ylabel('P(\zeta_{cg})')
xline(1,':','LineWidth',1.2);

% calculate beta
for m = 1:length(Dphi_m)
    
    % translate PDF -> zeta_cg
    z_e     = Pdt * (1:length(P{m})) - K{1}(m);     % equiv to t_e - <N> (?)
    
    % trim
    z_f     = [z_e(z_e>1) 1];
    P_f     = [P{m}(z_e>1) 0];
    
    % colormap/plot
    m_i             = ceil(256 * m/length(Dphi_m));
    C               = map(m_i,:);
    plot(z_e, P{m}, 'Color', C);
    F               = fill(z_f, P_f, C, 'EdgeColor', 'none');
    F.FaceAlpha     = 0.1;
    
    % integrate over tail for beta
    beta(m) = trapz(z_f,P_f)
    
end

% save beta data (for multi-beta)
if input('Save beta? [Y/N]: ','s') == 'Y'

    prefix = input('Save prefix: ','s');
    writematrix( Dphi_m, ['beta_' prefix '_x.dat'] )
    writematrix( beta, ['beta_' prefix '_y.dat'] )
end

%%
% potential
function v = v(phi, phi_in, eps, eta)
    
    global v0
    
    if phi < phi_in
        
        v = v0;
    else
        
        v = v0*(1 + sqrt(2*eps)*(phi-phi_in) + eta/2*(phi-phi_in)^2);
    end
    
end