clc
close all
clear 
%% Script developed for use two piezo
% name convenction:
% - H: analitic function, computed with formulas
% - FRF: experimental measured valuesù
% - w_cap = w_oc_oc
% - w_i: freq naturali
% - w_1: freq ottime piezo 1
% - w_2: freq ottime piezo 2


% data of the bar 
beam_defm;

% load experimantal data
FRF_sc_sc_double_piezo = load ("FRF/FRF_sc_sc_double_piezo.mat");
FRF_oc_oc_double_piezo = load ("FRF/FRF_oc_oc_double_piezo.mat");
FRF_rl_rl_double_piezo_first = load ("FRF/FRF_rl_rl_double_piezo_first");

% rename correctly the variable
FRF_sc_sc = FRF_sc_sc_double_piezo.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl;
FRF_oc_oc = FRF_oc_oc_double_piezo.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl;
FRF_rl_rl = FRF_rl_rl_double_piezo_first.Data1_MT_FRF_H1_2Zplus_1Zplus_Ampl;


csi_i = beam.xi.sc_sc(1:2); % natural damping of the system
% with different damping change the peak, but the intersection point with
% the short circuit response does not change
csi1 = 0.005;  
csi2 = 0.005;  
k = [0.24, 0.18];
k1= k(1);
k2= k(2);

% freq used for the plot
freq = linspace(1,500, length(FRF_sc_sc)); 

% trasform freq in pulse
w = freq' * 2 * pi;

% w_i = beam.nf.sc_sc(1:2); % w_sc_sc
% w_i = [21.1297, 116.538];
% w_cap = beam.nf.oc_oc(1:2); % w_oc_oc
w_i = 2*pi* beam.nf.sc_sc_double_piezo; % w_sc_sc
%w_i = [21.1297*2*pi, 116.538 * 2*pi]; % trovate in matlab
w_cap = 2*pi* beam.nf.oc_oc(1:2); % w_oc_oc

%optimal values if consider two modes uncorrelated 
w_e_1_opt = (w_cap(1));
w_e_2_opt = (w_cap(2));

csi_e_1_opt = sqrt(3)/2 * sqrt ((w_cap(1)^2 - w_i(1)^2) / (w_cap(1)^2 + w_i(1)^2));
csi_e_2_opt = sqrt(3)/2 * sqrt ((w_cap(2)^2 - w_i(2)^2) / (w_cap(2)^2 + w_i(2)^2));

% computed with w_e_1 is a pulse not a period
L1_opt = 1 ./ (w_e_1_opt^2 * (beam.Cp.C11));
L2_opt = 1 ./ (w_e_2_opt^2 * (beam.Cp.C22));

% L in freq
% L1_opt = L1_opt/(2*pi)^2;
% L2_opt = L2_opt/(2*pi)^2;


% R1_opt_1 = 2 * csi_e_1_opt * sqrt (L1_opt / beam.Cp.C11);
% R2_opt_1 = 2 * csi_e_2_opt * sqrt (L2_opt / beam.Cp.C22);
R1_opt = (2 * csi_e_1_opt) / (beam.Cp.C11 * w_e_1_opt);
R2_opt = (2 * csi_e_2_opt) / (beam.Cp.C22 * w_e_2_opt);

% R in freq
% R1_opt = R1_opt * 2*pi;
% R2_opt = R2_opt * 2*pi;

% analitic H_sc_sc
H_sc_sc = 1./(w_i.^2 + 1i.*2.*csi_i .*w_i.*w - w.^2);
H_oc_oc = 1./(w_cap.^2 + 1i.*2.*csi_i .*w_cap.*w - w.^2);
% experimental measure of velocity force : moltiplication of all H(w)
% analitic by jw
H_sc_sc = 1i * w .* H_sc_sc;
H_oc_oc = 1i * w .* H_oc_oc;

% analitic H_rl_rl
H_rl_rl = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, ...
                                  [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                                  L1_opt, L2_opt, R1_opt, R2_opt, ...
                                  beam.k.k1(1:2), beam.k.k2(1:2), ones(1,2));      
% experimental measure of velocity7force : moltiplication of all H(w)
% analitic by jw
H_rl_rl = 1i * w .* H_rl_rl;

%% Curve fitting

% === Funzione modello della FRF ===
H_model = @(phi_squared, x) abs(phi_squared(1) .* H_sc_sc(:,1) + ...
                                phi_squared(2) .* H_sc_sc(:,2));

% === Curve fitting sulla magnitudo della FRF ===
phi_init = [0.1, 0.1]; % Stima iniziale dei pesi modali
phi_opt = lsqcurvefit(H_model, phi_init, double(freq), double(abs(FRF_sc_sc)));

H_sc_sc_fitted = abs(phi_opt(1) * H_sc_sc(:,1) + ...
                     phi_opt(2) * H_sc_sc(:,2));

% Plot results 
disp('Mode shapes relativi nel punto di misura:');
disp(sqrt(phi_opt));


% H_rl_rl_fitted = abs(phi_opt(1) * H_rl_rl(:,1) + ...
%                      phi_opt(2) * H_rl_rl(:,2));

H_rl_rl_fitted = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, ...
                                        [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                                        L1_opt, L2_opt, R1_opt, R2_opt, ...
                                        beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);      
% experimental measure of velocity force : moltiplication of all H(w)
% analitic by jw
H_rl_rl_fitted = abs(sum(1i .* w .* H_rl_rl_fitted,2));


% %% optimization of L and R values
% tot = 5;
% csi_i_optimization = zeros(1,2);
% 


%% H_rl_rl, L optimization 
% In this section we otpimize the H_rl_rl
% We will use the fitted H_rl_rl function, so with implemented the mode
% shape inside

% Variabili per memorizzare il miglior risultato
tot = 2;  % Intervallo di ricerca
counter = 1;
csi_i_optimization = zeros(1,2); % for optimization put mechanical damping == 0

best_param1 = NaN;
best_param2 = NaN;
min_error = Inf;

% Loop sui due parametri
% uguale ordinata di intersezione tra rl-rl e sc-sc
% for ii = (w_e_1_opt - tot) : 0.01 : (w_e_1_opt + tot)
%     for jj = (w_e_2_opt - tot) : 0.01 : (w_e_2_opt + tot)
%         % computed with w_e_1 is a pulse not a period
%         Lii = 1 ./ (ii^2 * (beam.Cp.C11));
%         Ljj = 1 ./ (jj^2 * (beam.Cp.C22));
% 
%         % Rii = 0; %(2 * csi_e_1_opt) / (beam.Cp.C11 * ii);
%         % Rjj = 0; %(2 * csi_e_2_opt) / (beam.Cp.C22 * jj);
%         Rii = R1_opt; 
%         Rjj = R2_opt;
% 
%         func_sc_sc = @(w) abs(sum(1 ./(w_i.^2 - w.^2), 2));
%         func_rl_rl = @(w) abs(sum(double_piezo_reson_FRF (w, w_i, w_cap, csi_i_optimization, [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], Lii, Ljj, Rii, Rjj, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt), 2));
%         x_intersection = findIntersections(func_sc_sc, func_rl_rl, w);
%         if (length(x_intersection) < 6)
%             continue;
%         end
%         locs = x_intersection([1, 2, 5, 6]); %points 3 and 4 are in anti-resonance
%         pks = arrayfun (func_sc_sc, locs);
%         total_error = abs(pks(1) - pks(2)) + abs(pks(3) - pks(4));
% 
% 
%         if total_error < min_error
%             min_error = total_error;
%             best_param1 = ii;
%             best_param2 = jj;
%             % best_FRF_values = FRF_values;  % Salva la funzione ottimale
%             best_FRF_values = arrayfun (func_rl_rl, w);  % Salva la funzione ottimale
%             best_locs = locs(1:4);
%             best_pks = pks(1:4);
%             % hystory(:,counter) = FRF_values;
%             % counter = counter +1;
%         end
%     end
% end

% uguale altezza dei picchi (coppia per coppia)
for ii = (w_e_1_opt - tot) : 0.1 : (w_e_1_opt + tot)
    for jj = (w_e_2_opt - tot) : 0.1 : (w_e_2_opt + tot)
% for ii = (w_e_1_opt - tot) : (w_e_1_opt + tot)
%     for jj = (w_e_2_opt - tot) : (w_e_2_opt + tot)
        % Calcolo dei parametri L e R
        Lii = 1 ./ (ii^2 * (beam.Cp.C11));
        Ljj = 1 ./ (jj^2 * (beam.Cp.C22));

        Rii = R1_opt; 
        Rjj = R2_opt; 

        % Valutazione della funzione su tutto il range di frequenze
        FRF_values = arrayfun(@(w) abs(sum(double_piezo_reson_FRF(w, w_i, w_cap, csi_i_optimization, ...
                              [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                              Lii, Ljj, Rii, Rjj, beam.k.k1(1:2), beam.k.k2(1:2), phi_opt))), w);

        % Trova i primi 4 picchi
        % [pks, locs] = findpeaks(FRF_values, w, 'SortStr', 'descend');
        [pks, locs] = findpeaks(FRF_values, w);

        % Controlla se ci sono almeno 4 picchi
        if length(pks) < 4
            continue;
        end

        % Seleziona i primi 4 picchi
        pks = pks(1:4);

        % Calcola l'errore (differenza tra i primi due e tra il terzo e quarto)
        error1 = abs(pks(1) - pks(2));  % Differenza tra picco 1 e 2
        error2 = abs(pks(3) - pks(4));  % Differenza tra picco 3 e 4
        total_error = error1 + error2;  % Errore totale da minimizzare

        % Aggiorna i migliori parametri se l'errore è minore
        if total_error < min_error
            min_error = total_error;
            best_param1 = ii;
            best_param2 = jj;
            best_FRF_values = FRF_values;  % Salva la funzione ottimale
            best_locs = locs(1:4);
            best_pks = pks(1:4);
            hystory(:,counter) = FRF_values;
            counter = counter +1;
        end
    end
end

% Output dei risultati
fprintf('Best parameters found:\n');
fprintf('w_e1 = %.3f, w_e2 = %.3f, error = %.3f\n', best_param1 / (2*pi), best_param2 / (2*pi), min_error);

% From w_e1, w_e2 compute the optima values for L1, L2
L1_opt_new = 1 ./ (best_param1^2 * (beam.Cp.C11));
L2_opt_new = 1 ./ (best_param2^2 * (beam.Cp.C22));


%% damped with previous optimal resistances h_rl_rl computation
H_rl_rl_damped = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, ...
                                        [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                                        L1_opt_new, L2_opt_new, R1_opt, R2_opt, ...
                                        beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);     
H_rl_rl_damped = sum(1i .* w .* H_rl_rl_damped, 2);

%% R optimization
tot = 0.1; 
tot2 = 0.05;
min_error = inf;
csi1_cand = (csi_e_1_opt - tot): 0.1 : (csi_e_1_opt + tot);
csi2_cand = (csi_e_2_opt - tot) : 0.1 : (csi_e_2_opt + tot);


for ii = unique(max(zeros(1, length(csi1_cand)), csi1_cand))
    for jj = unique (max(zeros(1, length(csi2_cand)), csi2_cand))
% for ii = 0:0.1:(sqrt(3)/3)
%     for jj = 0:0.1:(sqrt(3)/3)

        Rii = 2 * ii * sqrt (L1_opt_new / beam.Cp.C11); 
        Rjj = 2 * jj * sqrt (L2_opt_new / beam.Cp.C22);
        
        % Valutazione della funzione su tutto il range di frequenze
        FRF_values_damped = arrayfun(@(w) abs(sum(double_piezo_reson_FRF(w, w_i, w_cap, csi_i_optimization, ...
                              [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                              L1_opt_new, L2_opt_new, Rii, Rjj, ...
                              beam.k.k1(1:2), beam.k.k2(1:2), phi_opt))), w);

        % Trova i primi 4 picchi
        [pks, locs] = findpeaks(FRF_values_damped, w);
        
        % Controlla se ci sono almeno 4 picchi
        if length(pks) < 4
            continue;
        end
        
        % Seleziona i primi 4 picchi
        pks = pks(1:4);

        % Calcola l'errore (differenza tra i primi due e tra il terzo e quarto)
        error1 = abs(pks(1) + pks(2));  % Differenza tra picco 1 e 2
        error2 = abs(pks(3) + pks(4));  % Differenza tra picco 3 e 4
        total_error = error1 + error2;  % Errore totale da minimizzare

        % Aggiorna i migliori parametri se l'errore è minore
        if total_error < min_error
            min_error = total_error;
            best_param1 = ii;
            best_param2 = jj;
            best_FRF_value_damped = FRF_values_damped;  % Salva la funzione ottimale
            best_locs = locs(1:4);
            best_pks = pks(1:4);
            % hystory(:,counter) = FRF_values;
            % counter = counter +1;
        end
    end
end
%% damped h_rl_rl with optimal resistances computation

R1_opt_new = 2 * best_param1 * sqrt (L1_opt_new / beam.Cp.C11);
R2_opt_new = 2 * best_param2 * sqrt (L2_opt_new / beam.Cp.C22);

H_rl_rl_damped_opt = double_piezo_reson_FRF (w, w_i, w_cap, csi_i, ...
                                        [beam.Cp.C11, beam.Cp.C12], [beam.Cp.C21, beam.Cp.C22], ...
                                        L1_opt_new, L2_opt_new, R1_opt_new, R2_opt_new, ...
                                        beam.k.k1(1:2), beam.k.k2(1:2), phi_opt);     
H_rl_rl_damped_opt = sum(1i .* w .* H_rl_rl_damped_opt, 2);