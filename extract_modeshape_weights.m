function [mode_shapes, FRF_weighted] = extract_modeshape_weights(FRF_exp, freq_exp, FRF_ana_single, freq_ana)
    % Input:
    % - FRF_exp: Matrice delle FRF sperimentali (sensori x frequenze)
    % - freq_exp: Vettore frequenze FRF sperimentali
    % - FRF_ana_single: Matrice delle FRF analitiche single mode (modi x frequenze)
    % - freq_ana: Vettore frequenze FRF analitiche
    % 
    % Output:
    % - mode_shapes: Mode shapes sperimentali estratti
    % - FRF_weighted: FRF analitica combinata con i pesi sperimentali

    % 1. Estrazione dei mode shapes con SVD
    [U, S, ~] = svd(FRF_exp, 'econ'); % Decomposizione SVD della FRF sperimentale
    mode_shapes = U; % I mode shapes sono nelle colonne di U
    
    % Normalizzazione mode shapes
    mode_shapes = mode_shapes ./ vecnorm(mode_shapes, 2, 1);
    
    % 2. Proiezione delle FRF analitiche single mode sui mode shapes sperimentali
    % Calcoliamo i coefficienti pesanti i singoli modi analitici
    weights = mode_shapes' * FRF_exp;  % Pesi basati sulla proiezione
    
    % 3. Costruzione della FRF analitica combinata
    FRF_weighted = FRF_ana_single' * weights; % Ricostruzione con pesi

    % % 4. Plot della FRF ricostruita rispetto alla sperimentale
    % figure;
    % subplot(2,1,1);
    % semilogy (freq_exp, abs(FRF_exp), "LineWidth",0.6);
    % %imagesc(freq_exp, 1:size(FRF_exp,1), abs(FRF_exp)); 
    % title('FRF Sperimentale');
    % xlabel('Frequenza [Hz]'); ylabel('Sensori');
    % colorbar;
    % 
    % subplot(2,1,2);
    % semilogy (freq_exp, abs(FRF_weighted), "LineWidth",0.6);
    % % imagesc(freq_exp, 1:size(FRF_weighted,1), abs(FRF_weighted)); 
    % title('FRF Ricostruita con Mode Shapes');
    % xlabel('Frequenza [Hz]'); ylabel('Sensori');
    % colorbar;
    % 
    % % Output della funzione
    % disp('Mode Shapes Sperimentali estratti e FRF analitica combinata calcolata.');
end
