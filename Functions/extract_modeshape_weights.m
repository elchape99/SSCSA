function x_intersections = findIntersections(f, g, x_range, num_points)
    % Trova le intersezioni tra due funzioni f(x) e g(x) in un intervallo dato
    % f, g: function handles delle due funzioni
    % x_range: [x_min, x_max] intervallo di ricerca
    % num_points: numero di punti per la ricerca iniziale

    % Generazione di punti iniziali per la ricerca delle radici
    x_vals = linspace(x_range(1), x_range(2), num_points);
    y_vals = f(x_vals) - g(x_vals);
    
    % Trova gli intervalli in cui cambia segno (possibili intersezioni)
    sign_changes = find(y_vals(1:end-1) .* y_vals(2:end) < 0);
    
    x_intersections = [];
    
    % Usa fzero per trovare gli zeri con i punti iniziali trovati
    for i = 1:length(sign_changes)
        x_guess = x_vals(sign_changes(i));
        x_sol = fzero(@(x) f(x) - g(x), x_guess);
        x_intersections = [x_intersections, x_sol]; %#ok<AGROW>
    end
end
