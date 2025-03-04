%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LUCAS ANDRADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x) -x * (1 - x.^2); % Definindo as fçs f, g e g'
g = @(x) 1 - x.^2;
g_prime = @(x) -2 * x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parâmetros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 0.5;  % VI
M=10; % numero de simulacoes
t0 = 0;    % t0
t_end = 10; % tfinal
dt_values = [0.1, 0.01]; % tamanhos de h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dt_idx = 1:length(dt_values)  % for para h = 0.1 e 0.01 
    dt = dt_values(dt_idx);

    figure; %para plotar em imagens diferentes 

    for sim = 1:M %% laço para o no de simulacoes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gerar movimento browniano W(t)
        N = floor((t_end - t0) / dt);
        W = sqrt(dt) * cumsum(randn(1, N+1));
        %%%%%%%%%%%%%%%%%%%%%%% Milstein %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [t, X_milstein] = milstein_method(f, g, g_prime, x0, t0, t_end, dt);
        %%%%%%%%%%%%%%%%%%%%%% Solucao exata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X_exact = ((1 + x0) * exp( W) + x0 - 1) ./ ((1 + x0) * exp(2 *W) + 1 - x0);
        %%%%%%%%%%%%%%%%%%%%% Plot do gráfico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(5, 2, sim); % Organizar os gráficos em um grid de 5x2
        plot(t, X_milstein, 'DisplayName', sprintf('Trajetória Milstein (Sim %d)', sim));
        hold on;
        plot(t, X_exact, 'DisplayName', 'Solução Exata');
        xlabel('Tempo');
        ylabel('X(t)');
        title(sprintf('Simulação %d com dt = %.2f', sim, dt));
        legend;
        grid on;
    end
end
