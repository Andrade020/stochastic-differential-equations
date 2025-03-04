%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LUCAS ANDRADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% nessa versao eu repiito a mesma trajetória para a solução Exata- 
%%% Verificar Questao2_b_alternativa !! lá eu apenas uso o mesmo incremento
f = @(x) -x * (1 - x.^2); % Definindo as fçs f, g e g'
g = @(x) 1 - x.^2;
g_prime = @(x) -2 * x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parâmetros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M= 10; % numero de simulacoes
x0 = 0.5;  % V_inic
t0 = 0;    % t0
t_end = 10; % tfinal
dt_fino = 0.01; % Passo de tempo mais fino  % tamanhos de h
dt_grosso = 0.1; % Passo de tempo mais grosso  % tamanhos de h

% gera browniano W com os dois h's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_fino = floor((t_end - t0) / dt_fino);
W_fino = sqrt(dt_fino) * cumsum(randn(1, N_fino + 1));
indices_grosso = 1:(dt_grosso/dt_fino):length(W_fino);
W_grosso = W_fino(indices_grosso);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  solução exata    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_exact_fino = ((1 + x0) * exp(W_fino) + x0 - 1) ./ ((1 + x0) * exp(2 *W_fino) + 1 - x0);
X_exact_grosso = ((1 + x0) * exp(W_grosso) + x0 - 1) ./ ((1 + x0) * exp(2 *W_grosso) + 1 - x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para plotar separado: 
for dt = [dt_grosso, dt_fino]
    figure;

     %separa o browniano%%%%%%%%%%%%%%%%%%%%
    if dt == dt_fino
        W = W_fino;
        X_exact = X_exact_fino;
    else
        W = W_grosso;
        X_exact = X_exact_grosso;
    end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M simulacoes (padrao M=10 )
    for sim = 1:M
        % Chamando a função milstein_method
        [t, X_milstein] = milstein_method(f, g, g_prime, x0, t0, t_end, dt);

        %%%%%%%%%%%%%%%%%%%%% Gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(5, 2, sim); % Organizar os gráficos em um grid de 5x2
        plot(t, X_milstein, 'DisplayName', sprintf('Trajetória Milstein (Sim %d)', sim));
        hold on;
        plot(t, X_exact, 'DisplayName', 'Solução Exata');
        xlabel('Tempo');
        ylabel('X(t)');
        title(sprintf('Simulação %d com dt = %.2f', sim, dt));
        legend;
        grid on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end
