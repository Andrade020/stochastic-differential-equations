%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LUCAS ANDRADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% inspirada no vídeo: https://www.youtube.com/watch?v=Nx2_btCEfrE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, X] = milstein_method(f, g, g_prime, x0, t0, t_end, dt)
    N = floor((t_end - t0) / dt);  % no de passos
    t = linspace(t0, t_end, N+1);
    X = zeros(1, N+1);
    X(1) = x0;

    for i = 2:N+1
        dW = sqrt(dt) * randn;  % Incremento Browniano 
        X(i) = X(i-1) + f(X(i-1)) * dt + g(X(i-1)) * dW + 0.5 * g(X(i-1)) * g_prime(X(i-1)) * (dW.^2 - dt);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% no matlab surpreendentemente devemos definir as funcoes potumamente %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINIÇÕES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = f(x)
    val = -x * (1 - x.^2);
end

function val = g(x)
    val = 1 - x.^2;
end

function val = g_prime(x)
    val = -2 * x;
end
