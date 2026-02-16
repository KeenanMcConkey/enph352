% ODE function for the rate equations
% 16x1 column vector (Ni, Nk)
function dNdt = op_ode(t, N)
    dNdt = zeros(16, 1);
    
    % Intensity I is global
    global I;       % Intesity [mW/cm^2]
    
    % Relevant parameters from live script
    hBar = 1.054572e-34;                    % Reduced Planck's constant [J*s]
    c0   = 2.997925e8;                      % Speed of light [m/s]
    e0   = 8.85418782e-12;                  % Permitivity of free space [F/m]
    tau = 27.68e-9;                         % Lifetime of the excited manifold [s]
    omega0 = 2 * pi * 377.10746338e12;      % Resonant frequency of the transition [Hz]

    % Oscillator strength matrix
    fik = [ 1/12 1/12 0 1/2 1/4 1/12 0 0 ;
        1/12 0 1/12 0 1/4 1/3 1/4 0  ;    
        0 1/12 1/12 0 0 1/12 1/4 1/2 ;
        1/2 0 0 1/3 1/6 0 0 0        ;
        1/4 1/4 0 1/6 1/12 1/4 0 0   ;
        1/12 1/3 1/12 0 1/4 0 1/4 0  ;
        0 1/4 1/4 0 0 1/4 1/12 1/6   ;
        0 0 1/2 0 0 0 1/6 1/3        ];
  
    % LCP selection matrix
    Cik = [ 0 1 0 0 0 1 0 0 ;
        0 0 1 0 0 0 1 0 ;
        0 0 0 0 0 0 0 1 ;
        1 0 0 0 1 0 0 0 ;
        0 1 0 0 0 1 0 0 ;
        0 0 1 0 0 0 1 0 ;
        0 0 0 0 0 0 0 1 ;
        0 0 0 0 0 0 0 0 ];    
    
    % dNi/dt equations
    for i = 1:8
        for k = 9:16
            kk = k-8; % Actual k value
            dNdt(i) = dNdt(i) + ...
                      6 * pi * c0^2 / hBar / omega0^3 * 10*I * fik(i, kk) * Cik(i, kk) * (N(k) - N(i)) - ...
                      1 / tau * fik(i, kk) * N(i);
        end
    end
    
    % dNk/dt equations
    for k = 9:16
        kk = k-8; % Actual k value
        for i = 1:8
            dNdt(k) = dNdt(k) + ...
                      6 * pi * c0^2 / hBar / omega0^3 * 10*I * fik(i, kk) * Cik(i, kk) * (N(i) - N(k)) + ...
                      1 / tau * fik(i, kk) * N(i);
        end
    end
end