function [ generation_hist ] = pso_alg( fit_func, dim, var_min, var_max, pop_size, max_gen )
%PSO_ALG optimization algorithm
generation_hist = zeros(1, max_gen);
disp('pso optimizing...');
%% initialization:
% velocity limitations:
w=1;              % Inertia Weight
wdamp = 0.99;     % Inertia Weight Damping Ratio
vel_max = repmat(0.1*(var_max - var_min), 1, dim);
vel_min = -vel_max;
C1 = 2;
C2 = 2;
gbest.sol = [];
gbest.fit = Inf;
particle.sol = [];
particle.pbest.sol = [];
particle.pbest.fit = Inf;
particle.vel = [];
pso_particles = repmat(particle, 1, pop_size);
for i=1:pop_size
    pso_particles(i).sol = (var_max - var_min) .* rand(1, dim) + var_min;
    pso_particles(i).vel = zeros(1, dim);
    pso_particles(i).pbest.sol = [];
    pso_particles(i).pbest.fit = Inf;
end

%% main loop:
for itr=1:max_gen
    for i=1:pop_size % for each particle
        par_fit = fit_func(pso_particles(i).sol);
        % check and update pbest:
        if par_fit < pso_particles(i).pbest.fit
            pso_particles(i).pbest.fit = par_fit;
            pso_particles(i).pbest.sol = pso_particles(i).sol;
        end
        % check and update gbest:
        if par_fit < gbest.fit
            gbest.fit = par_fit;
            gbest.sol = pso_particles(i).sol;
        end
    end
    
    for i=1:pop_size
        % update vel:
        pso_particles(i).vel = w*pso_particles(i).vel + ...
            C1 * rand(1, dim) .* (pso_particles(i).pbest.sol - pso_particles(i).sol) + ...
            C2 * rand(1, dim) .* (gbest.sol - pso_particles(i).sol);
        % check velocity:
        pso_particles(i).vel = max(min(pso_particles(i).vel, vel_max), vel_min);
        % velocity mirror effect:
        is_out = (pso_particles(i).sol < var_min | pso_particles(i).sol > var_max);
        pso_particles(i).vel(is_out) = -pso_particles(i).vel(is_out);
        % update solution:
        pso_particles(i).sol = pso_particles(i).sol + pso_particles(i).vel;
        % check solution limits:
        pso_particles(i).sol = max(min(pso_particles(i).sol, repmat(var_max, 1, dim)), ...
            repmat(var_min, 1, dim));
    end
    %log:
    generation_hist(itr) = gbest.fit;
    % update discount:
    w = w*wdamp;
end

end

