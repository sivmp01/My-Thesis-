%% verify_case3.m  -- run after the solver finishes
%clear; clc;

% Load results dumped by your model 
load('Case3_TCVaR_new.mat');  % provides: n, Ns, pi_s, delta_t, B, cap_price, pen_short,
                         % t_h, t_h_n, R_opt, S1v, Pel1v, S2_all, Pel2_all, Price_all,
                         % sshort_all, T_min, T_max, f_idx

% Do prices differ across scenarios?
mP = mean(Price_all,1);% average price per scenario
dP12 = max(abs(Price_all(:,1)-Price_all(:,2)));
dP13 = max(abs(Price_all(:,1)-Price_all(:,3)));
dP23 = max(abs(Price_all(:,2)-Price_all(:,3)));
fprintf('Mean price per scenario: %g  %g  %g\n', mP);
fprintf('Max ΔPrice across scen: [1-2]=%g  [1-3]=%g  [2-3]=%g\n', dP12, dP13, dP23);

% Do realized powers differ across scenarios?
mPel2 = mean(Pel2_all,1);
dQ12 = max(abs(Pel2_all(:,1)-Pel2_all(:,2)));
dQ13 = max(abs(Pel2_all(:,1)-Pel2_all(:,3)));
dQ23 = max(abs(Pel2_all(:,2)-Pel2_all(:,3)));
fprintf('Mean Pel2 per scenario: %g  %g  %g\n', mPel2);
fprintf('Max ΔPel2 across scen:  [1-2]=%g  [1-3]=%g  [2-3]=%g\n', dQ12, dQ13, dQ23);

% Is Stage-2 basically matching Stage-1?
dBase = max(max(abs(Pel2_all - Pel1v)));
fprintf('Max |Pel2 - Pel1| across all scen/hours: %g\n', dBase);

% Re-state key constants used in the model 
Res_min   = 1000;   % kW
Res_max   = 3000;   % kW  (set to your chosen max)
S_max3    = 150;    % kWh
alpha     = 0.8;    % same alpha as in the model
tol       = 1e-6;

% 1) Dimensions & basic sanity
assert(isvector(R_opt) && numel(R_opt)==n, 'R_opt size mismatch');
assert(all(size(Pel2_all)==[n Ns]), 'Pel2_all size mismatch');
assert(all(size(Price_all)==[n Ns]), 'Price_all size mismatch');

% 2) Reserve bounds
viol_upper = max( max(0, R_opt - min(Res_max, B*Pel1v) - 1e-5) );
viol_min   = max( max(0, (R_opt(R_opt>tol) < (Res_min - 1e-5))) );
fprintf('Reserve upper bound worst viol (kW): %.3g\n', viol_upper);
fprintf('Reserve min (1 MW) ok hours: %d / %d\n', sum(R_opt<=tol | R_opt>=Res_min-1e-5), n);

% 3) Realized power never exceeds commitment (per scenario)
viol_pel = max(0, max(max(Pel2_all - (Pel1v + 1e-5))));
fprintf('Pel2 <= Pel1 worst viol (kW): %.3g\n', viol_pel);

% 4) Shortfall logic: sshort >= 0 AND sshort >= Res_RA - B*(Pel1 - Pel2)
lhs = sshort_all;                                           % n x Ns
rhs = R_opt - B*(Pel1v - Pel2_all);                         % n x Ns (implicit expansion)
viol_short_nonneg = max(0, -min(lhs(:)));
viol_short_logic  = max(0, max(max(rhs - lhs - 1e-5)));
fprintf('sshort >= 0 worst viol: %.3g\n', viol_short_nonneg);
fprintf('sshort >= Res - B*(Pel1-Pel2) worst viol (kW): %.3g\n', viol_short_logic);

% 5) SoC limits
viol_S1_lo = max(0, max(0 - S1v));
viol_S1_hi = max(0, max(S1v - S_max3));
viol_S2_lo = max(0, max(0 - S2_all, [], 'all'));
viol_S2_hi = max(0, max(S2_all - S_max3, [], 'all'));
fprintf('S1 in [0,Smax] worst viol (lo/hi): %.3g / %.3g\n', viol_S1_lo, viol_S1_hi);
fprintf('S2 in [0,Smax] worst viol (lo/hi): %.3g / %.3g\n', viol_S2_lo, viol_S2_hi);

% Overall PASS/FAIL
ok = (viol_upper<1e-5) && all(R_opt<=Res_max+1e-5) ...
     && (viol_min==0) && (viol_pel<1e-5) ...
     && (viol_short_nonneg<1e-5) && (viol_short_logic<1e-5) ...
     && (viol_S1_lo<1e-5) && (viol_S1_hi<1e-5) && (viol_S2_lo<1e-5) && (viol_S2_hi<1e-5);
if ok
    fprintf('\n ALL CHECKS PASS.\n');
else
    fprintf('\n Some checks failed, mentioned below\n');
end

%Reconstruct scenario costs & CVaR summary (energy + shortfall penalty) 
% Energy cost per scenario: sum_t [ B * Pel2(t,w) * delta_t * Price(t,w) ]
E_cost = zeros(Ns,1);
PEN    = zeros(Ns,1);
for w = 1:Ns
    E_cost(w) = sum( B * Pel2_all(:,w) .* delta_t .* Price_all(:,w) );
    PEN(w)    = pen_short * sum( sshort_all(:,w) .* delta_t / 1000 );  % €/MWh * (kW*h)/1000
end
SC = E_cost + PEN;  % scenario costs (no capacity revenue here)

% Compute VaR_alpha and CVaR_alpha from (SC, pi_s)
[SCs, idx] = sort(SC, 'ascend');
pis = pi_s(:); pis = pis(idx);
cum = cumsum(pis);
VaR = SCs(find(cum >= alpha, 1, 'first'));
% Discrete CVaR (include all scenarios >= VaR)
mask = (SCs >= VaR - 1e-12);
CVaR = sum(SCs(mask).*pis(mask)) / sum(pis(mask));

fprintf('\n--- TCVaR summary (energy+penalty only) ---\n');
disp(table((1:Ns)', SC, 'VariableNames', {'Scenario','Cost_EUR'}));
fprintf('VaR_%.2f (EUR) = %.2f\n', alpha, VaR);
fprintf('CVaR_%.2f (EUR) = %.2f\n', alpha, CVaR);
