lambda = 1550e-9;
k0 = 2 * pi / lambda;

waveguides = {
struct('name', 'Waveguide (a)', 'n_cl_TE', 1.5120, 'n_c_TE', 1.5375, 'n_s_TE', 1.4446, ...
'n_cl_TM', 1.5110, 'n_c_TM', 1.5370, 'n_s_TM', 1.4446),
struct('name', 'Waveguide (b)', 'n_cl_TE', 1.4446, 'n_c_TE', 2.1381, 'n_s_TE', 1.4446, ...
'n_cl_TM', 1.4446, 'n_c_TM', 2.2213, 'n_s_TM', 1.4446)
};

h_range = linspace(0.1e-6, 10e-6, 500);

f_TE = @(n_eff, h, m, n_c, n_cl, n_s, lambda) ...
( (2*pi/lambda) * sqrt(n_c^2 - n_eff.^2) * h - m*pi - ...
atan( (2*pi/lambda)*sqrt(n_eff.^2 - n_cl^2) ./ ((2*pi/lambda)*sqrt(n_c^2 - n_eff.^2)) ) - ...
atan( (2*pi/lambda)*sqrt(n_eff.^2 - n_s^2) ./ ((2*pi/lambda)*sqrt(n_c^2 - n_eff.^2)) ) );

f_TM = @(n_eff, h, m, n_c, n_cl, n_s, lambda) ...
( (2*pi/lambda) * sqrt(n_c^2 - n_eff.^2) * h - m*pi - ...
atan( (n_c^2/n_cl^2) * ((2*pi/lambda)*sqrt(n_eff.^2 - n_cl^2) ./ ((2*pi/lambda)*sqrt(n_c^2 - n_eff.^2)) ) ) - ...
atan( (n_c^2/n_s^2) * ((2*pi/lambda)*sqrt(n_eff.^2 - n_s^2) ./ ((2*pi/lambda)*sqrt(n_c^2 - n_eff.^2)) ) ) );

for w = 1:length(waveguides)
wg = waveguides{w};
disp(['=== Analyzing: ', wg.name, ' ===']);
n_eff_TE = nan(6, length(h_range));
n_eff_TM = nan(6, length(h_range));
for i = 1:length(h_range)
h = h_range(i);
for m = 0:5
f_mode_TE = @(n_eff) f_TE(n_eff, h, m, wg.n_c_TE, wg.n_cl_TE, wg.n_s_TE, lambda);
try
n_eff_TE(m+1, i) = fzero(f_mode_TE, [max(wg.n_cl_TE, wg.n_s_TE)+1e-6, wg.n_c_TE-1e-6]);
catch
n_eff_TE(m+1, i) = NaN;
end
f_mode_TM = @(n_eff) f_TM(n_eff, h, m, wg.n_c_TM, wg.n_cl_TM, wg.n_s_TM, lambda);
try
n_eff_TM(m+1, i) = fzero(f_mode_TM, [max(wg.n_cl_TM, wg.n_s_TM)+1e-6, wg.n_c_TM-1e-6]);
catch
n_eff_TM(m+1, i) = NaN;
end
end
end
figure;
hold on;
colors = lines(6);
for m = 0:5
valid_TE = ~isnan(n_eff_TE(m+1, :));
valid_TM = ~isnan(n_eff_TM(m+1, :));
plot(h_range(valid_TE)*1e6, n_eff_TE(m+1, valid_TE), '-', 'Color', colors(m+1, :), 'LineWidth', 1.5, 'DisplayName', ['TE, m=', num2str(m)]);
plot(h_range(valid_TM)*1e6, n_eff_TM(m+1, valid_TM), '--', 'Color', colors(m+1, :), 'LineWidth', 1.5, 'DisplayName', ['TM, m=', num2str(m)]);
end
hold off;
xlabel('Core Thickness h (\mum)');
ylabel('Effective Index n_{eff}');
title(['Dispersion Curves: ', wg.name]);
legend('Location', 'best');
grid on;
h_min_single_mode_TE = h_range(find(~isnan(n_eff_TE(1, :)), 1, 'first'));
h_max_single_mode_TE = h_range(find(~isnan(n_eff_TE(2, :)), 1, 'first'));
h_min_single_mode_TM = h_range(find(~isnan(n_eff_TM(1, :)), 1, 'first'));
h_max_single_mode_TM = h_range(find(~isnan(n_eff_TM(2, :)), 1, 'first'));
disp(['TE Single-mode propagation thickness range: ', num2str(h_min_single_mode_TE*1e6), ' µm ~ ', num2str(h_max_single_mode_TE*1e6), ' µm']);
disp(['TM Single-mode propagation thickness range: ', num2str(h_min_single_mode_TM*1e6), ' µm ~ ', num2str(h_max_single_mode_TM*1e6), ' µm']);
if ~isempty(h_max_single_mode_TE)
h_selected_TE = h_max_single_mode_TE;
else
h_selected_TE = h_range(end);
end
f_fund_TE = @(n_eff) f_TE(n_eff, h_selected_TE, 0, wg.n_c_TE, wg.n_cl_TE, wg.n_s_TE, lambda);
n_eff_selected_TE = fzero(f_fund_TE, [max(wg.n_cl_TE, wg.n_s_TE)+1e-6, wg.n_c_TE-1e-6]);
p_TE = k0 * sqrt(n_eff_selected_TE^2 - wg.n_s_TE^2);
q_TE = k0 * sqrt(n_eff_selected_TE^2 - wg.n_cl_TE^2);
min_substrate_TE = 1/p_TE;
min_cladding_TE = 1/q_TE;
if ~isempty(h_max_single_mode_TM)
h_selected_TM = h_max_single_mode_TM;
else
h_selected_TM = h_range(end);
end
f_fund_TM = @(n_eff) f_TM(n_eff, h_selected_TM, 0, wg.n_c_TM, wg.n_cl_TM, wg.n_s_TM, lambda);
n_eff_selected_TM = fzero(f_fund_TM, [max(wg.n_cl_TM, wg.n_s_TM)+1e-6, wg.n_c_TM-1e-6]);
p_TM = k0 * sqrt(n_eff_selected_TM^2 - wg.n_s_TM^2);
q_TM = k0 * sqrt(n_eff_selected_TM^2 - wg.n_cl_TM^2);
min_substrate_TM = 1/p_TM;
min_cladding_TM = 1/q_TM;
h_mid_TE = (h_min_single_mode_TE + h_max_single_mode_TE) / 2;
f_fund_mid_TE = @(n_eff) f_TE(n_eff, h_mid_TE, 0, wg.n_c_TE, wg.n_cl_TE, wg.n_s_TE, lambda);
n_eff_mid_TE = fzero(f_fund_mid_TE, [max(wg.n_cl_TE, wg.n_s_TE)+1e-6, wg.n_c_TE-1e-6]);
beta_mid_TE = k0 * n_eff_mid_TE;
h_mid_TM = (h_min_single_mode_TM + h_max_single_mode_TM) / 2;
f_fund_mid_TM = @(n_eff) f_TM(n_eff, h_mid_TM, 0, wg.n_c_TM, wg.n_cl_TM, wg.n_s_TM, lambda);
n_eff_mid_TM = fzero(f_fund_mid_TM, [max(wg.n_cl_TM, wg.n_s_TM)+1e-6, wg.n_c_TM-1e-6]);
beta_mid_TM = k0 * n_eff_mid_TM;
disp('--- Based on upper bound of single-mode propagation range ---');
disp(['TE: h_selected = ', num2str(h_selected_TE*1e6), ' µm, n_eff = ', num2str(n_eff_selected_TE)]);
disp(['TE Minimum substrate thickness (1/p): ', num2str(min_substrate_TE*1e6), ' µm']);
disp(['TE Minimum cladding thickness (1/q): ', num2str(min_cladding_TE*1e6), ' µm']);
disp(['TM: h_selected = ', num2str(h_selected_TM*1e6), ' µm, n_eff = ', num2str(n_eff_selected_TM)]);
disp(['TM Minimum substrate thickness (1/p): ', num2str(min_substrate_TM*1e6), ' µm']);
disp(['TM Minimum cladding thickness (1/q): ', num2str(min_cladding_TM*1e6), ' µm']);
disp('--- Based on midpoint of single-mode propagation range ---');
disp(['TE: h_mid = ', num2str(h_mid_TE*1e6), ' µm, n_eff_mid = ', num2str(n_eff_mid_TE), ', β = ', num2str(beta_mid_TE), ' rad/m']);
disp(['TM: h_mid = ', num2str(h_mid_TM*1e6), ' µm, n_eff_mid = ', num2str(n_eff_mid_TM), ', β = ', num2str(beta_mid_TM), ' rad/m']);
disp(' ');
end
