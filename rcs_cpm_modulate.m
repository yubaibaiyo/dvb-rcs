function [symbols, cpm_info] = rcs_cpm_modulate(alpha, params)

if ~isfield(params, 'cpm_mod_index_h') || isempty(params.cpm_mod_index_h)
    params.cpm_mod_index_h = 1/3;
end
if ~isfield(params, 'cpm_alpha_rc') || isempty(params.cpm_alpha_rc)
    params.cpm_alpha_rc = 1;
end
if ~isfield(params, 'cpm_samples_per_sym') || isempty(params.cpm_samples_per_sym)
    params.cpm_samples_per_sym = 32;
end
if ~isfield(params, 'cpm_M') || isempty(params.cpm_M), params.cpm_M = 4;  end
if ~isfield(params, 'cpm_L') || isempty(params.cpm_L), params.cpm_L = 2;  end
if ~isfield(params, 'cpm_ph') || isempty(params.cpm_ph), params.cpm_ph = 3; end

h   = params.cpm_mod_index_h;
sps = params.cpm_samples_per_sym;
alpha_rc = params.cpm_alpha_rc;
alpha = alpha(:).';

% q(t)：0 (t<0)；RC t/4 - sin(pi t)/(4pi)、REC t/4 (0<=t<=2Ts)；t>2Ts 时 0.5
q_rc_func  = @(t) (t/4 - sin(pi * t) / (4 * pi)) .* (t >= 0 & t <= 2) + 0.5 * (t > 2);
q_rec_func = @(t) (t/4)                          .* (t >= 0 & t <= 2) + 0.5 * (t > 2);
q_func     = @(t) alpha_rc * q_rc_func(t) + (1 - alpha_rc) * q_rec_func(t);
g_rc_func  = @(t) (1/4) * (1 - cos(pi * t)) .* (t > 0 & t < 2);
g_rec_func = @(t) (1/4) * (t > 0 & t < 2);
g_func     = @(t) alpha_rc * g_rc_func(t) + (1 - alpha_rc) * g_rec_func(t);

phase = cpm_phase_trace(alpha, h, sps, q_func);
symbols = exp(1j * phase);

cpm_info.n_symbols       = numel(alpha);
cpm_info.alpha_values    = alpha;
cpm_info.phase_trace     = phase;
cpm_info.phase_start     = phase(1);            % 整突发时恒为 0（§7.3.7.2.1）
cpm_info.Vn_end          = mod(sum(alpha), params.cpm_ph);
cpm_info.h               = h;
cpm_info.M               = params.cpm_M;
cpm_info.L               = params.cpm_L;
cpm_info.mh              = mod(params.cpm_mod_index_h * params.cpm_ph, params.cpm_ph);
cpm_info.ph              = params.cpm_ph;
cpm_info.alpha_rc        = alpha_rc;
cpm_info.samples_per_sym = sps;
cpm_info.q_func          = q_func;
cpm_info.g_func          = g_func;
end


% ======================================================================
function phase = cpm_phase_trace(alpha, h, sps, q_func)
alpha = alpha(:).';
n = numel(alpha);
csum = [0, cumsum(alpha)];              % csum(m) = sum(alpha(1:m-1))
phase = zeros(n * sps, 1);
for s = 1:n
    for k = 1:sps
        t = (s - 1) + (k - 1) / sps;              % 绝对时间（单位：符号周期 Ts）
        j0 = max(1, ceil(t - 1));                 % t-(j-1) <= 2 的最早符号下标
        acc = 0;
        for j = j0:s
            acc = acc + alpha(j) * q_func(t - (j - 1));
        end
        acc = acc + 0.5 * csum(j0);               % j < j0 的符号：q = 0.5
        phase((s - 1) * sps + k) = 2 * pi * h * acc;
    end
end
end
