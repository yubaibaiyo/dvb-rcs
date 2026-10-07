function [symbols, cpm_info] = rcs_cpm_modulation(bits, params, terminate)

if ~isfield(params, 'cpm_M') || isempty(params.cpm_M), params.cpm_M = 4; end
if ~isfield(params, 'cpm_L') || isempty(params.cpm_L), params.cpm_L = 2; end
if ~isfield(params, 'cpm_ph') || isempty(params.cpm_ph), params.cpm_ph = 3; end
if ~isfield(params, 'cpm_mh') || isempty(params.cpm_mh)
    params.cpm_mh = 1;
end
if ~isfield(params, 'cpm_mod_index_h') || isempty(params.cpm_mod_index_h)
    params.cpm_mod_index_h = 1/3;
end
if ~isfield(params, 'cpm_samples_per_sym') || isempty(params.cpm_samples_per_sym)
    params.cpm_samples_per_sym = 32;
end
if nargin < 3 || isempty(terminate)
    terminate = false;
end

M  = params.cpm_M;
ph = params.cpm_ph;
sps = params.cpm_samples_per_sym;

%% 比特 -> 符号（表 7-23 / 7-24），并记录 2 比特标签
bits = bits(:);
m = log2(M);
n_symbols = ceil(numel(bits) / m);
pad = n_symbols * m - numel(bits);
if pad > 0
    bits = [bits; zeros(pad, 1)];
end
bit_groups = reshape(bits, m, n_symbols).';
label = bit_groups(:, 1) * 2 + bit_groups(:, 2);

alpha_data = rcs_cpm_symbols(bits, params);

%% 尾符号（表 7-25~7-28）：本段从累加器 0 起算
Vn = mod(sum(alpha_data), ph);
use_term = terminate;
if use_term
    [term_values, n_term] = rcs_cpm_termination(Vn, ph);
    alpha_all = [alpha_data, term_values];
else
    term_values = zeros(0, 1);
    n_term = 0;
    alpha_all = alpha_data;
end

%% 相位累加（数据 + 尾符号连续）
[symbols, mi] = rcs_cpm_modulate(alpha_all, params);
n_term_used = n_term;
n_symbols_total = numel(alpha_all);

if n_term_used > 0
    term_phase   = mi.phase_trace(n_symbols * sps + 1:end);
    term_symbols = exp(1j * term_phase);
else
    term_phase   = [];
    term_symbols = [];
end

%% 输出
cpm_info.n_symbols       = n_symbols;
cpm_info.n_symbols_total = n_symbols_total;
cpm_info.n_term          = n_term_used;
cpm_info.n_bits          = numel(bits);
cpm_info.data_labels     = label(:).';
cpm_info.term_values     = term_values;
cpm_info.term_alpha      = term_values;
cpm_info.term_phase      = term_phase;
cpm_info.term_symbols    = term_symbols;
cpm_info.phase_trace     = mi.phase_trace;
cpm_info.alpha_values    = alpha_all;
cpm_info.Vn              = Vn;
cpm_info.Vn_end          = mi.Vn_end;
cpm_info.mapping         = sprintf('h=%d/%d: %s', params.cpm_mh, ph, ...
                            ternary(abs(params.cpm_mod_index_h - 1/3) < 1e-12, ...
                                    'table 7-23', 'table 7-24 (Gray)'));
cpm_info.h               = params.cpm_mod_index_h;
cpm_info.M               = M;
cpm_info.L               = params.cpm_L;
cpm_info.mh              = params.cpm_mh;
cpm_info.ph              = ph;
cpm_info.alpha_rc        = mi.alpha_rc;
cpm_info.samples_per_sym = sps;
cpm_info.q_func          = mi.q_func;
cpm_info.g_func          = mi.g_func;
end


% ======================================================================
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
