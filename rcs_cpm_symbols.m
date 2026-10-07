function alpha = rcs_cpm_symbols(bits, params)

if ~isfield(params, 'cpm_M') || isempty(params.cpm_M)
    params.cpm_M = 4;
end
if ~isfield(params, 'cpm_mod_index_h') || isempty(params.cpm_mod_index_h)
    params.cpm_mod_index_h = 1/3;
end

M = params.cpm_M;
if M ~= 4
    error('rcs_cpm_symbols:M', '仅实现 M = 4（CC-CPM，表 7-23/7-24）');
end
m = log2(M);
bits = bits(:);
n_symbols = ceil(numel(bits) / m);
pad = n_symbols * m - numel(bits);
if pad > 0
    bits = [bits; zeros(pad, 1)];           % §7.3.7.2.2 末尾补 0 比特
end
bit_groups = reshape(bits, m, n_symbols).';
label = bit_groups(:, 1) * 2 + bit_groups(:, 2);      % MSB, LSB

if abs(params.cpm_mod_index_h - 1/3) < 1e-12
    tab = [-3, -1, 1, 3];                   % 表 7-23
else
    tab = [-3, -1, 3, 1];                   % 表 7-24（Gray）
end
alpha = tab(label + 1).';
alpha = alpha(:).';
end
