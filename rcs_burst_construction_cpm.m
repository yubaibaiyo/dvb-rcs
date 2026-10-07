function [burst_alpha, burst_info] = rcs_burst_construction_cpm(data_bits, params)
data_bits = data_bits(:);

ph = params.cpm_ph;
M  = params.cpm_M;

%% UW 符号（前导 = 中导）
uw_bits = hex_to_bits(params.uw_hex);
uw_sym  = rcs_cpm_symbols(uw_bits, params).';

%% 数据分两段（按表 A-3 的 Data #1 bit length）
n1 = params.cpm_data1_len;
if n1 > numel(data_bits), n1 = numel(data_bits); end
d1 = data_bits(1:n1);
d2 = data_bits(n1+1:end);
sym1 = rcs_cpm_symbols(d1, params).';
sym2 = rcs_cpm_symbols(d2, params).';

%% 逐段组装（累加器跨段连续，终止后归零）
burst_alpha = uw_sym;
acc = mod(sum(uw_sym), ph);

burst_alpha = [burst_alpha; sym1];
acc = mod(acc + sum(sym1), ph);
[tail1, n_term1] = rcs_cpm_termination(acc, ph);
Vn1 = mod(acc, ph);
burst_alpha = [burst_alpha; tail1(:)];
acc = mod(acc + sum(tail1), ph);                 % 终止后应为 0
Vn_after_term1 = acc;

burst_alpha = [burst_alpha; uw_sym];
acc = mod(acc + sum(uw_sym), ph);

tail2 = []; n_term2 = 0; Vn2 = NaN; Vn_after_term2 = NaN;
if ~isempty(sym2)
    burst_alpha = [burst_alpha; sym2];
    acc = mod(acc + sum(sym2), ph);
    [tail2, n_term2] = rcs_cpm_termination(acc, ph);
    Vn2 = mod(acc, ph);
    burst_alpha = [burst_alpha; tail2(:)];
    acc = mod(acc + sum(tail2), ph);
    Vn_after_term2 = acc;
end

%% 统计
burst_info.n_uw_symbols      = numel(uw_sym);
burst_info.n_preamble_bits   = numel(uw_bits);
burst_info.n_midamble_bits   = numel(uw_bits);
burst_info.n_data_bits       = numel(data_bits);
burst_info.n_data1_symbols   = numel(sym1);
burst_info.n_data2_symbols   = numel(sym2);
burst_info.tail1             = tail1;
burst_info.tail2             = tail2;
burst_info.n_term1           = n_term1;
burst_info.n_term2           = n_term2;
burst_info.Vn1               = Vn1;
burst_info.Vn2               = Vn2;
burst_info.Vn_after_term1    = Vn_after_term1;
burst_info.Vn_after_term2    = Vn_after_term2;
burst_info.total_symbols     = numel(burst_alpha);
burst_info.total_bits        = numel(burst_alpha) * log2(M);

end


% ======================================================================
function bits = hex_to_bits(hex_str)
%HEX_TO_BITS 每 4 比特一个十六进制字符，高位在前
hex_str = upper(strtrim(char(hex_str)));
n = numel(hex_str);
bits = zeros(n * 4, 1);
for i = 1:n
    v = double(hex_str(i));
    if     v >= 48 && v <= 57,  d = v - 48;          % '0'-'9'
    elseif v >= 65 && v <= 70,  d = v - 55;          % 'A'-'F'
    else
        error('rcs_burst_construction_cpm:hex', '十六进制串含非法字符 ''%s''', hex_str(i));
    end
    bits((i-1)*4+1:i*4) = bitget(d, 4:-1:1).';
end
end
