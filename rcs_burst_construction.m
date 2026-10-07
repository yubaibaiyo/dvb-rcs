function [burst_signal, burst_info] = rcs_burst_construction(symbols, params)

symbols = symbols(:);
n_data  = numel(symbols);

n_pre   = params.preamble_len;
n_post  = params.postamble_len;
n_blk   = params.pilot_block_len;
n_ps    = params.pilot_sum;
period  = params.pilot_period;

if n_data ~= params.payload_symbols
    warning('rcs_burst_construction:payload', ...
        ['净荷符号数 %d 与参考波形表规定的 %d 不一致（waveform_id=%d）：' ...
         '突发长度将偏离表值'], n_data, params.payload_symbols, params.waveform_id);
end

n_uw_sym = n_pre + (n_ps > 0) * n_blk + n_post;      % 无导频时 UW 里不含导频块
uw_sym   = uw_to_symbols(params.uw_hex, n_uw_sym, params);

preamble_sym = uw_sym(1:n_pre);

k = n_pre;
if n_ps > 0
    pilot_sym = uw_sym(k+1 : k+n_blk);
    k = k + n_blk;
else
    pilot_sym = zeros(0, 1);
end

postamble_sym = uw_sym(k+1 : k+n_post);



chunk = period - n_blk;                              % 每段净荷符号数
if n_ps > 0 && chunk <= 0
    error('rcs_burst_construction:pilot', ...
        'pilot_period=%d 不大于 pilot_block_len=%d，无法插入导频', period, n_blk);
end

burst_signal = preamble_sym;
n_pre_now = numel(burst_signal);
data_indices = zeros(0, 1);
pilot_indices = zeros(0, 1);
n_pilot_used = 0;
k = 1;                                               % 净荷游标

for b = 1:n_ps
    nn = min(chunk, n_data - k + 1);
    if nn <= 0
        break;                                       % 数据不足，不再插导频
    end
    burst_signal = [burst_signal; symbols(k:k+nn-1)];               
    data_indices = [data_indices; (numel(burst_signal)-nn+1 : numel(burst_signal)).'];  
    k = k + nn;
    if n_ps > 0
        burst_signal = [burst_signal; pilot_sym];                    
        pilot_indices = [pilot_indices; (numel(burst_signal)-n_blk+1 : numel(burst_signal)).'];  
        n_pilot_used = n_pilot_used + 1;
    end
end

if k <= n_data                                        % 末段剩余数据
    burst_signal = [burst_signal; symbols(k:end)];                   
    data_indices = [data_indices; (numel(burst_signal)-(n_data-k) : numel(burst_signal)).'];    
end

burst_signal = [burst_signal; postamble_sym];

if params.mod_order == 1
    % pi/2-BPSK: apply the outer rotation ONCE over the whole assembled burst so
    % that n is the absolute symbol index in the burst (7.3.7.1.4.1, table 7-16)
    n_idx = (0:numel(burst_signal)-1).';
    burst_signal = burst_signal .* exp(1j*(n_idx*pi/2 + pi/4));
end
n_post_now = numel(burst_signal);

%% ---------------- 位置索引与统计 ----------------
burst_info.preamble_idx  = (1:n_pre_now).';
burst_info.pilot_indices = pilot_indices;
burst_info.data_indices  = data_indices;
burst_info.postamble_idx = (n_post_now-n_post+1 : n_post_now).';
burst_info.n_preamble    = n_pre_now;
burst_info.n_pilots      = n_pilot_used;
burst_info.n_data        = n_data;
burst_info.n_postamble   = n_post;
burst_info.total_symbols = numel(burst_signal);
burst_info.preamble      = preamble_sym;
burst_info.pilot_block   = pilot_sym;
burst_info.postamble     = postamble_sym;

if numel(burst_signal) ~= params.burst_symbol_length
    warning('rcs_burst_construction:length', ...
        ['突发长度 %d 与参考波形表规定的 %d 不符（waveform_id=%d，' ...
         'pre=%d payload=%d pilot=%d x %d post=%d）'], ...
        numel(burst_signal), params.burst_symbol_length, params.waveform_id, ...
        n_pre, n_data, n_pilot_used, n_blk, n_post);
end

end


% ======================================================================
function sym = uw_to_symbols(uw_hex, n_sym, params)
m = params.mod_order;
bits = hex_to_bits(uw_hex);
n_bits = n_sym * m;
if numel(bits) < n_bits
    error('rcs_burst_construction:uw', ...
        ['UW 序列太短：uw_hex=%d 个十六进制字符 = %d 比特，需要 %d 比特' ...
         '（%d 个符号 x %d 比特）'], numel(uw_hex), numel(bits), n_bits, n_sym, m);
end
sym = rcs_modulation(bits(1:n_bits), params, false);
sym = sym(1:n_sym);
end

% ======================================================================
function bits = hex_to_bits(hex_str)
hex_str = upper(strtrim(char(hex_str)));
n = numel(hex_str);
bits = zeros(n * 4, 1);
for i = 1:n
    v = double(hex_str(i));
    if     v >= 48 && v <= 57,  d = v - 48;         
    elseif v >= 65 && v <= 70,  d = v - 55;          
    else
        error('rcs_burst_construction:hex', '十六进制串含非法字符 ''%s''', hex_str(i));
    end
    bits((i-1)*4+1:i*4) = bitget(d, 4:-1:1).';
end
end
