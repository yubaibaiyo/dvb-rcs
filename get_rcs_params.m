function params = get_rcs_params(burst_type, coding_scheme, waveform_cfg)

if nargin < 3 || isempty(waveform_cfg)
    waveform_cfg = struct();
end
if ~isfield(waveform_cfg, 'waveform_id') || isempty(waveform_cfg.waveform_id)
    waveform_cfg.waveform_id = 1;
end

%% ---------------- 编码方案 -> 参考波形表 ----------------
params.burst_type    = upper(burst_type);
params.coding_scheme = lower(coding_scheme);
switch params.coding_scheme
    case 'turbo'
        params.mod_type = 'linear';
        fmt = 'A1';
    case 'concatenated'
        params.mod_type = 'cpm';
        fmt = 'A3';
    otherwise
        error('get_rcs_params:scheme', 'coding_scheme 只能为 concatenated 或 turbo');
end

%% ---------------- 查内置参考波形表 ----------------
T   = rcs_waveform_table(fmt);
wid = waveform_cfg.waveform_id;
k   = find([T.id] == wid, 1);
if isempty(k)
    error('get_rcs_params:waveform', '表 %s 无 waveform_id = %d（可用: %s）', fmt, wid, mat2str([T.id]));
end
w = T(k);
params.waveform_id = wid;

%% ---------------- 码率：直接查表所得 ----------------
params.code_rate_str = w.code_rate;              % 表里的码率字符串原样使用
rp = sscanf(params.code_rate_str, '%d/%d');      % 仅把表里的 'n/m' 拆成分子/分母
params.code_rate_num = rp(1);
params.code_rate_den = rp(2);
params.code_rate     = rp(1) / rp(2);

%% ---------------- 逐表取参数 ----------------
if strcmp(params.mod_type, 'cpm')
    params.mod_order           = [];             % CC-CPM 无线性调制阶数
    params.mod_name            = 'CC-CPM (M=4)';
    params.mapping             = w.cc_type;
    params.cpm_M               = w.M;
    params.cpm_L               = w.L;
    params.cpm_mh              = w.mh;
    params.cpm_ph              = w.ph;
    params.cpm_mod_index_h     = w.h;
    params.cpm_alpha_rc        = w.alpha_rc;
    params.cpm_rolloff         = 0.4;
    params.cpm_samples_per_sym = 32;
    params.cpm_data1_len       = w.data1_bits;
    params.cpm_data2_len       = w.data2_bits;
    params.cpm_term_bits       = w.term_bits;
    params.cpm_term_symbols    = w.term_symbols;
    params.uw_hex              = w.uw_hex;
    params.preamble_len        = 64;
    params.cpm_midamble_len    = 64;
    params.pilot_period        = 0;
    params.pilot_block_len     = 0;
    params.pilot_sum           = 0;
    params.payload_symbols     = (w.data1_bits + w.data2_bits) / log2(w.M);
    params.burst_symbol_length = w.burst_symbols;
    params.payload_bytes       = [];             % 表 A-3 给的是比特数，无字节列
else
    params.mod_order           = w.mod_order;    % 直接查表
    params.mod_name            = w.mapping;
    params.mapping             = w.mapping;
    params.payload_bytes       = w.payload_bytes;% 表 A-1 的净荷字节数（= 原始数据长度/8）
    params.preamble_len        = w.preamble_len;
    params.postamble_len       = w.postamble_len;
    params.pilot_period        = w.pilot_period;
    params.pilot_block_len     = w.pilot_block;
    params.pilot_sum           = w.pilot_sum;
    params.payload_symbols     = w.payload_symbols;
    params.burst_symbol_length = w.burst_symbols;
    params.uw_hex              = w.uw_hex;
end

%% ---------------- 允许逐项覆盖（仅突发构造相关）----------------
ovr = {'preamble_len', 'postamble_len', 'pilot_period', 'pilot_block_len', ...
       'pilot_sum', 'payload_symbols', 'uw_hex'};
for i = 1:numel(ovr)
    f = ovr{i};
    if isfield(waveform_cfg, f) && ~isempty(waveform_cfg.(f))
        params.(f) = waveform_cfg.(f);
    end
end

%% ---------------- CRC 类型----------------
switch params.burst_type
    case 'TRF'
        params.crc_type = 'crc32';
    otherwise
        params.crc_type = 'crc16';
end

%% ---------------- 编码器输入长度（由参考波形表的净荷长度反推）----------------
% 突发不带任何头部：编码器输入 = 原始数据 + CRC
%   A-1（TC-LM）给的是净荷字节数  -> fec_input_bits = payload_bytes*8
%   A-3（CC-CPM）直接给比特数      -> fec_input_bits = 表中的 fec_input_bits
if strcmp(params.mod_type, 'cpm')
    params.fec_input_bits = w.fec_input_bits;
else
    params.fec_input_bits = w.payload_bytes * 8;
end
switch params.crc_type
    case 'crc32'
        params.crc_bits = 32;
    otherwise
        params.crc_bits = 16;
end
params.input_len = params.fec_input_bits - params.crc_bits;

%% ---------------- 卷积码（CC-CPM，7.3.5.2）----------------
% 母码率固定 1/2；K∈{3,4} 与生成多项式 (5,7)o / (15,17)o 一一绑定（查表）
params.conv_mother_rate = 1/2;
params.conv_K    = 3;
params.conv_g1   = 5;      % 1 + x^2（八进制字面值）
params.conv_g2   = 7;      % 1 + x + x^2
params.conv_cc_type = '(5,7)o';
if strcmp(params.mod_type, 'cpm') && w.K == 4
    params.conv_K  = 4;
    params.conv_g1 = 15;   % 1 + x + x^3
    params.conv_g2 = 17;   % 1 + x + x^2 + x^3
    params.conv_cc_type = '(15,17)o';
end
params.conv_tail_bits = params.conv_K - 1;     % 2 或 3（§7.3.5.2.2）
params.conv_params.g           = [params.conv_g1, params.conv_g2];
params.conv_params.g1          = params.conv_g1;
params.conv_params.g2          = params.conv_g2;
params.conv_params.K           = params.conv_K;
params.conv_params.mother_rate = 1/2;
params.conv_params.target_rate = params.code_rate_str;   % 直接来自参考波形表

%% ---------------- Turbo（TC-LM，7.3.5.1）----------------
params.turbo_K        = 5;              % 4 级寄存器
params.turbo_n_states = 16;
params.turbo_polys    = {23 35 27};     % 反馈 23o、Y 35o、W 27o
params.turbo_params.poly        = params.turbo_polys;
params.turbo_params.n_states    = params.turbo_n_states;
params.turbo_params.target_rate = params.code_rate_str;  % 直接来自参考波形表
params.turbo_params.N_inner     = [];
%% ---------------- TC-LM turbo 交织/打孔参数（表 A-1 行内列，§7.3.5.1.1 / §7.3.5.1.3）----------------
if strcmp(params.mod_type, 'linear')
    params.turbo_params.N         = w.payload_bytes * 4;   % couple 数 N = 净荷比特/2
    params.turbo_params.P         = w.P;
    params.turbo_params.Q0        = w.Q0;
    params.turbo_params.Q1        = w.Q1;
    params.turbo_params.Q2        = w.Q2;
    params.turbo_params.Q3        = w.Q3;
    params.turbo_params.y_period  = w.y_period;
    params.turbo_params.y_pattern = w.y_pattern;
    params.turbo_params.w_period  = w.w_period;
    params.turbo_params.w_pattern = w.w_pattern;
end

%% ---------------- MODCOD 组合校验（§7.3.7.1.4）----------------
if ~isempty(params.mod_order) && params.mod_order == 3 && params.code_rate < 2/3
    error('get_rcs_params:modcod', '8PSK 要求码率 ≥ 2/3');
end
if ~isempty(params.mod_order) && params.mod_order == 4 && params.code_rate < 3/4
    error('get_rcs_params:modcod', '16QAM 要求码率 ≥ 3/4');
end
end
