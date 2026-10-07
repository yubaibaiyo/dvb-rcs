function turbo_out = rcs_turbo_encode(data, turbo_params)

data = data(:);
K = numel(data);
if mod(K, 2) ~= 0
    error('rcs_turbo_encode:K', 'FEC 输入比特数 K=%d 不是偶数（K = 2N）', K);
end
N = K / 2;
if isfield(turbo_params, 'N') && ~isempty(turbo_params.N)
    N = turbo_params.N;
    if 2*N ~= K
        error('rcs_turbo_encode:N', 'turbo_params.N=%d 与 K/2=%d 不符', N, K/2);
    end
end
if mod(N, 4) ~= 0
    error('rcs_turbo_encode:N4', 'N=%d 必须是 4 的倍数（§7.3.5.1.0）', N);
end

P  = turbo_params.P;
Q0 = turbo_params.Q0; Q1 = turbo_params.Q1; Q2 = turbo_params.Q2; Q3 = turbo_params.Q3;

%% ---------------- 交织（(7-2)~(7-4)）----------------
perm = zeros(1, N);                                  % perm(j+1) = Π(j)
for j = 0:N-1
    switch mod(j, 4)
        case 0, qj = 0;
        case 1, qj = 4*Q1;
        case 2, qj = 4*Q0*P + 4*Q2;
        case 3, qj = 4*Q0*P + 4*Q3;
    end
    perm(j+1) = mod(P*j + qj + 3, N);
end
if numel(unique(perm)) ~= N
    error('rcs_turbo_encode:perm', '交织  不是置换（N=%d, P=%d, Q=%d/%d/%d/%d）', N, P, Q0, Q1, Q2, Q3);
end

data_int = zeros(K, 1);
for j = 0:N-1
    i = perm(j+1);
    a = data(2*i+1); b = data(2*i+2);
    if mod(j, 2) == 0                            % j 偶数：交换 A/B
        tmp = a; a = b; b = tmp;
    end
    data_int(2*j+1) = a;
    data_int(2*j+2) = b;
end

%% ---------------- 循环初态（表 7-13）----------------
T13 = table_7_13();
r = mod(N, 15);
if r == 0
    error('rcs_turbo_encode:T13', 'N mod 15 = 0（N=%d）在表 7-13 中没有对应行', N);
end
[~, ~, s1] = crsc_core(data,     0);
[~, ~, s2] = crsc_core(data_int, 0);
C1 = T13(r, s1 + 1);
C2 = T13(r, s2 + 1);

[Y1, W1, ~] = crsc_core(data,     C1);
[Y2, W2, ~] = crsc_core(data_int, C2);

y_period  = turbo_params.y_period;
w_period  = turbo_params.w_period;
y_pattern = char(turbo_params.y_pattern);
w_pattern = char(turbo_params.w_pattern);

Z1 = [puncture_stream(Y1, y_period, y_pattern); puncture_stream(W1, w_period, w_pattern)];
Z2 = [puncture_stream(Y2, y_period, y_pattern); puncture_stream(W2, w_period, w_pattern)];
if numel(Z1) ~= numel(Z2)
    error('rcs_turbo_encode:Z', 'Z1 (%d) 与 Z2 (%d) 长度不等', numel(Z1), numel(Z2));
end

% 输出顺序（§7.3.7.1.4：all couples of systematic bits (A,B) are transmitted first,
% followed by ... (Z1,Z2)；符号索引表 N->Z1,0, N+1->Z2,0, N+2->Z1,1 ...）：
% [系统比特(A,B 自然序); (Z1,0,Z2,0), (Z1,1,Z2,1), ...]
parity = zeros(2*numel(Z1), 1);
parity(1:2:end) = Z1;
parity(2:2:end) = Z2;

turbo_out = [data; parity];
end


% ======================================================================
function [Y, W, Slast] = crsc_core(data, init_state)
%CRSC_CORE 每个 couple 移位一次的双二进制递归系统卷积编码器核心
K = numel(data);
N = K/2;
S = zeros(1, 4);
for k = 4:-1:1
    S(k) = mod(floor(init_state / 2^(4-k)), 2);     % S1..S4
end
Y = zeros(N, 1);
W = zeros(N, 1);
for j = 1:N
    A = data(2*j-1);
    B = data(2*j);
    u = mod(A + B + S(3) + S(4), 2);                % tap 1：A、B 与反馈抽头 S3,S4
    Y(j) = mod(u + S(1) + S(2) + S(4), 2);          % 1 + x1 + x2 + x4
    W(j) = mod(u + S(2) + S(3) + S(4), 2);          % 1 + x2 + x3 + x4
    S = [u, S(1:3)];                                % 移位：一次/couple
end
Slast = S(1)*8 + S(2)*4 + S(3)*2 + S(4);
end


% ======================================================================
function out = puncture_stream(stream, period, pattern)
%PUNCTURE_STREAM 按周期 pattern 作用在流上（'1' 保留，'0' 删除）
stream = stream(:);
pattern = strtrim(char(pattern));
if isempty(pattern) || period <= 0
    out = stream; return;
end
keep = false(numel(stream), 1);
for k = 0:numel(stream)-1
    ch = pattern(mod(k, period) + 1);               % 周期性重复图样
    keep(k+1) = (ch == '1');
end
out = stream(keep);
end


% ======================================================================
function T = table_7_13()
%TABLE_7_13 表 7-13：初始循环状态（行 = N mod 15 = 1..14，列 = 末状态 0..15）
T = [ 0 14  3 13  7  9  4 10 15  1 12  2  8  6 11  5; ...
      0 11 13  6 10  1  7 12  5 14  8  3 15  4  2  9; ...
      0  8  9  1  2 10 11  3  4 12 13  5  6 14 15  7; ...
      0  3  4  7  8 11 12 15  1  2  5  6  9 10 13 14; ...
      0 12  5  9 11  7 14  2  6 10  3 15 13  1  8  4; ...
      0  4 12  8  9 13  5  1  2  6 14 10 11 15  7  3; ...
      0  6 10 12  5  3 15  9 11 13  1  7 14  8  4  2; ...
      0  7  8 15  1  6  9 14  3  4 11 12  2  5 10 13; ...
      0  5 14 11 13  8  3  6 10 15  4  1  7  2  9 12; ...
      0 13  7 10 15  2  8  5 14  3  9  4  1 12  6 11; ...
      0  2  6  4 12 14 10  8  9 11 15 13  5  7  3  1; ...
      0  9 11  2  6 15 13  4 12  5  7 14 10  3  1  8; ...
      0 10 15  5 14  4  1 11 13  7  2  8  3  9 12  6; ...
      0 15  1 14  3 12  2 13  7  8  6  9  4 11  5 10];
end
