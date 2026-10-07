function symbols = rcs_modulation(bits, params, rotate)

bits = bits(:);
m = params.mod_order;
if nargin < 3 || isempty(rotate)
    rotate = true;
end
n_symbols = ceil(length(bits) / m);
pad_len = n_symbols * m - length(bits);
if pad_len > 0
    bits = [bits; zeros(pad_len, 1)];
end
bit_groups = reshape(bits, m, n_symbols).';

% systematic bit count K = 2N (from the turbo parameters when available)
K = [];
if isfield(params, 'turbo_params') && isstruct(params.turbo_params) && ...
   isfield(params.turbo_params, 'N') && ~isempty(params.turbo_params.N)
    K = 2 * params.turbo_params.N;
end
if ~isempty(K) && K > numel(bits)
    K = [];
end

switch m
    case 1  % pi/2-BPSK (table 7-16)
        u = 1 - 2*bit_groups(:,1);   % +1 / -1, |u| = 1 (constellation radius 1)
        if rotate
            n_idx = (0:n_symbols-1).';
            symbols = u .* exp(1j*(n_idx*pi/2 + pi/4));
        else
            symbols = u;
        end

    case 2  % QPSK (table 7-17)
        I = 1 - 2*bit_groups(:,1);
        Q = 1 - 2*bit_groups(:,2);
        symbols = (I + 1j*Q)/sqrt(2);

    case 3  % 8PSK (tables 7-18 / 7-19)
        ph = (1:2:15).' * pi/8;
        lut = exp(1j*ph);
        gray_map = [0 0 0; 0 0 1; 1 0 1; 1 0 0; 1 1 0; 1 1 1; 0 1 0; 0 1 1];
        if isempty(K)
            symbols = zeros(n_symbols, 1);
            for i = 1:n_symbols
                idx = find(all(gray_map == bit_groups(i, :), 2));
                symbols(i) = lut(idx);
            end
        else
            sys = bits(1:K);
            par = bits(K+1:end);
            si = 1; pi_i = 1;
            symbols = zeros(n_symbols, 1);
            for i = 1:n_symbols
                if pi_i <= numel(par)
                    u0 = par(pi_i); pi_i = pi_i + 1;
                    u1 = bget(sys, si); si = si + 1;
                    u2 = bget(sys, si); si = si + 1;
                else
                    u0 = bget(sys, si); si = si + 1;
                    u1 = bget(sys, si); si = si + 1;
                    u2 = bget(sys, si); si = si + 1;
                end
                idx = find(all(gray_map == [u0 u1 u2], 2));
                symbols(i) = lut(idx);
            end
        end

    case 4  % 16QAM (tables 7-20 / 7-21)
        pam = [-1; 1; -3; 3]/sqrt(10);
        if isempty(K)
            symbols = zeros(n_symbols, 1);
            for i = 1:n_symbols
                b = bit_groups(i, :);
                symbols(i) = pam(b(3)*2 + b(4) + 1) + 1j*pam(b(1)*2 + b(2) + 1);
            end
        else
            sys = bits(1:K);
            par = bits(K+1:end);
            n_sym = n_symbols;
            % omit a trailing symbol that would carry only a solitary parity bit (7.3.7.1.4.4)
            if n_sym >= 1 && (n_sym - 1) * 3 >= K && mod(K + numel(par), 4) == 1
                n_sym = n_sym - 1;
            end
            symbols = zeros(n_sym, 1);
            si = 1; pj = 1;
            for i = 1:n_sym
                if pj <= numel(par)
                    uQ1 = par(pj);      pj = pj + 1;
                    uQ0 = bget(sys, si); si = si + 1;
                    uI1 = bget(sys, si); si = si + 1;
                    uI0 = bget(sys, si); si = si + 1;
                else
                    uQ1 = bget(sys, si); si = si + 1;
                    uQ0 = bget(sys, si); si = si + 1;
                    uI1 = bget(sys, si); si = si + 1;
                    uI0 = bget(sys, si); si = si + 1;
                end
                Qv = pam(uQ1*2 + uQ0 + 1);
                Iv = pam(uI1*2 + uI0 + 1);
                symbols(i) = Iv + 1j*Qv;
            end
        end

    otherwise
        error('modulation type error');
end
end


% ======================================================================
function v = bget(vec, k)
%BGET vec(k) with 0 padding when k runs past the end (7.3.7.1.4 null filling)
if k <= numel(vec)
    v = vec(k);
else
    v = 0;
end
end
