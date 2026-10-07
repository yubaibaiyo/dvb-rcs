function [term_values, n_term] = rcs_cpm_termination(Vn, ph)

switch ph
    case 7
        table = [0 0 0; 3 3 0; 3 2 0; 3 1 0; 3 0 0; 2 0 0; 1 0 0];   % 表 7-25
    case 5
        table = [0 0 0; 3 1 0; 3 0 0; 2 0 0; 1 0 0];                  % 表 7-26
    case 4
        table = [0 0; 3 0; 2 0; 1 0];                                 % 表 7-27
    case 3
        table = [0 0; 2 0; 1 0];                                      % 表 7-28
    otherwise
        error('rcs_cpm_termination:ph', ...
              'ph = %d 无相位终止表（表 7-25~7-28 仅覆盖 7/5/4/3）', ph);
end
Vn = mod(Vn, ph);
term_values = table(Vn + 1, :);
n_term = numel(term_values);
end
