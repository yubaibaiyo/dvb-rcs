function [data_with_crc, crc_value] = rcs_crc32_encode(data)

data = data(:);

taps = [26 23 22 16 12 11 10 8 7 5 4 2 1];     

reg = ones(1, 32);                              

for i = 1:numel(data)
    feedback = mod(data(i) + reg(32), 2);

    new_reg = zeros(1, 32);
    new_reg(1) = feedback;                      
    new_reg(2:32) = reg(1:31);                  
    new_reg(taps + 1) = mod(new_reg(taps + 1) + feedback, 2);

    reg = new_reg;
end

crc_value = flipud(reg(:));
data_with_crc = [data; crc_value];

end
