function [data_with_crc, crc_value] = rcs_crc16_encode(data)


data = data(:);

reg = zeros(1, 16);

for i = 1:numel(data)
    feedback = mod(data(i) + reg(16), 2);

    new_reg = zeros(1, 16);
    new_reg(1) = feedback;                  
    new_reg(2:16) = reg(1:15);              
    new_reg(16) = mod(new_reg(16) + feedback, 2);   
    new_reg(3)  = mod(new_reg(3)  + feedback, 2);   

    reg = new_reg;
end

crc_value = flipud(reg(:));
data_with_crc = [data; crc_value];

end
