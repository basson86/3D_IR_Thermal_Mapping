function [Temp_map]=TempConvert(Im_ir)

m0 = -53.771;
m1 = 0.0045575;
m2 = -1.1621e-07;
m3 = 1.692e-012;
m4 = -9.9176e-18;

Temp_map=m0 + m1.*Im_ir  + m2 .* Im_ir.^2 + m3 .* Im_ir.^3 + m4 .*Im_ir.^4;

