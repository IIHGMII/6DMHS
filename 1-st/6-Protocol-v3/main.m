clc;clear;
K_rice = 10;

parfor (i = 1:20 , 4)
    [P_out(i)] = test_6DMA_Proposed(K_rice);
end





