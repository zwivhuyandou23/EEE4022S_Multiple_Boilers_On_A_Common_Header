function [Qsp1, Qsp2] = allocator(Q_star1, Q_star2, Q_bias)
% Simple allocator: splits bias dQ between two boilers
w1 = 0.9;
w2 = 0.1;
Qsp1 = Q_star1 + w1 * Q_bias;
Qsp2 = Q_star2 + w2 * Q_bias;
end