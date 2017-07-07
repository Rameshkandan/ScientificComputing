function [y] = Fun_ExpEuler(y_0, A, tau)
    y = y_0 + tau * A * y_0;

