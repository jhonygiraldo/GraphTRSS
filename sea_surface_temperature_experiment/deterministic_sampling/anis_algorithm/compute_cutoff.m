function [ omega ] = compute_cutoff( L_k, k, S )
%COMPUTE_CUTOFF
%   AUTHOR: Aamir Anis, USC
%   This function computes the cutoff frequency for a given
%   sampling set

% % %
% PARAMETER DESCRIPTION
% 
% INPUT
% L: kth power of Laplacian
% S: Sampling set, input as a logical vector
% k: Power of Laplacian while computing cutoff, higher the order,
% greater the accuracy, but the complexity is also higher.
% 
% OUTPUT
% omega: cutoff frequency of the given set
% 
% % %

% Symmetrization
L_k = 0.5*(L_k+L_k.');

% compute minimum eigen-pair: efficient way
[~,omega] = eigs(L_k(~S,~S),1,'sm');
omega = abs(omega)^(1/k);

