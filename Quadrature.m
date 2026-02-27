    classdef Quadrature < handle
        properties
            N % Polynomial degree
            Points % Gauss-Lobatto-Legendre points
            Weights % Quadrature weights
            Polynomials % Associated Legendre polynomials
        end
        
        methods
            % Constructor
            function obj = Quadrature(N)
                if nargin < 1 || N < 1
                    error('Quadrature: Invalid input. N must be >= 1.');
                end
                obj.N = N;
                [obj.Points, obj.Weights, obj.Polynomials] = obj.computeNodesAndWeights(N);
            end
            
            % Method to compute points and weights
            function [x, w, P] = computeNodesAndWeights(~, N)
                % N1 = N + 1 for convenience
                N1 = N + 1;
    
                % Initial guess for points: Chebyshev-Gauss-Lobatto
                x = cos(pi * (0:N) / N)';
    
                % Storage for Legendre polynomials
                P = zeros(N1, N1);
    
                % Iterative Newton-Raphson refinement
                xold = 2 * ones(size(x)); % Initial dummy value
                while max(abs(x - xold)) > eps
                    xold = x;
    
                    % Evaluate Legendre polynomials and derivatives
                    P(:, 1) = 1; % P_0
                    P(:, 2) = x; % P_1
    
                    for k = 2:N
                        P(:, k + 1) = ((2 * k - 1) * x .* P(:, k) - (k - 1) * P(:, k - 1)) / k;
                    end
    
                    % Newton-Raphson step
                    x = xold - (x .* P(:, N1) - P(:, N)) ./ (N1 * P(:, N1));
                end
    
                % Compute weights
                w = 2 ./ (N * N1 * P(:, N1).^2);
                x = -x;
            end
            
            % Method to display points and weights
            function displayPointsAndWeights(obj)
                disp('Gauss-Lobatto-Legendre Points:');
                disp(obj.Points);
                disp('Quadrature Weights:');
                disp(obj.Weights);
            end
            
            % Method to evaluate a function at quadrature points
            function values = evaluateFunction(obj, funcHandle)
                if nargin < 2 || ~isa(funcHandle, 'function_handle')
                    error('Quadrature: Invalid function handle input.');
                end
                values = funcHandle(obj.Points);
            end
        end
    end
