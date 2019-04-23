% First Assignment for the Linear Programming File
% Authors:
% Alejandro Rodriguez Orozco - 
% Miguel Gonzalez Borja - 155766

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                                 Runtime                                  %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First example (Matousek, pg 57)
% A = [-1, 1, 1, 0, 0;
%       1, 0, 0, 1, 0;
%       0, 1, 0, 0, 1];
% 
% b = [1;
%      3;
%      2];
% 
% c = [1;
%      1;
%      0;     
%      0;
%      0];
% 
% sbasis = [3, 4, 5];
% sbfs = [0, 0, 1, 3, 2];

% Unbounded Example (Matousek, pg 61)
A = [1, -1, 1, 0;
      -1, 1, 0, 1];

b = [1;
     2];

c = [1;
     0;
     0;     
     0];

sbasis = [3, 4];
sbfs = [0, 0, 1, 2];


[bound, obasis, obfs oval] = phaseTwo(A, b, c, sbasis, sbfs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                              Main Functions                              %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs)
    % maximise c^T x
    % subject to Ax = b, x >= 0, b >=0
    %
    % Input:
    % A mxn matrix with m <= n and rank of A is m
    % b column vector with m rows
    % c column vector with n rows
    % sbasis a vector of size m of indices of column vectors for a feasible basis for this problem from which to start the simplex method
    % sbfs a vector of size n which is the basic feasible solution corresponding to this basis
    %
    % Output:
    % bound = 0 if the problem is unbounded (there is no optimal solution)
    % bound = 1 if the problem is bounded (there is an optimal solution)
    % obasis = a vector of size m of indices of column vectors which gives an
    % optimal feasible basis for the problem if the problem is bounded
    % obfs = a vector of size n which is the optimal basic feasible solution corresponding to this optimal basis if the problem is bounded
    % oval = the objective value of this optimal basic feasible solution (if the problem is bounded)

    % Save the size of A for future use
    [m, n] = size(A);

    % Define basic_vars and null_vars as vectors with the indices for the
    % basic and non basic vars.
    basic_vars = sbasis;
    null_vars = setdiff(1:n, sbasis);

    % Calculate the tableau T as a (m+1)x(n-m+1) matrix, arranged as follows:
    %   |p   Q|
    %   |z0  r|
    % Where:
    %   p  = inv(A_B)*b  (or sbfs(basic_vars).T)
    %   Q  = -inv(A_B)*A_N
    %   z0 = c_B.T*p
    %   r  = c_N + (Q.T*c_B)
    T = zeros(m+1, n-m+1);
    Q = -A(:, basic_vars)\A(:, null_vars);

    T(1:m, 1) = transpose(sbfs(basic_vars));
    T(1:m, 2:(n-m+1)) = Q;
    T(m+1, 1) = dot(c(basic_vars), sbfs(basic_vars));
    T(m+1, 2:(n-m+1)) = transpose(c(null_vars) + transpose(Q)*c(basic_vars));
    
    T

    [ub, op, deg, candidates] = tableauStatus(T);
    
    % Check if the problem is unbounded
    if ub
        bound = 0;
        obasis = [];
        obfs = [];
        oval = -inf;
        return
    end
    
    while ~op
        
        % TODO: add steepest edge step if non degenerate bfs
        [T, enter, leave] = BlandsRuleStep(T, candidates);
        aux = null_vars(enter-1);
        
        null_vars(enter-1) = basic_vars(leave);
        basic_vars(leave) = aux;
        % TODO: make a pretty print for tableaus
        T
        basic_vars
        null_vars
        [ub, op, deg, candidates] = tableauStatus(T);
        
        if ub
            bound = 0;
            obasis = [];
            obfs = [];
            oval = -inf;
            return
        end
    end
    
    bound = 1;
    obasis = sort(basic_vars);
    obfs = zeros(1, n);
    for i = 1:m
        obfs(basic_vars(i)) = T(i, 1);
    end
    oval = T(m+1, 1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                           Auxiliary Functions                            %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [unbounded, optimal, degenerate, candidates] = tableauStatus(T)
    % Status of the tableau
    % Input:
    % T - tableau to check
    % Output:
    % unbounded - 1 if unbounded, else 0
    % optimal - 1 if optimal, else 0
    % degenerate - 1 if degenerate bfs, else 0
    % candidates - (k)x(3) matrix, where k is the number of possible
    % enetring varialbles for the simplex. Each row contains, in order, the
    % index of the entering variable and the index of the leaving variable 
    [m, n] = size(T);
    unbounded = 0;
    optimal = 1;
    degenerate = 0;
    candidates = [];
    
    % unbounded, optimal and candidates check
    for j = 2:n
        if T(m, j) > 0
            % Entering candidate found, find corresponding leaving var
            optimal = 0;
            min = inf;
            min_index = -1;
            for i = 1:m-1
                if T(i, j) < 0 && -T(i, 1)/T(i, j) < min
                    min = -T(i, 1)/T(i, j);
                    min_index = i;
                end
            end
            
            candidates = [candidates; j, min_index];
            
            if min_index == -1
                unbounded = 1;
                break;
            end
        end
    end
    
    % degenerate check
    for i = 1:m-1
        if T(i, 1) == 0
            degenerate = 1;
            break;
        end
    end
end

function [T] = nextTableau(T, col, row)
    % Calculates the next tableau
    % Input:
    % T - (m)x(n) matrix with the previous tableau, arranged as stated in
    %     phaseTwo
    % col - index of the row representing the entering variable, must be a value in [2, n]
    % row - index of the row representing the leaving variable, value in [1, m-1]
    % Output:
    % T - The next tableau
    % 
    
    [m, n] = size(T);
    % Save the previous column
    aux_col = T(:, col);
    
    %Replace the entering column with [0, 0, ..., -1, 0, ..., 0] (-1 in leaving row)
    T(:, col) = zeros(m, 1);
    T(row, col) = -1;
    
    %Divide the leaving row by the coeficient of the entering variable.
    T(row, :) = -T(row, :)/aux_col(row);
    for i = 1:m
        if i ~= row
            T(i, :) = T(i, :) + aux_col(i)*T(row, :);
        end
    end
    
end

function[T, enter, leave] = BlandsRuleStep(T, candidates)
    % TODO: check if valid, otherwise check for min index (might need more args)
    enter = candidates(1, 1);
    leave = candidates(1, 2);
    T = nextTableau(T, enter, leave);
end


