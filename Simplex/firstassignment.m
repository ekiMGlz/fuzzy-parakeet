% First Assignment for the Linear Programming Class
% Authors:
% Alejandro Rodriguez Orozco - 
% Miguel Gonzalez Borja - 155766

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                                 Runtime                                  %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global verbose setting, set to 1 to print additional info during each
% funciton, set to 0 to ommit
global verbose
verbose = 1;

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
% A = [1, -1, 1, 0;
%       -1, 1, 0, 1];
% 
% b = [1;
%      2];
%
% c = [1;
%      0;
%      0;     
%      0];
%
% sbasis = [3, 4];
% sbfs = [0, 0, 1, 2];

% Degenerate Example (Matousek, pg 62)
%A = [-1, 1, 1, 0;
%      1, 0, 0, 1];

%b = [0;
%     2];

%c = [0;
%     1;
%     0;     
%     0];

%sbasis = [3, 4];
%sbfs = [0, 0, 0, 2];

% No bfs example (Matousek, pg 64)
A = [1, 3, 1;
     0, 2, 1];

b = [4;
     2];

c = [1;
     2;
     0];

% Ejemplo Andrea
%A = [-1, -1];
%b = [0];
%c = [1; 
%     0];

% [bound, obasis, obfs oval] = phaseTwo(A, b, c, sbasis, sbfs);
[status, obasis, obfs, oval] = bothPhases(A, b, c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                              Main Functions                              %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[nvac, basis, bfs] = phaseOne(A, b, c)
    % maximise c^T x
    % subject to Ax = b, x >= 0, b >=0
    %
    % Input:
    % A mxn matrix with m <= n and rank of A is m
    % b column vector with m rows
    % c column vector with n rows
    %
    % Output:
    % nvac = 0 if the feasible set is empty
    % nvac = 1 if the feasible set is non-empty
    % basis = a vector of size m of indices of column vectors for a feasible basis for
    % the problem if the feasible set is non-empty
    % bfs = a vector of size n of the basic feasible solution corresponding to this
    % basis (if the feasible set is non-empty)
    
    global verbose
    
    % Save the size of A for future use
    [m, n] = size(A);
    
    % Create a new auxiliary linear program
    % max sum -z_i
    % subject to Dx=b
    % With D = (A | I_m) and introducing new m new variables zi
    D = [A, eye(m)];
    sbasis = (n+1):(n+m);
    sbfs = [zeros(1, n), transpose(b)];
    z = [zeros(n, 1); -ones(m, 1)];
    
    if verbose
        fprintf("Creating auxiliary linear program...\n")
    end
    
    % Solve the auxiliary linear program
    [bound, obasis, obfs, oval] = phaseTwo(D, b, z, sbasis, sbfs);
    
    if verbose
        fprintf("Auxiliary linear program solved, checking feasibility...\n")
    end
    
    % Check if the original probleam was feasible
    if oval == 0
        basis = setdiff(obasis, (m+1):(m+n));
        
        k = length(basis);
        % Check the new basis to make sure it has m elements
        if k<m
            if verbose
                fprintf("Auxiliary basis contains correction variables, adjusting initial basis...\n");
            end
            
            % Get non basic variables and add enough of them to our new
            % basis to get m variables
            null_vars = setdiff(1:m, basis);
            basis = sort([basis, null_vars(1:(m-k))]);
        end
        
        bfs = obfs(1:n);
        nvac = 1;
    else
        basis = [];
        bfs = [];
        nvac = 0;
    end
    
end

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
    %          optimal feasible basis for the problem if the problem is bounded
    % obfs = a vector of size n which is the optimal basic feasible solution corresponding to this optimal basis if the problem is bounded
    % oval = the objective value of this optimal basic feasible solution (if the problem is bounded)

    global verbose
    
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
    
    
    % Check the tableau
    [bound, optimal, degenerate, candidates] = tableauStatus(T);
    
    if verbose
        fprintf("Initial Tableau:\n")
        printTableau(T, basic_vars, null_vars);
    end
    
    % If the tableau shows a possible unbounded soln, exit
    if ~bound
        if verbose
            fprintf("Unbounded solution detected...\n")
        end
        obasis = [];
        obfs = [];
        oval = inf;
        return
    end
    
    % Loop until an optimal solution is reached
    while ~optimal
        
        % TODO: add different steps if non degenerate bfs
        % Take a step using Bland's Rule
        [T, basic_vars, null_vars] = BlandStep(T, basic_vars, null_vars, candidates);
        
        if verbose
            fprintf("New Tableau:\n")
            printTableau(T, basic_vars, null_vars);
        end
        
        % Check the tableau
        [bound, optimal, degenerate, candidates] = tableauStatus(T);
        
        % If the tableau shows a possible unbounded soln, exit
        if ~bound
            if verbose
                fprintf("Unbounded solution detected...\n")
            end
            obasis = [];
            obfs = [];
            oval = inf;
            return
        end
    end
    
    if verbose
        fprintf("Optimal solution found...\n")
    end
    
    % Get the correct values for the optimal basis, bfs and value
    obasis = sort(basic_vars);
    obfs = zeros(1, n);
    for i = 1:m
        obfs(basic_vars(i)) = T(i, 1);
    end
    oval = T(m+1, 1);
    
end

function[status, obasis, obfs, oval] = bothPhases(A, b, c)
    % maximise c^T x
    % subject to Ax = b, x >= 0, b >=0
    %
    % Input:
    % A mxn matrix with m <= n and rank of A is m
    % b column vector with m rows
    % c column vector with n rows
    %
    % Output:
    % status = -1 if the feasible set is empty
    % status = 0 if the feasible set is non-empty but the problem is unbounded 
    %          (there is no optimal solution)
    % status = 1 if the problem is bounded (there is an optimal solution)
    % obasis = a vector of size m of indices of an optimal feasible basis for the 
    %          problem if the feasible set is non-empty and the problem is bounded (in terms of a 
    %          set of indices of column vectors)
    % obfs = a vector of size n which is the optimal basic feasible solution 
    %        corresponding to this optimal basis if the feasible set is non-empty and the problem is bounded
    % oval = the objective value of this optimal basic feasible solution (if the 
    %        feasible set is non-empty and the problem is bounded)
    
    global verbose
    
    if verbose
        fprintf("Verbose Mode On\n\n")
        fprintf("Obtaining starting basic feasible soln...\n")
    end
    
    % Obtain a bfs for the problem
    [nvac, sbasis, sbfs] = phaseOne(A, b, c);
    
    % Check feasibility
    if ~nvac
        if verbose
            fprintf("Feasible set is empty, exiting with status -1\n")
        end
        status = -1;
        obasis = [];
        obfs = [];
        oval = inf;
        return
    end
    
    if verbose
        fprintf("Starting bfs found:\n")
        sbfs
        sbasis
        fprintf("Moving on to Phase Two...\n")
    end
    
    % If feasible, solve the problem
    [bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs);
    
    % Check bound
    if ~bound
        if verbose
            fprintf("Unbounded problem, exiting with status 0\n")
        end
        status = 0;
        obasis = [];
        obfs = [];
        oval = inf;
        return
    end
    
    if verbose
        fprintf("Solution found, exiting with status 1\n")
    end
    % Problem solved
    status = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                          %%%
%%%                           Auxiliary Functions                            %%%
%%%                                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bound, optimal, degenerate, candidates] = tableauStatus(T)
    % Status of the tableau
    % Input:
    %   T - tableau to check
    % Output:
    %   bounded - 1 if bounded, else 0
    %   optimal - 1 if optimal, else 0
    %   degenerate - 1 if degenerate bfs, else 0
    %   candidates - (k)x(3) matrix, where k is the number of possible
    %                enetring varialbles for the simplex. Each row contains, in order, the
    %                index of the column representing the entering variable, 
    %                the index of the row representing the leaving variable,
    %                and the new value for the entering variable.
    [m, n] = size(T);
    bound = 1;
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
            
            candidates = [candidates; j, min_index, min];
            
            if min_index == -1
                bound = 0;
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
    %   T - (m)x(n) matrix with the previous tableau, arranged as stated in
    %       phaseTwo
    %   col - index of the row representing the entering variable, 
    %         must be a value in [2, n]
    %   row - index of the row representing the leaving variable, 
    %         value in [1, m-1]
    % Output:
    %   T - The next tableau
    % 
    
    [m, n] = size(T);
    % Save the previous column
    aux_col = T(:, col);
    
    % Replace the entering column with [0, 0, ..., -1, 0, ..., 0] (-1 in leaving row)
    T(:, col) = zeros(m, 1);
    T(row, col) = -1;
    
    % Divide the leaving row by the coeficient of the entering variable.
    T(row, :) = -T(row, :)/aux_col(row);
    % Adjust the rest of the rows to the new variable
    for i = 1:m
        if i ~= row
            T(i, :) = T(i, :) + aux_col(i)*T(row, :);
        end
    end
    
end

function[T, basic_vars, null_vars] = BlandStep(T, basic_vars, null_vars, candidates)
    % Takes a step using Bland's Rule
    % Input:
    %   T - (m)x(n) matrix with the previous tableau, arranged as stated in
    %       phaseTwo
    %   basic_vars - index of the row representing the entering variable, 
    %                must be a value in [2, n]
    %   null_vars - index of the row representing the leaving variable, 
    %               value in [1, m-1]
    %   candidates - (k)x(3) matrix containing possible candidates for 
    %                entering the basis
    % Output:
    %   T - The next tableau
    %   basic_vars - vector with the index for the basic variables
    %   null_vars - vector with the index for the non basic variables

    global verbose
    
    % Decide on the entering var
    % Since the columns of the tableau are not necesarilly sorted with respect 
    % to the index of each variable, the minimum index must be manually found
    enter = candidates(1, 1);
    leave = candidates(1, 2);
    min = null_vars(enter-1);
    for r = transpose(candidates(2:end, :))
        if null_vars(r(1)-1) < min
            min = null_vars(r(1)-1);
            enter = r(1);
            leave = r(2);
        end
    end
    
    % Generate the next tableau
    T = nextTableau(T, enter, leave);
    
    if verbose
        fprintf("Pivot using Bland's Rule:\n")
        fprintf("\tEntering Variable: x_%d\n", null_vars(enter-1))
        fprintf("\tLeaving Variable: x_%d\n", basic_vars(leave))
    end
    
    % Swap the entering and leaving variables in basic_vars and null_vars
    aux = null_vars(enter-1);
    null_vars(enter-1) = basic_vars(leave);
    basic_vars(leave) = aux;
end

function[] = printTableau(T, basic_vars, null_vars)
    
    S = "      \tp\t\t" + sprintf("x_%d\t\t", null_vars) + "\n";
    [m, n] = size(T);
    for i = 1:m-1
        S = S + sprintf("x_%d =\t", basic_vars(i)) + sprintf("%0.2f\t", T(i, :)) + "\n";
    end
    
    S = S + "z   =\t" + sprintf("%0.2f\t", T(end, :)) + "\n";
    fprintf(S)
end