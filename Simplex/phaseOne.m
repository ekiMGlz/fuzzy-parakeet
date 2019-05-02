% Simplex Phase One function
% Authors:
% Alejandro Rodriguez Orozco - 
% Miguel Gonzalez Borja - 155766

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
    
    % Set debug to 1 to print additional info during the execution of the
    % function. Set to 0 do suppress 
    debug = 1;
    
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
    
    if debug
        fprintf("Creating auxiliary linear program...\n")
    end
    
    % Solve the auxiliary linear program
    [bound, obasis, obfs, oval] = phaseTwo(D, b, z, sbasis, sbfs);
    
    if debug
        fprintf("Auxiliary linear program solved, checking feasibility...\n")
    end
    
    % Check if the original probleam was feasible
    if oval == 0
        basis = setdiff(obasis, (m+1):(m+n));
        
        k = length(basis);
        % Check the new basis to make sure it has m elements
        if k<m
            if debug
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

