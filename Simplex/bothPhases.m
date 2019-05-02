% Full Simplex function
% Authors:
% Alejandro Rodriguez Orozco - 
% Miguel Gonzalez Borja - 155766

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
    
    % Set debug to 1 to print additional info during the execution of the
    % function. Set to 0 do suppress 
    debug = 1;
    
    if debug
        fprintf("Debug Mode On\n\n")
        fprintf("Obtaining starting basic feasible soln...\n")
    end
    
    % Obtain a bfs for the problem
    [nvac, sbasis, sbfs] = phaseOne(A, b, c);
    
    % Check feasibility
    if ~nvac
        if debug
            fprintf("Feasible set is empty, exiting with status -1\n")
        end
        status = -1;
        obasis = [];
        obfs = [];
        oval = inf;
        return
    end
    
    if debug
        fprintf("Starting bfs found:\n")
        sbfs
        sbasis
        fprintf("Moving on to Phase Two...\n")
    end
    
    % If feasible, solve the problem
    [bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs);
    
    % Check bound
    if ~bound
        if debug
            fprintf("Unbounded problem, exiting with status 0\n")
        end
        status = 0;
        obasis = [];
        obfs = [];
        oval = inf;
        return
    end
    
    if debug
        fprintf("Solution found, exiting with status 1\n")
    end
    % Problem solved
    status = 1;
end