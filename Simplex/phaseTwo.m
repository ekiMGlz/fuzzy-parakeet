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
    % If it is, obasis and obfs will be the witness
    if ub
        bound = 0;
        witness = null_vars(candidates(end, 1)-1);
        obasis = [basic_vars, witness];
        sbfs(witness) = inf;
        obfs = sbfs;
        oval = inf;
        return
    end

end

function [T] = nextTableau(T, enter, leave)
    % Calculates the next tableau
    % Input:
    % T - (m)x(n) matrix with the previous tableau, arranged as stated in
    %     phaseTwo
    % enter - index of the entering variable, must be a value in [2, n]
    % leave - index of the leaving variable, value in [1, m-1]
    % Output:
    % T - The next tableau
    % 
    
end

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
    % index of the entering variable, the index of the leaving variable and
    % the new value of the entering variable
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
            
            candidates = [candidates; j, min_index, min];
            
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