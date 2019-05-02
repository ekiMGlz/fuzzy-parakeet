% Examples

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
% A = [1, 3, 1;
%      0, 2, 1];
% 
% b = [4;
%      2];
% 
% c = [1;
%      2;
%      0];

% Ejemplo Andrea
A = [-1, -1];
b = [0];
c = [1; 
     0];

% [bound, obasis, obfs oval] = phaseTwo(A, b, c, sbasis, sbfs);
[status, obasis, obfs, oval] = bothPhases(A, b, c);