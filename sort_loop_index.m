function[loop_list_sorted]=sort_loop_index(loop_list)
%SORT_LOOP_INDEX Sort the loop indices such that the smallest index comes
%first.
%
% LOOP_LIST_SORTED=SORT_LOOP_INDEX(LOOP_LIST)
% The input LOOP_LIST has the form of the output of functions from the FBLDetect
% function suite, or at least one column named "loop" that contains a cell
% array with each entry capturing the order of the node indices (integer values) forming the loop.
%
%Example
%   J = [-1 0 0 -1; 1 -1 0 1; 0 1 -1 0; 0 0 1 -1];
%   loop_t = find_loops(J);
%   loop_t.loop{1}
%   loop_t_sorted=sort_loop_index(loop_t);
%   loop_t_sorted.loop{1}
% 
% See also FIND_LOOPS, FIND_LOOPS_NOSCC, FIND_LOOPS_VSET

function[outvec]=sort_vector(inputvec)
    [a,b]=min(inputvec);
    outvec=[inputvec(b:end),inputvec(2:b)];
end

loop_list_sorted_loops=cellfun(@sort_vector,loop_list.loop,'UniformOutput',false);

loop_list_sorted=loop_list;
%replace the loop lists by sorted loop lists
loop_list_sorted.loop=loop_list_sorted_loops;
end

