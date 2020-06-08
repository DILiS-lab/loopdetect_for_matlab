function[loop_ind]=find_edges(loop_list, source_node, target_node)
% FIND_EDGES Find loops in a loop list that contain a certain direct 
% regulation (edge).
%
% LOOP_IND = FIND_EGDES(LOOP_LIST, SOURCE_NODE, TARGET_NODE) 
% The function returns a list of indices LOOP_IND within the provided
% list of feedback loops LOOP_LIST of the loops containing the edge from 
% SOURCE_NODE to TARGET_NODE. 
%
% Example 
% This call would extract the indices of all loops in which variable 2 is
% regulated by variable 1:
% loop_edge_ind=find_edges(loop_list,1,2)
% loop_list(loop_edge_ind,:) %returns the loops containing the regulation
% from variable 2 to variable 1
% 
% See also: find_loops(), find_loops_noscc(), find_loops_vset()


%search source_node and its index in all loops, determine the following
%node index and check whether this is the target_node
loops_follow_source_node_is_target_node=cellfun(@(z) isequal(z(find(z==source_node,1)+1),target_node),loop_list.loop);

%this outputs a logical array of length loop_list whose entries are 
% - 1 if the loop contains source_node and the consecutive node is
% target_node
% - 0 if the loop does not contain source_node or if the loop contains 
% source_node and the consecutive node is not target_node 

%transfer the logical array to indices of the loops (i.e. indices where
%the previous array is 1
loop_ind=find(loops_follow_source_node_is_target_node);
    
end