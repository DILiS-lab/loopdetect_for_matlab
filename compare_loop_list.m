function[ind_a_id,ind_a_switch,ind_a_notin,ind_b_id,ind_b_switch]=compare_loop_list(loop_list_a,loop_list_b)
% COMPARE_LOOP_LIST Compare two lists of loops. 
%
% COMPARE_LOOP_LIST analyzes which loops in loop_list_a are shared between both lists, have 
% switched sign or do not occur in loop_list_b. Loop lists should be outputs of
% functions such as find_loops(), find_loops_noscc(), find_loops_vset(). 
% Internally, the function performs sorting of the loops (relying on sort_loop_index) in
% order to enable comparison.
%
% IND_A_ID = COMPARE_LOOP_LIST(LOOP_LIST_A,LOOP_LIST_B) returns the indices
% of loops in LOOP_LIST_A that are also in LOOP_LIST_B (same loop and same
% sign).
%
% [IND_A_ID,IND_A_SWITCH,IND_A_NOTIN] = COMPARE_LOOP_LIST(LOOP_LIST_A,LOOP_LIST_B)
% returns in addition the indices of loops in LOOP_LIST_A that occur in
% LOOP_LIST_B but with switched sign (positive feedback loop instead of
% negative feedback loop or vice versa), and the indices of the loops that
% do not occur in LOOP_LIST_B
%
% [IND_A_ID,IND_A_SWITCH,IND_A_NOTIN,IND_B_ID,IND_B_SWITCH] = COMPARE_LOOP_LIST(LOOP_LIST_A,LOOP_LIST_B)
% returns also the indices of the loops in LOOP_LIST_B, IND_B_ID, corresponding to the
% loops in LOOP_LIST_B given by IND_A_ID, as well as the indices of the loops in LOOP_LIST_B,
% IND_B_SWITCH, that correspond to the
% loops in LOOP_LIST_B given by IND_A_SWITCH. Thus, LOOP_LIST_A(IND_A_ID,:) is exactly
% identical to LOOP_LIST_B(IND_B_ID,:) (also in the same order), up to
% changes in with which node index is started in the description of the
% loops. LOOP_LIST_A(IND_A_SWITCH,:) differs from
% LOOP_LIST_B(IND_B_SWITCH,:) only in the 'sign' column (and possibly in
% the first reported node index of loops).
%
% See also: FIND_LOOPS, FIND_LOOPS_NOSCC, FIND_LOOPS_VSET, SORT_LOOP_INDEX

% Examples
% length(ind_a_id)==height(loop_list_a) means all loops from list a are
% also in list b
% loop_list_a(ind_a_id(2),:) is the same loop as loop_list_b(ind_b_id(2),:)
% 



function[loop_list_sorted]=sort_loop_index(loop_list)
   %sort the loop indices to start with lowest index in all loops
    function[outvec]=sort_vector(inputvec)
        [a,b]=min(inputvec);
        outvec=[inputvec(b:end),inputvec(2:b)];
    end

    loop_list_sorted_loops=cellfun(@sort_vector,loop_list.loop,'UniformOutput',false);

    loop_list_sorted=loop_list;
    %replace the loop lists by sorted loop lists
    loop_list_sorted.loop=loop_list_sorted_loops;
end

%first: sort the loops in both lists for starting with the lowest index.
loop_list_a_sorted=sort_loop_index(loop_list_a);
loop_list_b_sorted=sort_loop_index(loop_list_b);

%check for each loop in loop_list1 whether there is a loop in loop_list2
%such that the loops are identical (maybe different sign)
is_in_index=cellfun(@(a) find(cellfun(@(z) isequal(a,z), loop_list_b_sorted.loop)),loop_list_a_sorted.loop,'UniformOutput',false);

%indices of these loops in list 1
indices_in_a=find(cellfun('length',is_in_index)==1);
%indices of these loops in list 2
indices_in_b=cell2mat(is_in_index(cellfun('length',is_in_index)==1));

%which of the loops have the same sign
ind_a_id=indices_in_a(loop_list_a.sign(indices_in_a)==...
    loop_list_b.sign(indices_in_b));
ind_b_id=indices_in_b(loop_list_a.sign(indices_in_a)==...
    loop_list_b.sign(indices_in_b));

%which of the loops have switching signs
ind_a_switch=indices_in_a(not(loop_list_a.sign(indices_in_a)==...
    loop_list_b.sign(indices_in_b)));
ind_b_switch=indices_in_b(not(loop_list_a.sign(indices_in_a)==...
    loop_list_b.sign(indices_in_b)));

%indices of loops in list a that do not occur in list b
ind_a_notin = find(cellfun('length',is_in_index)==0);

end