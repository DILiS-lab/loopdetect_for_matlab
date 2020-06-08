function[loop_rep,loop_rep_index,jac_rep,jac_rep_index]=find_loops_vset(func, vset, max_loop_num,compute_full_list)
%FIND_LOOPS_VSET Detect feedback loops from sets of variable values of
%interest.
%
% LOOP_REP = FIND_LOOPS_VSET(FUNC, VSET, MAX_LOOP_NUM) can be used as a wrapper for the
% functions numerical_jacobian_complex() and find_loops(). For the ODE function y'=FUNC(y),
% it returns the loop list at column vector of variable values y=VSET. 
% FUNC should be a function handle that depends on the variable values
% only; parameters and time have to be set to fixed values. 
% MAX_LOOP_NUM is an optional input parameter that restricts the maximal
% number of reported loops (default: 1e6). LOOP_REP is a cell array of length
% 1, LOOP_REP{1} returns the loop list as Matlab Table (formatted as in find_loops()).

% [LOOP_REP,LOOP_REP_INDEX]= FIND_LOOPS_VSET(FUNC, VSET)
% Multiple sets of variable values can be supplied to the function which reports
% loop lists corresponding to different classes of Jacobian matrices (those
% with similar sign structure). VSET can be, for example, the transposed 
% solution vector when solving the ODE system with
% one of the Matlab solvers, or concatenations of other sets of
% variable values. For an ODE system of size _n_, the dimension of VSET should
% be _nxk_ with _k_ being the number of different sets of variable values at
% which loops should be detected (i.e. k column vectors of length n are concatenated).
% LOOP_REP_INDEX has length k (one entry for each set of variable values in VSET)
% and captures as an integer index which loop list corresponds to which set
% of variable values. Please note that loop lists can be identical despite
% belonging to different sign classes of Jacobians.
%
% [LOOP_REP,LOOP_REP_INDEX]= FIND_LOOPS_VSET(FUNC, VSET,MAX_LOOP_NUM,COMPUTE_FULL_LIST)
% Syntax as before. If the logical input COMPUTE_FULL_LIST is set to false (default: true),
% the classes of Jacobians are further reduced and loop lists are not re-computed
% for Jacobians that clearly do not allow for altered loop lists. This is 
% the case if no new regulation appear and
% only signs of regulations are altered that are not member of any loop. 
% Loop lists can still be identical.
%
% [LOOP_REP,LOOP_REP_INDEX,JAC_REP,JAC_REP_INDEX]= FIND_LOOPS_VSET(_)
% Syntax as before. JAC_REP returns the different detected signed Jacobians,
% JAC_REP_INDEX gives the index of the signed Jacobian for each set of
% variable values. Only if COMPUTE_FULL_LIST is set to false, LOOP_REP can
% contain
% fewer elements than JAC_REP, otherwise both have the same number of elements. 
% 
% Examples
% Determine the loop list for a function at one set of variable values
%   klin=[165,0.044,0.27,550,5000,78,4.4,5.1];
%   knonlin=[0.3,2];
%   loop_list = find_loops_vset(@(x)func_POSm4(0,x,klin,knonlin),[1 2 1 1])
%
% Determine the loop list along the solution of an ODE system, report at 
% most 100000 loops and do not compute the loop list for all Jacobians.
%   sol=dlmread('li08_solution.tsv','delimiter','\t');
%   [loop_rep,loop_rep_index,jac_rep,jac_rep_index]=find_loops_vset(...
%       @(x)func(0,x), sol(2:end,:), 1e5, false)
% See also: FIND_LOOPS(), NUMERICAL_JACOBIAN_COMPLEX()


% Please note that reported loop lists might be identical if two sign
% switches occur that are both affecting always the same loops. 



if nargin<3
    max_loop_num=1e6;
end
if nargin<4
    compute_full_list=true;
end

%determine the initial Jacobian and loop list
J_ini=numerical_jacobian_complex(func,vset(:,1));

vset_count=size(vset,2); %how many different sets of variable values

%determine all different Jacobian matrices (signs only)
Jsig_tab=zeros(1,size(vset,2));
Jsig_rep={sign(J_ini)};
for i=1:vset_count
    J_temp=numerical_jacobian_complex(func,vset(:,i),1e-8);
    %check whether the new Jacobian is different from the Jacobians found earlier    
    for j=1:length(Jsig_rep) 
        disttemp=max(max(abs(sign(J_temp)-Jsig_rep{j})));
        if disttemp==0 %if another Jacobian is the same
            Jsig_tab(i)=j;
            break
        end
    end
    if Jsig_tab(i)==0 %if no other Jacobian was the same before
        Jsig_rep{end+1}=sign(J_temp); %save the Jacobian
        Jsig_tab(i)=length(Jsig_rep);
    end         
end
jac_rep=Jsig_rep;
jac_rep_index=Jsig_tab;

%determine loops for the different Jacobians
%loops should only differ if 
% - switch from zero no nonzero
% - switch of an entry that is an edge of a loop

loop_tab=zeros(1,size(vset,2));
loop_tab(Jsig_tab==1)=1; %set the loop index according to the Jacobian that belongs
loop_rep={find_loops(J_ini,max_loop_num)}; %J_ini equals Jsig_rep{1}

if compute_full_list %if we should compute the full loop list
    
    for i=2:length(Jsig_rep) %check each Jacobian (if there are more than 1)
        J_temp=Jsig_rep{i};
        loop_rep{end+1}=find_loops(J_temp,max_loop_num);
    end
    loop_tab=jac_rep_index; %loop identity equals Jacobian identity

else %if we should not compute the full loop list, more efficient (exclude Jacobians that cannot generate different loop lists)
  
    for i=2:length(Jsig_rep) %check each Jacobian (if there are more than 1)
        J_temp=Jsig_rep{i};
        %determine if there is a switch from zero to nonzero
        switch_to_nonzero=zeros(1,i-1);
        for j=1:(i-1) %compare with all earlier Jacobians
            switch_to_nonzero(j)=sum(sum(abs(J_temp(Jsig_rep{j}==0))))>0;
        end
        %determine if the sign switches of an edge that is contained in a loop
        change_in_a_loop=zeros(1,i-1);
        for j=1:(i-1) %compare with all earlier Jacobians
            if switch_to_nonzero(j)==0 %we only have to check this if we don't already have another change
                [target_node_ind, source_node_ind]=find(not(J_temp-Jsig_rep{j}==0));
                for k=1:length(target_node_ind)
                    change_in_a_loop(j)=not(isempty(find_edges(loop_rep{loop_tab(find(Jsig_tab==j,1))},... %the loop list corresponding to Jacobian j
                        source_node_ind(k),target_node_ind(k))));
                    if change_in_a_loop(j)==1 %in case there is a change in at least one edge
                        break
                    end
                end
            end
        end
        %check if there is any earlier Jacobian that is neither changed from zero to
        %nonzero for J_temp, nor any change in an edge of a loop is observed --> this
        %Jacobian will lead to exactly the same loops as the Jacobian of
        %interest, J_temp
        loop_ind_temp=find(switch_to_nonzero==0 & change_in_a_loop==0,1); %find only the first Jacobian that has the same profile
        if isempty(loop_ind_temp) %in case no Jacobian coincides, we compute the loops anew
            loop_rep{end+1}=find_loops(J_temp,max_loop_num);
            loop_tab(Jsig_tab==i)=length(loop_rep);
        else %in case the Jacobian is similar to another Jacobian 
            %(e.g. because the sign switch appeared in a nonzero entry that
            %does not generate any loop) we just use that Jacobian
            loop_tab(Jsig_tab==i)=loop_ind_temp;
        end   
    end
end

loop_rep_index=loop_tab;

end
        
