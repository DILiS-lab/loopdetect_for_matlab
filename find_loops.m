function [loop_tab] = find_loops(A_in,max_num_loops)
%FIND_LOOPS Compute feedback loops from Jacobian or adjacency matrix.
%
% loop_tab = FIND_LOOPS(A_in, max_num_loops) returns a Matlab
% table with three columns, 'loop', 'length', 'sign', that represents 
% feedback loops as derived from matrix A_in.
% 'loop' is a cell array returning the regulations, e.g. the output 
% 1 2 3 1 means that variable 1 regulates variable 2, variable 2 regulates 
% variable 3 and variable 3 regulates variable 1.
% 'length' is a positive integer indicating the number of variables forming
% the loop, i.e. the loop length.
% 'sign' is a value in {-1,1} and indicates whether the feedback loop is 
% negative or positive, i.e. whether the number of negative interactions
% forming the loop is uneven or even.
% A_in is a (Jacobian) nxn-matrix that captures the regulations of n 
% variables on each other. Thereby, entry a_ij represents the regulation of 
% variable i by variable j. a_ij should be non-zero only if variable j 
% regulates variable i. The effect of the regulation (inhibition vs.
% activation) should be encoded in the sign of a_ij. 
% max_num_loops is an optional parameter indicating the maximal number of
% feedback loops that are detected and returned. To limit runtime, it 
% defaults to 1e6.
% 
% FIND_LOOPS uses the decomposition into strongly connected components 
% by Tarjan's algorithm (function tarjan.m from Matlab File Exchange server)
% and furthermore relies on Johnson's algorithm as adapted from 
% find_elem_circuits.m by cmaes on github.
%
% Examples
% System of 4 variables with 4 self-loops (loops of length 1), all
% negative, 1 positive loop of length 3 and one negative loop of length 4
%   J = [-1 0 0 -1; 1 -1 0 1; 0 1 -1 0; 0 0 1 -1];
%   loop_t = find_loops(J)
%   loop_t.loop{1} % delivers the first loop
% Limiting the output to 4 loops:
%   loop_t = find_loops(J,4)
%
% $Author: kabaum $  $Date: 2020/04/06 09:46 $ $Revision: 0.1 $
% Copyright: GNU GPLv3 K. Baum 2020
%
% See also: FIND_LOOPS_NOSCC(), find_loops_vset()


%set maximal number of reported loops (plus self-loops)
if nargin < 2
    max_num_loops = 1e6;
end

%transpose (required for correct direction detection) and turn into sparse
%representation
%i.e. entry A_in(2,3) determines how species 2 is regulated by species 3
if ~issparse(A_in)
    A = sparse(A_in'); %this command alters the format of how the matrix is represented
else
    A=A_in';
end

%preparatory functions for the matrix

    function[self_loop_tab,A_red] = determine_and_remove_self_loops(C)
        %determine self-loops, output as table
        self_loops=sign(diag(C));
        self_loops_ind=find(~(self_loops==0));
        self_loops=self_loops(self_loops_ind);
        A_red=C-diag(diag(C)); %remove diagonal entries
        self_loop_tab=table(... %write the input into result table
            mat2cell(cat(2,self_loops_ind,self_loops_ind),ones(length(self_loops_ind),1)),...
            ones(length(self_loops_ind),1),...
            self_loops,...
            'VariableNames',{'loop','length','sign'});
    end

%detect nonviable_nodes and set the matrix rows & columns to zero
    function[A_red]=set_nonviable_nodes_to_zero(C)
        %determine those variables that have no ingoing or no outgoing
        %edges!
        %assume that the diagonal entries have been set to zero already!
        nonviable_nodes=find(any(C)+any(C')<2);
        A_red=C;
        A_red(nonviable_nodes,:)=0; %set nonviable rows to zero
        A_red(:,nonviable_nodes)=0; %set nonviable columns to zero
    end

%preprocess A
[self_loop_tab,A] = determine_and_remove_self_loops(A); %determine self-loops and eradicate diag entries
[A]=set_nonviable_nodes_to_zero(A);

%now advance with determining strongly connected components 
%function from Matlab Exchange, see licence
[comp_ids]=tarjan(A); %

%size of largest connected component will be used to define the blocking
%structures
max_comp_size=length(comp_ids{1}); %they are sorted by descending size, this is therefore the maximal size

%functions for Johnson's algorithm (cmaes!) of simple circuit detection
%require global variables stack, blocked, Blist, cycles
n = max_comp_size;
Blist = cell(n,1); %this has something to do with blocking nodes
blocked = false(1,n); %whether nodes are blocked or not
stack=[]; %important for keeping track
cycles = {}; %save the cycles in there, the function circuit alters this object
numcycles=size(self_loop_tab,1); %keep track of found number of cycles, add detected self loops!

    function unblock(u)
        blocked(u) = false;
        for w=Blist{u}
            if blocked(w)
                unblock(w)
            end
        end
        Blist{u} = [];
    end

    function f = circuit(v, s, C,ind_transfer,num_loops_max) %this function seems to detect paths from v to s in matrix C, saves them to "cycle"
        f = false;
        
        stack(end+1) = ind_transfer(v);
        blocked(v) = true;
        
        for w=find(C(v,:)) %determine non-zero elements in row v of C
            if w == s
                cycles{end+1} = [stack ind_transfer(s)];
                numcycles=numcycles+1; %keep track of number of cycles
                f = true;
            elseif ~blocked(w)
                if circuit(w, s, C, ind_transfer,num_loops_max) % if there is a path from w to s
                    f = true;
                end
            end
            if numcycles>=num_loops_max %leave function if we have enough FBLs detected
                return
            end
        end
        
        if f
            unblock(v)
        else
            for w = find(C(v,:))
                if ~ismember(v, Blist{w})
                    Bnode = Blist{w};
                    Blist{w} = [Bnode v];
                end
            end
        end
        
        stack(end) = [];
    end


%go through every component of length larger than 1 - they are sorted in
%descending order of their sizes!
i=1;
while i<=length(comp_ids)&&length(comp_ids{i})>1  
    
    A_temp=A(comp_ids{i},comp_ids{i}); %submatrix of the strongest connected component variables only
    
    s = 2; %start at node 2; we might go for starting at node n and revert the direction?
    while s <= size(A_temp,1) %go from s=2 up to length of A_temp
        
        %we could also go and find connected components for each of these
        %subgraphs! try in terms of runtime...
    
        % Subgraph of G induced by {1,...,s,}
        F = A_temp(1:s,1:s);
        blocked(1:s) = false;
        Blist(1:s) = cell(s,1);
        circuit(s, s, F,comp_ids{i},max_num_loops);
        if numcycles >= max_num_loops %leave the loop if enough cycles have been detected
            break
        end
        s = s + 1; %go one species more

    end
    if numcycles >= max_num_loops %leave the loop if enough cycles have been detected
        break
    end
i=i+1; %go one connected component further
end

%now we have all cycles or the maximally allowed number
%start preparing the table

%determine cycle signs

loop_length=cellfun(@length,cycles)-1;
function[sign_val]=loop_sign_func(cycinds, C)
    %determine sign from indices
    sign_val= mod(sum(sign(C(sub2ind(size(C),cycinds(1:(end-1)),cycinds(2:end))))<1),2)...
        *(-2)+1;
end
loop_sign=cellfun(@(x)loop_sign_func(x,A),cycles);

%cast the loop table and add self-loops such that at most max_num_loops are
%returned
loop_tab=cat(1,table(cycles',loop_length',loop_sign','VariableNames',{'loop','length','sign'}),self_loop_tab(1:min(height(self_loop_tab),max_num_loops-numcycles),:));

end
        
       