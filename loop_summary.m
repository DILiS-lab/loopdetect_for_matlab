function[tab]=loop_summary(fbl_loop_tab,col_val)
% LOOP_SUMMARY Compute the counts of a feedback loop list
%
% The loop list could be output e.g. from find_loops() and should be a
% table with columns named 'length' and 'sign' that are used to determine
% the counts.
%
% TAB = LOOP_SUMMARY(LOOP_LIST) returns a table TAB with the counts of the loops
% in  LOOP_LIST subdivided in three rows capturing
% the sign of the loops (all, positive, negative) and their
% lengths (1 to maximal length) in the columns
%
% TAB = LOOP_SUMMARY(LOOP_LIST,COL_VAL) allows for choosing the
% value that is spread over the columns. For COL_VAL being 'length', the
% same output as above is obtained, for COL_VAL being 'sign', the table is
% transposed.
%
% See also: find_loops(), find_loops_noscc(), find_loops_vset()

if nargin<2
    col_val='length';
end

%maximal loop length
max_loop_length=max(fbl_loop_tab.length);
tab=zeros(3,max_loop_length);


%for each loop length
for i=1:max_loop_length
    tab(2,i)=sum(fbl_loop_tab.sign(fbl_loop_tab.length==i)==-1); %negative loops of length i
    tab(3,i)=sum(fbl_loop_tab.sign(fbl_loop_tab.length==i)==1);%positive loops of length i
    tab(1,i)=sum(fbl_loop_tab.length==i); %total loops of length i
end
%generate length description
for k = 1:max_loop_length
    col_name_tab{k} = sprintf('%s_%d','length',k);
end

if strcmp(col_val,'length') %variables will be length
    tab=array2table(int64(tab),'VariableNames',col_name_tab,'RowNames',{'total' 'negative' 'positive'});
end
if strcmp(col_val,'sign') %variable will be sign
    tab=array2table(int64(tab'),'VariableNames',{'total' 'negative' 'positive'},'RowNames',col_name_tab);
end

end
    