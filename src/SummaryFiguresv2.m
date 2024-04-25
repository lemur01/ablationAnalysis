
% compute statistics and create figure
delta = 3
T = readtable('C:\Projects\Elena\Data\Results\DV1.csv')
%%
allT = table2array(T(:,2:6));
%%
%t = unique(allT(:,2))
col = 'cm'
emb = unique(allT(:,1))

% figure
% for dv =0:1
% for id = 1:length(emb)
%     emb(id)
%     id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
%     plot(allT(id,2), (allT(id,4)-allT(id,3))/(2*3), strcat(col(dv+1), '-')); hold on
% end
% end
%%
figure
for dv =0:1
for id = 1:length(emb)
    emb(id)
    id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
    plot(allT(id,2), (allT(id,5)-allT(id,4))/(2*delta), strcat(col(dv+1), '-')); hold on
end
end
xlabel('Time (frames)')
ylabel('Average displacements(pixels)')
%%
%%
tm =  unique(allT2(:,2));
for t = 1:length(tm)
    id = find(allT2(:,2)==tm(t));
    figure; boxplot((allT2(id,5)-allT2(id,4))/2,allT2(id,7));
    vals = (allT2(id,5)-allT2(id,4))/2;
    [h p] = ttest2(vals(find(allT2(id,7)==1)),vals(find(allT2(id,7)==0)))
    title(strcat('p = ', num2str(p), '. Time = ', num2str(tm(t))));
end    

