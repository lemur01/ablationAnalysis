
% compute statistics and create figure
delta = 3
T = readtable('..\data\Results\Summary.csv')
allT = table2array(T);
%%
col = 'cm'
emb = unique(allT(:,1))

figure
for dv =0:1
for id = 1:length(emb)
    emb(id)
    id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
    plot(allT(id,2), (allT(id,4)-allT(id,3))/(2*delta), strcat(col(dv+1), '-')); hold on
end
end
xlabel('Time (frames)')
ylabel('Average displacements(pixels)')
%%
% figure
% for dv =0:1
% for id = 1:length(emb)
%     emb(id)
%     id = find((allT(:,1)==emb(id)) & (allT(:,end)==dv));
%     plot(allT(id,2), (allT(id,6)-allT(id,5))/(2*delta), strcat(col(dv+1), '-')); hold on
% end
% end
% xlabel('Time (frames)')
% ylabel('Average displacements(pixels)')
%%
%%
tm =  unique(allT(:,2));
for t = 1:length(tm)
    id = find(allT(:,2)==tm(t));
    figure; boxplot((allT(id,4)-allT(id,3))/2,allT(id,7));
    vals = (allT(id,4)-allT(id,3))/2;
    [h p] = ttest2(vals(find(allT(id,7)==1)),vals(find(allT(id,7)==0)))
    title(strcat('T.test .p = ', num2str(p), '. Time = ', num2str(tm(t))));
    xticklabels({'AP', 'DV'})
    ylabel('displacement (pixels)')
    set(gca, 'FontSize',14)
end    
