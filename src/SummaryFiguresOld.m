% compute statistics and create figures
delta = 3;
T = readtable('SummaryNewAnalysis.csv')
%%
allT = table2array(T);

%%
%t = unique(allT(:,2))
% col = 'cm'
% emb = unique(allT(:,1))
% figure
% for dv =0:1
% for id = 1:length(emb)
%     emb(id)
%     id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
%     plot(allT(id,2), allT(id,3), strcat(col(dv+1), '-')); hold on
% end
% end
% %%
% figure
% for dv =0:1
% for id = 1:length(emb)
%     emb(id)
%     id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
%     plot(allT(id,2), (allT(id,6)-allT(id,5))/(2*3), strcat(col(dv+1), '-')); hold on
% end
% end
%%
figure
for dv =0:1
for id = 1:length(emb)
    emb(id)
    id = find((allT(:,1)==emb(id)) &(allT(:,end)==dv));
    plot(allT(id,2), (allT(id,4)-allT(id,3))/(2*3), strcat(col(dv+1), '-')); hold on
end
end
xlabel('Time (frames)')
ylabel('Average displacements(pixels)')
%%
allT2 = allT;
allT2(:,3:6) = allT(:,3:6)/delta;
writetable(array2table(allT2, 'VariableNames',{'Embryo','Time', 'mean_all_left','mean_all_right', 'mean_sel_left','mean_sel_right','AP-DV' }), 'SummaryNewAnalysis.csv')
%%
tm =  unique(allT2(:,2));
for t = 1:length(tm)
    id = find(allT2(:,2)==tm(t));
    figure; boxplot((allT2(id,4)-allT2(id,3))/2,allT2(id,7));
    vals = (allT2(id,6)-allT2(id,5))/2;
    [h p] = ttest2(vals(find(allT2(id,7)==1)),vals(find(allT2(id,7)==0)))
    title(strcat('p = ', num2str(p), '. Time = ', num2str(tm(t))));
    
end
