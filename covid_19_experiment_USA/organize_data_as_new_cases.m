clear all, close all, clc
load('covid_19_confirmed_cases.mat');
indx2delete = find(Position(:,1) == 0);
Data(indx2delete,:) = [];
Position(indx2delete,:) = [];
%%
new_cases = zeros(size(Data));
new_cases(:,1) = Data(:,1);
for i=2:size(Data,2)
    new_cases(:,i) = Data(:,i)-Data(:,i-1);
end
Data = new_cases;
save('covid_19_new_cases.mat','Data','Position');