load out3.mat;
out_1=out;%lambdac=100
load out3_10.mat;
out_1_10=out;
load out3_200.mat;
out_1_200=out;

load out4.mat;
out_2=out;
load out5.mat;
out_3=out;
% est5=out.est_ms_member;
% est5_2=out.est_ms_member2;

%% ******************************************************************************************************************************%
%比较不同传感器数下的OSPA
truth.K=100;
mCnt=length(out_1.est_ms_member);
sensors=3;
ospa_vals= zeros(truth.K,sensors,mCnt);
ospa_vals1= zeros(truth.K,sensors,mCnt);
ospa_vals_Graph= zeros(truth.K,sensors,mCnt);
ospa_c= 100;
ospa_p= 1;
for ii=1:mCnt
    ospa_vals(:,1,ii)=out_1.est_ms_member{ii}.ospa';
    ospa_vals(:,2,ii)=out_2.est_ms_member{ii}.ospa';
    ospa_vals(:,3,ii)=out_3.est_ms_member{ii}.ospa';
    
    ospa_vals1(:,1,ii)=out_1.est_ms_member1{ii}.ospa';
    ospa_vals1(:,2,ii)=out_2.est_ms_member1{ii}.ospa';
    ospa_vals1(:,3,ii)=out_3.est_ms_member1{ii}.ospa';
    
    ospa_vals_Graph(:,1,ii)=out_1.est_ms_member2{ii}.ospa';
    ospa_vals_Graph(:,2,ii)=out_2.est_ms_member2{ii}.ospa';
    ospa_vals_Graph(:,3,ii)=out_3.est_ms_member2{ii}.ospa';
end

ospa=sum(ospa_vals,3)/mCnt;
ospa1=sum(ospa_vals1,3)/mCnt;
ospa_Graph=sum(ospa_vals_Graph,3)/mCnt;

%画图
nb_methods = 3;

X = []; X1 = [];Y = [];  ACTid = [];

for s_idx = 1:sensors
    X   = [X; ospa(:,s_idx)];
    X1 = [X1; ospa1(:,s_idx)];
    Y   = [Y; ospa_Graph(:,s_idx)];
    
    ACTid =  [ACTid; s_idx*ones(truth.K,1)];
end

% my_labels = {'                 S=3',' ','                  S=4',' ','                  S=5','  '};
% my_labels = {'                        S=3','','','                        S=4','','','                        S=5','',''};
my_labels = {'','{\itS}=3','','','{\itS}=4','','','{\itS}=5','',''};

xylabel = repmat('xyz', sensors*truth.K,1);
hfig = figure;
h_box = boxplot([X;X1; Y], {repmat(ACTid,nb_methods,1),  xylabel(:)} ,'factorgap',  0, 'Widths', 0.6,  'FactorSeparator',[1]);
box on
set(gca,'ygrid','on') 

set(gca,'xtick',1:14, 'xticklabel', my_labels) 
a = get(gca,'XTickLabel');
set(gca,'xTickLabel', a,'fontsize', 35)

color = ['g','b', 'c','g','b', 'c','g','b', 'c','g','b', 'c', 'g','b', 'c'];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

legend(c(1:3), 'Original MS-MeMBer filter ','MS-MeMBer filter with group structure estimation', 'Proposed method');
ylabel('Average OSPA distance');
for ii=1:size(h_box,2) %  for each box in the box plot 
     set(h_box(1:4,ii),'linewidth',1);  % 1:4 are the wiskers
     set(h_box(6,ii),'linewidth', 2.5, 'color','k');  % 6 is the median       
end

%% ******************************************************************************************************************************%
%比较不同杂波密度下的OSPA
truth.K=100;
mCnt=length(out_1.est_ms_member);
nlambda=3;
ospa_vals= zeros(truth.K,nlambda,mCnt);
ospa_vals_Graph= zeros(truth.K,nlambda,mCnt);
ospa_c= 100;
ospa_p= 1;
for ii=1:mCnt
    ospa_vals(:,1,ii)=out_1_10.est_ms_member{ii}.ospa';
    ospa_vals(:,2,ii)=out_1.est_ms_member{ii}.ospa';
    ospa_vals(:,3,ii)=out_1_200.est_ms_member{ii}.ospa';
    
    ospa_vals1(:,1,ii)=out_1_10.est_ms_member1{ii}.ospa';
    ospa_vals1(:,2,ii)=out_1.est_ms_member1{ii}.ospa';
    ospa_vals1(:,3,ii)=out_1_200.est_ms_member1{ii}.ospa';
    
    ospa_vals_Graph(:,1,ii)=out_1_10.est_ms_member2{ii}.ospa';
    ospa_vals_Graph(:,2,ii)=out_1.est_ms_member2{ii}.ospa';
    ospa_vals_Graph(:,3,ii)=out_1_200.est_ms_member2{ii}.ospa';
end

ospa=sum(ospa_vals,3)/mCnt;
ospa1=sum(ospa_vals1,3)/mCnt;
ospa_Graph=sum(ospa_vals_Graph,3)/mCnt;

%画图
nb_methods = 3;

X = [];X1 = []; Y = [];  ACTid = [];

for s_idx = 1:nlambda
    X   = [X; ospa(:,s_idx)];
    X1   = [X1; ospa1(:,s_idx)];
    Y   = [Y; ospa_Graph(:,s_idx)];
    
    ACTid =  [ACTid; s_idx*ones(truth.K,1)];
end

% my_labels = {'                        \lambda_c=10','','','                          \lambda_c=100','','','                          \lambda_c=200','',''};
my_labels = {'','\lambda_c=10','','','\lambda_c=100','','','\lambda_c=200','',''};

xylabel = repmat('xyz', nlambda*truth.K,1);
hfig = figure;
h_box = boxplot([X; X1; Y], {repmat(ACTid,nb_methods,1),  xylabel(:)} ,'factorgap',  0, 'Widths', 0.6,  'FactorSeparator',[1]);
box on
set(gca,'ygrid','on') 

set(gca,'xtick',1:14, 'xticklabel', my_labels) 
a = get(gca,'XTickLabel');
set(gca,'xTickLabel', a,'fontsize', 35)

color = ['g','b', 'c','g','b', 'c','g','b', 'c','g','b', 'c', 'g','b', 'c'];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

legend(c(1:3), 'Original MS-MeMBer filter ','MS-MeMBer filter with group structure estimation', 'Proposed method');
ylabel('Average OSPA distance');
for ii=1:size(h_box,2) %  for each box in the box plot 
     set(h_box(1:4,ii),'linewidth',1);  % 1:4 are the wiskers
     set(h_box(6,ii),'linewidth', 2.5, 'color','k');  % 6 is the median       
end



