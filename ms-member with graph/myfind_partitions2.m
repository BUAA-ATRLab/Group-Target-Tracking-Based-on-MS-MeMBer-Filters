% S-D维分配算法求最优划分P,添加滑窗
function out = myfind_partitions2(model, filter, MS_MeMBer_predict, Z)   

%% 滑窗参数设置
window = filter.window_len;%窗长,1维假设+（window-1）维量测，window>=2
if window<2
   fprintf('\n 滑窗长度太小');
   exit
end
M_kk1 = MS_MeMBer_predict.M; % number of predicted Bernoulli components
win_init=1:window-1:model.obs.nb_sensors;
win_end=win_init+window-2;
win_end(win_end>model.obs.nb_sensors)=model.obs.nb_sensors;

%% 初始化
curr_paths = 1:M_kk1; % the paths of all the current subsets W
curr_paths=curr_paths';
% we initilize the path parameters with the predicted parameters 
curr_w  = log(MS_MeMBer_predict.r); % we assume the mass of the current component as initial subset weight 
curr_m  = MS_MeMBer_predict.m;
curr_P  = MS_MeMBer_predict.P;


%% 滑窗循环

for i=1:length(win_init)
    s_init=win_init(i);%窗口起始传感器标号
    s_end=win_end(i);%窗口结尾传感器标号
    %求得各种组合
    [curr_paths,curr_w,curr_m,curr_P]=MDA_window(Z,s_init,s_end,MS_MeMBer_predict,model, filter,curr_paths,curr_w,curr_m,curr_P);
        
    %MDA
    valcost=curr_w;
    minvalue=min(min(valcost));
    if isinf(minvalue)
        minvalue=-3000;
    end
    val=valcost-minvalue+1;
    [~,~,IC] = unique(curr_paths(:,1:s_init),'rows');
    curr_paths1=[IC curr_paths(:,s_init+1:end)];
    [gam1, ~] = k_best_assignment(curr_paths1, -val, filter.SD_max);
    temp_paths=[];
    for j=1:size(gam1,3)
        temp_paths=[temp_paths;gam1(:,:,j)];
    end
    [~,idx]=ismember(temp_paths,curr_paths1,'rows');
    curr_paths=curr_paths(idx,:);
    curr_w=curr_w(:,idx);
    curr_m=curr_m(:,idx);
    curr_P=curr_P(:,:,idx);    
end

%输出
log_alpha_P=zeros(1,size(gam1,3));
gam_temp=curr_paths(:,2:end);
idx0=find(gam_temp>M_kk1&gam_temp<=2*M_kk1);
gam_temp(idx0)=0;
idx1=find(gam_temp>2*M_kk1);
gam_temp(idx1)=gam_temp(idx1)-2*M_kk1;
curr_paths(:,2:end)=gam_temp;

for j=1:size(gam1,3);
    %% out
    out.sub{j}=curr_paths((j-1)*M_kk1+1:j*M_kk1,:);
    out.valCost{j}=curr_w(:,(j-1)*M_kk1+1:j*M_kk1);
    out.m{j}=curr_m(:,(j-1)*M_kk1+1:j*M_kk1);
    out.P{j}=curr_P(:,:,(j-1)*M_kk1+1:j*M_kk1);
    log_alpha_P(j)=sum(out.valCost{j});
end
out.alpha_P = exp( log_alpha_P - log( sum( exp(log_alpha_P) ) ) ); % normalize






end