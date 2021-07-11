function truth= gen_truth(model)

%variables
truth.K= 100;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.G= zeros(truth.K,1);            %ground truth for number of groups
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 11;
wturn = 2*pi/180;

xstart(:,1)  = [ -200+30; -10; 2000-20; 0; 0 ];        tbirth(1)  = 1;     tdeath(1)  = 80;
xstart(:,2)  = [ -200+30; -10; 2000+20; 0; 0 ];          tbirth(2)  = 1;    tdeath(2)  = 80;
xstart(:,3)  = [ -200; -10; 2000; 0; 0 ];           tbirth(3)  = 1;    tdeath(3)  = 80;
xstart(:,4)  = [ -200+30; -10; 2000+60; 0; 0 ];           tbirth(4)  = 1;    tdeath(4)  = 80;

xstart(:,5)  = [ 0; 10; 500; 6; 0 ];           tbirth(5)  = 1;    tdeath(5)  = truth.K+1;

xstart(:,6)  = [ -1000; 0; 300; 10; 0 ];           tbirth(6)  = 1;    tdeath(6)  = truth.K+1;
xstart(:,7)  = [ -1040; 0; 300; 10; 0 ];           tbirth(7)  = 1;    tdeath(7)  = truth.K+1;
xstart(:,8)  = [ -1020; 0; 340; 10; 0 ];           tbirth(8)  = 1;    tdeath(8)  = truth.K+1;

xstart(:,9)  = [ 300; 10; 2000; 0; 0 ];        tbirth(9)  = 21;     tdeath(9)  = truth.K+1;
xstart(:,10)  = [ 300+30; 10; 2000-30; 0; 0 ];          tbirth(10)  = 21;    tdeath(10)  = truth.K+1;
xstart(:,11)  = [ 300; 10; 2000-260; 10; 0 ];           tbirth(11)  = 21;    tdeath(11)  = truth.K+1;
%% generate the tracks

% % target 1 to 3 (group 1)
% targetnum=[1 2 3];
% targetstate = xstart(:,targetnum);
% for k=1:100
%     targetstate = gen_groupstate(targetstate, model);
%     truth.X{k}= [truth.X{k} targetstate];
%     truth.track_list{k} = [truth.track_list{k} targetnum];
%     truth.N(k) = truth.N(k) + 3;
% end
% 
% % target 4 (group 2)
% for targetnum=4
%     targetstate = xstart(:,targetnum);
%     for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
%         targetstate = gen_newstate(targetstate, model);
%         truth.X{k}= [truth.X{k} targetstate];
%         truth.track_list{k} = [truth.track_list{k} targetnum];
%         truth.N(k) = truth.N(k) + 1;
%      end
% end



for targetnum=[1 2 3]
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
     end
end

for targetnum=4
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):tbirth(targetnum)+60-1
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(4,:)=10;
    for k=tbirth(targetnum)+60:min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
     end
end


for targetnum=5
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
     end
end
for targetnum=6:8
    targetstate = xstart(:,targetnum);
    for k=1:10
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(5,:)=wturn/2;
    for k=11:30
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(5,:)=-wturn/2;
    for k=31:70
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(5,:)=wturn/2;
    for k=71:90
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(5,:)=0;
    for k=91:100
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
end

for targetnum=9:10
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
     end
end

for targetnum=11
    targetstate = xstart(:,targetnum);
    for k=tbirth(targetnum):tbirth(targetnum)+20-1
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
    end
    targetstate(4,:)=0;
    for k=tbirth(targetnum)+20:min(tdeath(targetnum),truth.K)
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        targetstate = gen_newstate(targetstate, model);
     end
end

truth.total_tracks= nbirths;

for k=1:truth.K
   if ~isempty(truth.X{k}) 
       truth.X{k}=truth.X{k}([1 3 2 4],:);
       P=repmat(diag([10 10 10 10]),1,1,truth.N(k));
       Adj = genere_groups_graph(truth.X{k}, P, 15, 15);
       Gk= connected_components(Adj);
       truth.G(k)=length(Gk);
   end
end





%plot ground truths
figure; truths= gcf; hold on;
[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 2],:);
    plot( Pt(1,:),Pt(2,:),'k-','MarkerSize',3,'LineWidth',1); 
    plot( Pt(1,1), Pt(2,1), 'ko','MarkerSize',4,'LineWidth',1);
    plot( Pt(1,(k_death(i)-k_birth(i)+1)), Pt(2,(k_death(i)-k_birth(i)+1)), 'k^','MarkerSize',4,'LineWidth',1);
end
axis equal;
axis([model.obs.range_c(1,1) model.obs.range_c(1,2) model.obs.range_c(2,1) model.obs.range_c(2,2)]); 
xlabel('x(m)');
% xlabel({'x(m)';'\o:start£»\Delta:end'});
ylabel('y(m)');
end

function [X_track,k_birth,k_death]= extract_tracks(X,track_list,total_tracks)

K= size(X,1); 
x_dim= size(X{K},1); 
k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end;
X_track= zeros(x_dim,K,total_tracks);
k_birth= zeros(total_tracks,1);
k_death= zeros(total_tracks,1);

max_idx= 0;
for k=1:K
    if ~isempty(X{k}),
        X_track(:,k,track_list{k})= X{k};
    end;
    if max(track_list{k})> max_idx, %new target born?
        idx= find(track_list{k}> max_idx);
        k_birth(track_list{k}(idx))= k;
    end;
    if ~isempty(track_list{k}), max_idx= max(track_list{k}); end;
    k_death(track_list{k})= k;
end;
end

function ca= makecolorarray(nlabels)
    lower= 0.1;
    upper= 0.9;
    rrr= rand(1,nlabels)*(upper-lower)+lower;
    ggg= rand(1,nlabels)*(upper-lower)+lower;
    bbb= rand(1,nlabels)*(upper-lower)+lower;
    ca.rgb= [rrr; ggg; bbb]';
    ca.lab= cell(nlabels,1);
    ca.cnt= 0;   
end