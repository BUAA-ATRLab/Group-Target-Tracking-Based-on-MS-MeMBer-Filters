%% plot OSPA

ospa=zeros(mcnt,truth.K);
ospa1=zeros(mcnt,truth.K);
ospa2=zeros(mcnt,truth.K);
for i=1:mcnt
   ospa2(i,:)= est_ms_member2{i}.ospa;
   ospa1(i,:)= est_ms_member1{i}.ospa;
   ospa(i,:)= est_ms_member{i}.ospa;
end
ave_ospa=mean(ospa,1);
ave_ospa1=mean(ospa1,1);
ave_ospa2=mean(ospa2,1);
figure;
plot(1:truth.K,ave_ospa2,'r','LineWidth',1);hold on
plot(1:truth.K,ave_ospa1,'b','LineWidth',1);hold on
plot(1:truth.K,ave_ospa,'k','LineWidth',1);hold on
legend('Proposed Method','MS-MeMBer filter with group structure estimation','Original MS-MeMBer filter','MarkerSize',15);
xlabel('time(s)');
ylabel('OSPA distance');

%% plot card
card2=zeros(mcnt,truth.K);
mean_card2=zeros(1,truth.K);
std_card2=zeros(1,truth.K);
card1=zeros(mcnt,truth.K);
mean_card1=zeros(1,truth.K);
std_card1=zeros(1,truth.K);
card=zeros(mcnt,truth.K);
mean_card=zeros(1,truth.K);
std_card=zeros(1,truth.K);
for k=1:mcnt
    card2(k,:)=est_ms_member2{k}.state_n_hat';
    card1(k,:)=est_ms_member1{k}.state_n_hat';
    card(k,:)=est_ms_member{k}.state_n_hat';
end
mean_card2=mean(card2);
std_card2=std(card2);
mean_card1=mean(card1);
std_card1=std(card1);
mean_card=mean(card);
std_card=std(card);
figure;
subplot(3,1,1);
stairs(1:truth.K,truth.N,'k','LineWidth',1);hold on
plot(1:truth.K,mean_card2,'r.','LineWidth',1);hold on
plot(1:truth.K,mean_card2+std_card2,'b-','LineWidth',1);hold on
plot(1:truth.K,mean_card2-std_card2,'b-','LineWidth',1);hold on
grid on;
legend('True','Mean','StDev');
title('Proposed Method');
xlabel('(a)');
axis([1 truth.K 5 13]);

subplot(3,1,2);
stairs(1:truth.K,truth.N,'k','LineWidth',1);hold on
plot(1:truth.K,mean_card1,'r.','LineWidth',1);hold on
plot(1:truth.K,mean_card1+std_card1,'b-','LineWidth',1);hold on
plot(1:truth.K,mean_card1-std_card1,'b-','LineWidth',1);hold on
grid on;
xlabel('(b)');
axis([1 truth.K 5 13]);
ylabel('The member cardinality');
title('MS-MeMBer filter with group structure estimation');

subplot(3,1,3);
stairs(1:truth.K,truth.N,'k','LineWidth',1);hold on
plot(1:truth.K,mean_card,'r.','LineWidth',1);hold on
plot(1:truth.K,mean_card+std_card,'b-','LineWidth',1);hold on
plot(1:truth.K,mean_card-std_card,'b-','LineWidth',1);hold on
grid on;
axis([1 truth.K 5 13]);
xlabel({'(c)','time(s)'});
title('Original MS-MeMBer filter');




