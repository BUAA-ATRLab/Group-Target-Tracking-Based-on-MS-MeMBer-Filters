function fHand = plot_results(model, truth, est, varargin)
% fHand = plot_results(model, truth, est, varargin)
% INPUT:
%      model - struc that contains kinematic and obs. model params. 
%      truth - struc that cointains groungtruth tracks
%      est   - struc that contains estimated tracks
%      varargin{1} - (optional) flag to plot track labels (default is true)
% OUTPUT:
%      fHand - figure handle
%
% Function used to plot estimated tracks on top of ground-truth tracks.
%
% Contact: augustin.saucan@mail.mcgill.ca

label_plot_flag = 1;
if ~isempty(varargin)
    label_plot_flag = varargin{1};
end

% verify that there are labels to plot
if ~isfield(est, 'state_l_hat')
    label_plot_flag = 0;
end

fHand = figure;
subplot(2,1,1) 
for t = 1:truth.K
    if ~isempty(truth.X{t})
    ht = plot(truth.X{t}(1,:), truth.X{t}(2,:),'dk','LineWidth',2, 'MarkerSize', 10);
    end
    hold on
    if ~isempty(est.state_X_hat{t})
        he = plot(est.state_X_hat{t}(1,:), est.state_X_hat{t}(2,:), '+r', 'LineWidth', 2, 'MarkerSize', 10);
        if label_plot_flag 
            txt_h = labelpoints( est.state_X_hat{t}(1,:), est.state_X_hat{t}(2,:), convert_label(est.state_l_hat{t}), 'N', 0.05, 1, 'FontSize', 10, 'Color', 'k' );
        end
    end
end
    axis([model.obs.xrange model.obs.yrange])
    ylabel('Y coordinate [m]', 'fontSize',25,'fontWeight','bold')
    xlabel('X coordinate [m]', 'fontSize',25,'fontWeight','bold')
    legend([ht he],{'True tracks', 'Estimated tracks'});

    grid on; box on; set(gca,'GridLineStyle','--');
    title(strcat('True target tracks and estimated tracks-', est.filter.title))
    set(gca,'FontSize',25,'fontWeight','bold');
subplot(2,1,2)   
    plot(1:truth.K, truth.N,'-k','LineWidth',3, 'MarkerSize', 15)
    hold on
    plot(1:truth.K, est.state_n_hat,'+r','LineWidth',3, 'MarkerSize', 10)
    axis([0 truth.K 1 15])
    ylabel('True and estimated cardinality','fontSize',25,'fontWeight','bold')
    xlabel('Time sample [s]','fontSize',25,'fontWeight','bold')
    legend('True cardinality', 'Estimated cardinality')
    grid on; box on; set(gca,'GridLineStyle','--');
    set(gca,'FontSize',25,'fontWeight','bold');
% set(findall(fHand,'type','text'),'fontSize',25,'fontWeight','bold')
% set(txt_h, 'fontSize',10, 'fontWeight', 'bold')


figure;
num_group=zeros(1,truth.K);
for t = 1:truth.K
    if ~isempty(est.graph_G_hat{t})
        num_group(t)=max(est.graph_G_hat{t});
    end
end
plot(1:truth.K, num_group,'-r','LineWidth',3, 'MarkerSize', 15);hold on
plot(1:truth.K, truth.G,'-k','LineWidth',3, 'MarkerSize', 15);hold on


%plot meas
figure; 
for k=1:truth.K
    if ~isempty(est.state_X_hat{k})
        plot(est.state_X_hat{k}(1,:), est.state_X_hat{k}(2,:), 'k.','MarkerSize', 4);hold on
    end
end
axis equal; axis([model.obs.range_c(1,1) model.obs.range_c(1,2) model.obs.range_c(2,1) model.obs.range_c(2,2)]); 
xlabel('x/m');
ylabel('y/m');

end


function label_cell = convert_label(label_array) 
% Function that converts the label matrix to a cell array of labels. Each
% cell contains the label of a target as a string of the type 'time_of_birth | birth_index'.  

label_cell = cell(1,length(label_array(1,:)));
for i = 1:length(label_array(1,:))
    label_cell{i} = strcat( num2str(label_array(1,i)), '|', num2str(label_array(2,i)));
end

end