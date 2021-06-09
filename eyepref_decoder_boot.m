function [Dprime_p_value_matrix,dec_acc_p_value_matrix,shuff_output] = eyepref_boot_comparison(x,y,condition_lables,showplots_dprime,showplots_dec_acc)
% x: predictor (i.e. Axon ratio)
% Y: data (i.e. AMPA ODI)
shuffle_n = 10000;
boot_n = 10000;

if size(y,1)>size(y,2)
    y = y';
end
pop_n = length(y);
if size(x,1)==pop_n
    x = x';
end

ybinary = y>0;

idx_contra = find(ybinary);
idx_ipsi = find(~ybinary);

% calculate measured dprime for each version of the predictor (e.g. Fr with different masks)
for ii = 1:size(x,1)
    Dprime(ii) = dprime(x(ii,idx_contra), x(ii,idx_ipsi));
    dec_acc(ii) = eye_dec_acc(x(ii,idx_contra), x(ii,idx_ipsi));
end

% get paired bootstraped sets of all predictor values
clear bootstrp_dprime
for i = 1:boot_n
    contra_samp = idx_contra(randi([1,length(idx_contra)],length(idx_contra),1));
    ipsi_samp = idx_ipsi(randi([1,length(idx_ipsi)],length(idx_ipsi),1));
    
    for ii = 1:size(x,1)
        Dprime_bootstrp(ii,i) = dprime(x(ii,contra_samp), x(ii,ipsi_samp));
        dec_acc_bootstrp(ii,i) = eye_dec_acc(x(ii,contra_samp), x(ii,ipsi_samp));
    end
end

% get the shuffled distribution for each predictor value
for i = 1:shuffle_n
    shuffle_idx = randperm(length(ybinary),length(ybinary));
    ybinary_shuf = ybinary(shuffle_idx);

    idx_contra_shuf = find(ybinary_shuf);
    idx_ipsi_shuf = find(~ybinary_shuf);

    for ii = 1:size(x,1)
        Dprime_shuf(ii,i) = dprime(x(ii,idx_contra_shuf), x(ii,idx_ipsi_shuf));
        dec_acc_shuf(ii,i) = eye_dec_acc(x(ii,idx_contra_shuf), x(ii,idx_ipsi_shuf));
    end
end

if size(x,1)>1
    H0 = 0; % nul hypothesis value
    for ii = 1:size(x,1) %% condition 1
        for jj = 1:size(x,1) %% condition 2
            if ii==jj
                % difference between shuffled data and observed value
                Dprime_shuf_diff = (Dprime_shuf(ii,:)-Dprime(ii));
                dec_acc_shuf_diff = (dec_acc_shuf(ii,:)-dec_acc(ii));
                
                % 2 sided test (number of differences greater than zero),(number of differences smaller than zero)
                % Dprime_p_value_matrix(ii,jj) = min([mean(Dprime_bott_diff>=H0), mean(Dprime_bott_diff<=H0)])*2;

                Dprime_p_value_matrix(ii,jj) = min([mean(Dprime_shuf_diff>=H0), mean(Dprime_shuf_diff<=H0)])*2;
                dec_acc_value_matrix(ii,jj) = min([mean(dec_acc_shuf_diff>=H0), mean(dec_acc_shuf_diff<=H0)])*2;
            else
                % difference between bootstraped data and observed value
                Dprime_bott_diff = (Dprime_bootstrp(ii,:)-Dprime_bootstrp(jj,:));
                dec_acc_bott_diff = (dec_acc_bootstrp(ii,:)-dec_acc_bootstrp(jj,:));
                
                % 2 sided test (number of differences greater than zero),(number of differences smaller than zero)
                % see  http://qed.econ.queensu.ca/working_papers/papers/qed_wp_1127.pdf for equation and explanation
                Dprime_p_value_matrix(ii,jj) = min([mean(Dprime_bott_diff>=H0), mean(Dprime_bott_diff<=H0)])*2;
                dec_acc_p_value_matrix(ii,jj) = min([mean(dec_acc_bott_diff>=H0), mean(dec_acc_bott_diff<=H0)])*2;
            end
        end
    end
else
    Dprime_p_value_matrix = nan;
    dec_acc_p_value_matrix = nan;
end

if showplots_dprime
    % summary
    figure; set(gcf,'color','w')
    if size(x,1)>1
        subplot(1,2,1)
    end
    for ii = 1:size(x,1)
        distributionPlot(Dprime_bootstrp(ii,:)','xValues',ii,'showMM',0,'histOpt',0,'color',[0.5 0.5 0.5]); hold on
        scatter(ii,mean(Dprime_bootstrp(ii,:)),'k','filled'); 
        errorbar(ii,mean(Dprime_bootstrp(ii,:)),std(Dprime_bootstrp(ii,:)),'k','LineWidth',2);
        scatter(ii,Dprime(ii),'^r','filled'); hold on
        
        distributionPlot(Dprime_shuf(ii,:)','xValues',ii,'showMM',0,'histOpt',0,'color',[0.5 0 0]);
        scatter(ii,mean(Dprime_shuf(ii,:)),'r','filled')
        errorbar(ii,mean(Dprime_shuf(ii,:)),std(Dprime_shuf(ii,:)),'r','LineWidth',2);
    end
    xlim([0 size(x,1)+1])
    xticks([1 : size(x,1)]);
    set(gca,'TickDir','out')
    ylabel('dprime'); ylim('auto')
    legend('off')
    xticklabels(condition_lables); xtickangle(45)
    axis square
    
    if size(x,1)>1
        ax = subplot(1,2,2);
        imagesc(Dprime_p_value_matrix);
        colormap(ax,flip(cbrewer('div','RdYlGn',100),1))
        xticks([1:size(x,1)])
        yticks([1:size(x,1)])
        caxis([0.001 0.05])
        axis image
        set(gca, 'XAxisLocation', 'top')
        set(gca, 'TickDir', 'out')
        colorbar
        title('pairwise P-values')
        xticklabels(condition_lables); xtickangle(45)
        yticklabels(condition_lables)
        suptitle('Dprime statistics')
    end
end
if showplots_dec_acc 
    % summary decoding accuracy 
    figure; set(gcf,'color','w')
    if size(x,1)>1
        subplot(1,2,1)
    end
    for ii = 1:size(x,1)
        distributionPlot(dec_acc_bootstrp(ii,:)'*100,'xValues',ii,'showMM',0,'histOpt',0,'color',[0.5 0.5 0.5]); hold on
        scatter(ii,mean(dec_acc_bootstrp(ii,:))*100,'k','filled'); 
        errorbar(ii,mean(dec_acc_bootstrp(ii,:))*100,std(dec_acc_bootstrp(ii,:))*100,'k','LineWidth',2);
        scatter(ii,dec_acc(ii)*100,'^r','filled'); hold on
        
        distributionPlot(dec_acc_shuf(ii,:)'*100,'xValues',ii,'showMM',0,'histOpt',0,'color',[0.5 0 0]);
        scatter(ii,mean(dec_acc_shuf(ii,:))*100,'r','filled')
        errorbar(ii,mean(dec_acc_shuf(ii,:))*100,std(dec_acc_shuf(ii,:))*100,'r','LineWidth',2);
    end
    xlim([0 size(x,1)+1])
    xticks([1 : size(x,1)]);
    set(gca,'TickDir','out')
    ylabel('decoding accuracy (%)'); ylim('auto')
    legend('off')
    xticklabels(condition_lables); xtickangle(45)
    axis square
    
    if size(x,1)>1
        ax = subplot(1,2,2);
        imagesc(dec_acc_p_value_matrix);
        colormap(ax,flip(cbrewer('div','RdYlGn',100),1))
        xticks([1:size(x,1)])
        yticks([1:size(x,1)])
        caxis([0.001 0.05])
        axis image
        set(gca, 'XAxisLocation', 'top')
        set(gca, 'TickDir', 'out')
        colorbar
        title('pairwise P-values')
        xticklabels(condition_lables); xtickangle(45)
        yticklabels(condition_lables)
        suptitle('decoding accuracy')
    end
end

shuff_output.Dprime = Dprime;
shuff_output.dec_acc = dec_acc;
shuff_output.Dprime_bootstrp = Dprime_bootstrp;
shuff_output.dec_acc_bootstrp = dec_acc_bootstrp;
shuff_output.Dprime_shuf = Dprime_shuf;
shuff_output.dec_acc_shuf = dec_acc_shuf;

end

function dec_ac = eye_dec_acc(x_contra,x_ipsi)
xmean_contra = mean(x_contra);
xstd_contra = std(x_contra);
xmean_ipsi = mean(x_ipsi);
xstd_ipsi = std(x_ipsi);

full_pop = [x_contra,x_ipsi]; 
eye_pref_pop = cat(2,ones(1,length(x_contra)),zeros(1,length(x_ipsi)));

for i = 1:length(full_pop)
    prob_contra = normpdf(full_pop(i),xmean_contra,xstd_contra);
    prob_ipsi = normpdf(full_pop(i),xmean_ipsi,xstd_ipsi);
    
    prediction(i) = prob_contra>prob_ipsi;
end

dec_ac = sum((prediction == eye_pref_pop))./length(full_pop); % percent predictions correct

end

function Dprime = dprime(contra, ipsi)
xmean_contra = nanmean(contra);
xmean_ipsi = nanmean(ipsi);
xstd_contra = nanstd(contra);
xstd_ipsi = nanstd(ipsi);
n_contra = length(contra);
n_ipsi = length(ipsi);
pooled_std = sqrt([(n_contra-1)*xstd_contra^2 + (n_ipsi-1)*xstd_ipsi^2]/[n_contra + n_ipsi - 2]);
Dprime = (xmean_contra-xmean_ipsi)/pooled_std;

end






