% Columns in 'data' matrices (other columns are blank):
% 1. Trial # for given participant
% 3. Certain points
% 4. Gain points
% 5. Loss points
% 7. Chosegamble (1/0)
% 8. Outcome (= Certain points, Gain points, or Loss points)
% 9. RT (seconds)
% 13. Time elapsed (from experiment start)
% 14. Updated total score (starts at 1000 at experiment start)
% 15. Trial outcome (1 = safe, 2 = loss, 3 = won)

load('risktaskdata_20220907.mat') %replace with actual data
set(0,'defaultfigurecolor',[1 1 1]); set(0,'defaultaxesfontsize',16);
set(0,'defaultlinelinewidth',2); set(0,'defaultaxeslinewidth',2);

NUM_PARTICIPANTS = length(alldata);
NUM_TRIALS = length(alldata(1).data);
ALPHA = 0.05;
INIT_SCORE = 1000;

for i = 1 : NUM_PARTICIPANTS

    individual_data = alldata(i).data;

    certain_trials = individual_data(individual_data(:,7) == 0);
    safe_ev = sum(individual_data(certain_trials,3));
    
    gamble_trials = individual_data(individual_data(:,7) == 1);
    gamble_ev = mean(sum(individual_data(gamble_trials,4:5)));
    
    alldata(i).final_ev = safe_ev + gamble_ev;

end

avg_ev = mean([alldata.final_ev]);
avg_final_score = INIT_SCORE + avg_ev

for i = 1 : NUM_PARTICIPANTS

    individual_data = alldata(i).data;   

    gamble_trials = individual_data(individual_data(:,7) == 1);
    won_gambles = sum(individual_data(gamble_trials,15) == 3);

    alldata(i).percent_gambles_won = 100 * won_gambles / length(gamble_trials);
    alldata(i).percent_gambles_chosen = 100 * length(gamble_trials) / NUM_TRIALS;

end

overall_avg_percent_gambles_won = mean([alldata.percent_gambles_won])

figure
scatter([alldata.percent_gambles_chosen],[alldata.percent_gambles_won], 'filled')
xlim([0 100])
ylim([0 100])
xlabel('% of gambles chosen')
ylabel('% of gambles won')
title('Percent of gambles chosen vs Percent of gambles won')
hold on
yline(50)
There seems to be no correlation between percentage of gambles chosen and percentage of gambles won. The percentages of gambles won generally hover around 50% regardless of the actor's attitude toward risk.

for i = 1 : NUM_PARTICIPANTS

    individual_data = alldata(i).data;

    gain_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) == 0);
    num_gain_gambles = sum(individual_data(gain_trials,7) == 1);
    alldata(i).percent_risk_in_gain = 100 * num_gain_gambles / length(gain_trials);

    loss_trials = individual_data(individual_data(:,4) == 0 & individual_data(:,5) < 0);
    num_loss_gambles = sum(individual_data(loss_trials,7) == 1);
    alldata(i).percent_risk_in_loss = 100 * num_loss_gambles / length(loss_trials);

    mixed_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) < 0);
    num_mixed_gambles = sum(individual_data(mixed_trials,7) == 1);
    alldata(i).percent_risk_in_mixed = 100 * num_mixed_gambles / length(mixed_trials);

end

figure
subplot(1,3,1)
scatter([alldata.percent_risk_in_gain], [alldata.percent_risk_in_mixed], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken, gain")
ylabel("% risks taken, mixed")
[r,p] = corr([alldata.percent_risk_in_gain]', [alldata.percent_risk_in_mixed]');
title("Gain vs Mixed")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))

subplot(1,3,2)
scatter([alldata.percent_risk_in_gain], [alldata.percent_risk_in_loss], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken, gain")
ylabel("% risks taken, loss")
[r,p] = corr([alldata.percent_risk_in_gain]', [alldata.percent_risk_in_loss]');
title("Gain vs Loss")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))

subplot(1,3,3)
scatter([alldata.percent_risk_in_mixed], [alldata.percent_risk_in_loss], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken, mixed")
ylabel("% risks taken, loss")
[r,p] = corr([alldata.percent_risk_in_mixed]', [alldata.percent_risk_in_loss]');
title("Mixed vs Loss")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))
In the gains vs mixed correlation (r = 0.395, p = 0.162), we fail to reject the null hypothesis, as p > .05.
In the gains vs loss correlation (r = 0.409, p = 0.147), we fail to reject the null hypothesis, as p > .05.
In the mixed vs loss correlation (r = 0.736, 0.003), we reject the null hypothesis, as p < .05.
The strongest correlation is in the final comparison between mixed and loss conditions. It is much weaker in the first two conditions, gains vs mixed and gains vs loss.

for i = 1 : NUM_PARTICIPANTS

    individual_data = alldata(i).data;

    gain_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) == 0);

    gain_trials_1 = gain_trials(1:length(gain_trials) / 2, :);
    num_gain_gambles_1 = sum(individual_data(gain_trials_1,7) == 1);
    alldata(i).percent_risk_in_gain_1 = 100 * num_gain_gambles_1 / length(gain_trials_1);
    
    gain_trials_2 = gain_trials(length(gain_trials) / 2 + 1 : end, :);
    num_gain_gambles_2 = sum(individual_data(gain_trials_2,7) == 1);
    alldata(i).percent_risk_in_gain_2 = 100 * num_gain_gambles_2 / length(gain_trials_2);

    loss_trials = individual_data(individual_data(:,4) == 0 & individual_data(:,5) < 0);

    loss_trials_1 = loss_trials(1:length(loss_trials) / 2, :);
    num_loss_gambles_1 = sum(individual_data(loss_trials_1,7) == 1);
    alldata(i).percent_risk_in_loss_1 = 100 * num_loss_gambles_1 / length(loss_trials_1);
    
    loss_trials_2 = loss_trials(length(loss_trials) / 2 + 1 : end, :);
    num_loss_gambles_2 = sum(individual_data(loss_trials_2,7) == 1);
    alldata(i).percent_risk_in_loss_2 = 100 * num_loss_gambles_2 / length(loss_trials_2);

    mixed_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) < 0);
    
    mixed_trials_1 = mixed_trials(1:length(mixed_trials) / 2, :);
    num_mixed_gambles_1 = sum(individual_data(mixed_trials_1,7) == 1);
    alldata(i).percent_risk_in_mixed_1 = 100 * num_mixed_gambles_1 / length(mixed_trials_1);
    
    mixed_trials_2 = mixed_trials(length(mixed_trials) / 2 + 1 : end, :);
    num_mixed_gambles_2 = sum(individual_data(mixed_trials_2,7) == 1);
    alldata(i).percent_risk_in_mixed_2 = 100 * num_mixed_gambles_2 / length(mixed_trials_2);

end

figure
subplot(1,3,1)
scatter([alldata.percent_risk_in_gain_1], [alldata.percent_risk_in_gain_2], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken in first half")
ylabel("% risks taken in second half")
[r,p] = corr([alldata.percent_risk_in_gain_1]', [alldata.percent_risk_in_gain_2]');
title("Gain")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))

subplot(1,3,2)
scatter([alldata.percent_risk_in_loss_1], [alldata.percent_risk_in_loss_2], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken in first half")
ylabel("% risks taken in second half")
[r,p] = corr([alldata.percent_risk_in_loss_1]', [alldata.percent_risk_in_loss_2]');
title("Loss")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))

subplot(1,3,3)
scatter([alldata.percent_risk_in_mixed_1], [alldata.percent_risk_in_mixed_2], 'filled')
xlim([0 100])
ylim([0 100])
xlabel("% risks taken in first half")
ylabel("% risks taken in second half")
[r,p] = corr([alldata.percent_risk_in_mixed_1]', [alldata.percent_risk_in_mixed_2]');
title("Mixed")
subtitle(sprintf("r = %.3f, p = %.3f", r, p))
In the first vs second half of gain gambles correlation (r = 0.579, p = 0.030), we reject the null hypothesis, as p < .05.
In the first vs second half of loss gambles correlation (r = 0.396, p = 0.161), we fail to reject the null hypothesis, as p > .05.
In the first vs second half of mixed gambles correlation (r = 0.698, 0.006), we reject the null hypothesis, as p < .05.
The correlations are relatively strong in the gain and mixed correlations, indicating higher risk-related stability in those conditions. That is, the temporal point of the game has little bearing on the level of risk-taking when considering gains or mixed gambles. It appears as though there is less stability in the loss condition, implying perhaps that as the game progresses, a player may become desperate –– at a certain point, they may be left with a small amount of points and thus might be less sensitive to magnitude and instead more focused on avoiding certain loss. 

for i = 1 : NUM_PARTICIPANTS

    individual_data = alldata(i).data;

    safe_trials = individual_data(individual_data(:,7) == 0,:);
    gain_safe = safe_trials(safe_trials(:,4) > 0 & safe_trials(:,5) == 0, 9);
    loss_safe = safe_trials(safe_trials(:,4) == 0 & safe_trials(:,5) < 0, 9);
    mixed_safe = safe_trials(safe_trials(:,4) > 0 & safe_trials(:,5) < 0, 9); 

    risky_trials = individual_data(individual_data(:,7) == 1,:);
    gain_risky = risky_trials(risky_trials(:,4) > 0 & risky_trials(:,5) == 0, 9);
    loss_risky = risky_trials(risky_trials(:,4) == 0 & risky_trials(:,5) < 0, 9);
    mixed_risky = risky_trials(risky_trials(:,4) > 0 & risky_trials(:,5) < 0, 9);

    alldata(i).gain_safe = median(gain_safe);
    alldata(i).loss_safe = median(loss_safe);
    alldata(i).mixed_safe = median(mixed_safe);

    alldata(i).gain_risky = median(gain_risky);
    alldata(i).loss_risky = median(loss_risky);
    alldata(i).mixed_risky = median(mixed_risky);

end

gain_safe_mean = mean([alldata.gain_safe]);
gain_safe_sem = get_sem([alldata.gain_safe]);

loss_safe_mean = mean([alldata.loss_safe]);
loss_safe_sem = get_sem([alldata.loss_safe]);

mixed_safe_mean = mean([alldata.mixed_safe]);
mixed_safe_sem = get_sem([alldata.mixed_safe]);

gain_risky_mean = mean([alldata.gain_risky]);
gain_risky_sem = get_sem([alldata.gain_risky]);

loss_risky_mean = mean([alldata.loss_risky]);
loss_risky_sem = get_sem([alldata.loss_risky]);

mixed_risky_mean = mean([alldata.mixed_risky]);
mixed_risky_sem = get_sem([alldata.mixed_risky]);

mean_rt = [gain_safe_mean, gain_risky_mean, mixed_safe_mean, mixed_risky_mean, loss_safe_mean, loss_risky_mean];
sem_rt = [gain_safe_sem, gain_risky_sem, mixed_safe_sem, mixed_risky_sem, loss_safe_sem, loss_risky_sem];

figure
barorder = {'gain safe','gain risky','mixed safe','mixed risky','loss safe','loss risky'};
bar(mean_rt);
set(gca,'xtick',1:6,'xticklabel',barorder);
hold on
er = errorbar(mean_rt, sem_rt);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
xlabel('Choice type');
ylabel('Response time (s)');
title('Mean of median response times');

p_gain = signrank([alldata.gain_safe], [alldata.gain_risky])
p_mixed = signrank([alldata.mixed_safe], [alldata.mixed_risky])
p_loss = signrank([alldata.loss_safe], [alldata.loss_risky])
Gain: p = 0.0580 (fail to reject null hypothesis, a = 0.05)
Loss: p = 0.0012 (reject null hypothesis, a = 0.05)
Mixed: p = 0.1040 (fail to reject null hypothesis, a = 0.05)
The primary feature to remark about this model is the outlier –– mean of median response times for the "loss risky" condition. This indicates significant hesitation when subjects are given the option of certain loss or a gamble between no loss and a potentially greater loss. This makes sense, as the subject is trying to figure out how to best avoid loss, while neither of these options are desirable.

tgaincutoff = [1.5 1.75 1.9 2.1 2.4   2.6 3 3.2 4 4.5 5.1]; %riskygain/safe ratio
tmixedcutoff = [0.4 0.6 0.8 0.9 1.1   1.3 1.6 2 2.8 3.5 5.1]; %-riskygain/riskyloss ratio
tlosscutoff = tgaincutoff; %riskyloss/safe ratio

NUM_COLS = size(alldata(1).data,2);
NUM_BINS = length(tgaincutoff) - 1;
gain_bin_matrix = zeros(NUM_PARTICIPANTS, NUM_BINS);
loss_bin_matrix = zeros(NUM_PARTICIPANTS, NUM_BINS);
mixed_bin_matrix = zeros(NUM_PARTICIPANTS, NUM_BINS);

for i = 1 : NUM_PARTICIPANTS
    
    individual_data = alldata(i).data;

    col_exact = zeros(NUM_TRIALS,1);
    col_binned = zeros(NUM_TRIALS,1);

    gain_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) == 0);
    gain_ratios = individual_data(gain_trials,4) ./ individual_data(gain_trials,3);
    col_exact(gain_trials) = gain_ratios;
    
    gain_bins = arrayfun(@(r) get_bin(r, tgaincutoff),gain_ratios);
    col_binned(gain_trials) = gain_bins;
 
    loss_trials = individual_data(individual_data(:,4) == 0 & individual_data(:,5) < 0);
    loss_ratios = individual_data(loss_trials,5) ./ individual_data(loss_trials,3);
    col_exact(loss_trials) = loss_ratios;

    loss_bins = arrayfun(@(r) get_bin(r, tlosscutoff),loss_ratios);
    col_binned(loss_trials) = loss_bins;

    mixed_trials = individual_data(individual_data(:,4) > 0 & individual_data(:,5) < 0);
    mixed_ratios = abs(individual_data(mixed_trials,4) ./ individual_data(mixed_trials,5));
    col_exact(mixed_trials) = mixed_ratios;

    mixed_bins = arrayfun(@(r) get_bin(r, tmixedcutoff),mixed_ratios);
    col_binned(mixed_trials) = mixed_bins;

    alldata(i).data(:,NUM_COLS+1) = col_exact;
    alldata(i).data(:,NUM_COLS+2) = col_binned;

    for b = 1 : NUM_BINS

        gain_bin = tgaincutoff(b);
        updated_gain_trials = alldata(i).data(gain_trials,:);
        binned_trials = updated_gain_trials(updated_gain_trials(:,NUM_COLS + 2) == gain_bin, :);
        percent_gain_gambles = sum(binned_trials(:,7) == 1) / size(binned_trials,1);
        gain_bin_matrix(i,b) = percent_gain_gambles;

        loss_bin = tlosscutoff(b);
        updated_loss_trials = alldata(i).data(loss_trials,:);
        binned_trials = updated_loss_trials(updated_loss_trials(:,NUM_COLS + 2) == loss_bin, :);
        percent_loss_gambles = sum(binned_trials(:,7) == 1) / size(binned_trials,1);
        loss_bin_matrix(i,b) = percent_loss_gambles;
        
        mixed_bin = tmixedcutoff(b);
        updated_mixed_trials = alldata(i).data(mixed_trials,:);
        binned_trials = updated_mixed_trials(updated_mixed_trials(:,NUM_COLS + 2) == mixed_bin, :);
        percent_mixed_gambles = sum(binned_trials(:,7) == 1) / size(binned_trials,1);
        mixed_bin_matrix(i,b) = percent_mixed_gambles;
        
    end

end

% gains
gain_means = mean(gain_bin_matrix);
gain_sem = get_sem(gain_bin_matrix);

loss_means = mean(loss_bin_matrix);
loss_sem = get_sem(loss_bin_matrix);

mixed_means = mean(mixed_bin_matrix);
mixed_sem = get_sem(mixed_bin_matrix);

figure
subplot(1,3,1)
plot(gain_means)
hold on
er = errorbar(gain_means, gain_sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
set(gca,'xtick',1:NUM_BINS,'xticklabel',{})
xlabel("Gamble value decile")
ylabel("% gain gambles chosen")
title("Gain")

subplot(1,3,2)
plot(loss_means)
hold on
er = errorbar(loss_means, loss_sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
set(gca,'xtick',1:NUM_BINS,'xticklabel',{})
xlabel("Gamble value decile")
ylabel("% loss gambles chosen")
title("Loss")

subplot(1,3,3)
plot(mixed_means)
hold on
er = errorbar(mixed_means, mixed_sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
set(gca,'xtick',1:NUM_BINS,'xticklabel',{})
xlabel("Gamble value decile")
ylabel("% mixed gambles chosen")
title("Mixed")

function b = get_bin(num,bins)
    for i = 2:length(bins)
        if num < bins(i)
            b = bins(i - 1);
            return
        end
    end
    b = -1;
end
function s = get_sem(M)
    s = std(M) / (sqrt(length(M)));
end
