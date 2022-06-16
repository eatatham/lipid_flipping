%% This code determines the coefficients of the fluorescence decay equations described in Falzone & Accardi, 2020

%read in current metadata and clean and organise for maximum usability
in_out = readtable('metadatafilepath');
file_names = string(in_out.FileName);
in_out.avg_finish = round(in_out.avg_finish, 1);
in_out.fast_exp_start = round(in_out.fast_exp_start, 1);
in_out.fast_exp_finish = round(in_out.fast_exp_finish, 1);
in_out.slow_exp_start = round(in_out.slow_exp_start, 1);
in_out.PL_exp_start = round(in_out.PL_exp_start, 1);

%% First, fit protein-free fluorescence decay to single exponential function to find LiPF at t=0 and gamma as described in Falzone & Accardi, 2020
%This code modifies the method of Falzone & Accardi to account for an
%unknown slow reaction seen in both the protein free and flipping assay
%data sets

empty_directory = "nameofdirectory";
empty_file = "nameoffile";
empty_data = readtable("filepath" + empty_directory + "/" + empty_file + ".xlsx");
match_empty_file = find(file_names == empty_file);

%make arrays of relevant data
F_empty = table2array(empty_data(:,2)); %5 Hz fluorescence
t_raw_empty = table2array(empty_data(:,1));   %5 Hz time

%remove any data past 600s
long = find(t_raw_empty > 600);
F_empty(long) = [];
t_raw_empty(long) = [];

%need to choose which 5s to average over for normalisation separately for each data set. Plot the data and make decision from there.
% figure(1); plot(t_raw_empty, F_empty, '.')
%%% insert chosen start and end %%%
avg_finish_empty = in_out.avg_finish(match_empty_file);
avg_start_empty = avg_finish_empty - 5;
avg_begin_empty = find(t_raw_empty == avg_start_empty);
avg_end_empty = find(t_raw_empty == avg_finish_empty);

F_max_empty = mean(F_empty(avg_begin_empty:avg_end_empty)); %max fluorescence
F_PF = F_empty./F_max_empty;                %5 Hz normalised fluorescence

%plot ln of 5 Hz fluorescence and decide range over which line is straight
%will assume this section is dominated by the fast reaction of diothionite
%reduction
% figure(2); plot(t_raw_empty, log(F_PF), '.')
%%% insert chosen start and end %%%
fast_exp_start = in_out.fast_exp_start(match_empty_file);
fast_exp_finish = in_out.fast_exp_finish(match_empty_file);
fast_exp_begin = find(t_raw_empty == fast_exp_start);
fast_exp_end = find(t_raw_empty == fast_exp_finish);

%pick out the straight line data from the larger data set
fast_fit_y = F_PF(fast_exp_begin:fast_exp_end);
fast_fit_x = (0:0.2:length(fast_fit_y)/5-0.2)'; %create time vector for the range of the fast exponential, where the beginning of the exponential is t=0

%Fit data in the time range of the straight section to single exponential from Falzone & Accardi, 2020
fast_fit_options = fitoptions('Method','NonlinearLeastSquares', 'Lower',[0 -Inf], 'Upper',[1 Inf], 'StartPoint',[0.5 0.1]);
fast_fit_type = fittype('Li_PF0 + (1 - Li_PF0).*exp(-gamma.*fast_fit_x)', 'dependent',{'fast_fit_y'}, 'independent',{'fast_fit_x'}, 'options',fast_fit_options);
[fast_fit, fast_fit_stats] = fit(fast_fit_x, fast_fit_y, fast_fit_type);

figure(3)
plot(fast_fit, fast_fit_x, fast_fit_y)
ylabel('F_{PF}')
xlabel('time, s')
title('Fast Reaction Single Exponential Fit (protein free)')

%save values of Li_PF0 and gamma as variables
results_fast_fit = coeffvalues(fast_fit);
Li_PF0 = results_fast_fit(1);
gamma = results_fast_fit(2);
fast_fit_ci = confint(fast_fit, 0.95);
%and update table of data
in_out.Li_PF0(match_empty_file) = Li_PF0;
in_out.gamma(match_empty_file) = gamma;
in_out.Li_PF0_ci95low(match_empty_file) = fast_fit_ci(1,1);
in_out.Li_PF0_ci95high(match_empty_file) = fast_fit_ci(2,1);
in_out.Li_PF0_rmse(match_empty_file) = fast_fit_stats.rmse;
in_out.gamma_ci95low(match_empty_file) = fast_fit_ci(1,2);
in_out.gamma_ci95high(match_empty_file) = fast_fit_ci(2,2);
in_out.gamma_rmse(match_empty_file) = fast_fit_stats.rmse;

%Extrapolate the fit to the rest of the time course of the protein free experiment
t1 = (0:0.2:600-fast_exp_start)';                             %create new 5Hz time vector so the exponential can start at t=0
t1 = round(t1, 1);
fast_empty = Li_PF0 + (1 - Li_PF0).*exp(-gamma.*t1);

%Subtract fast reaction from raw data, leaving only the unknown slow
%reaction
F_PF_exp = F_PF(fast_exp_begin:length(F_PF));    %pick out the normalised raw data where the reactions are happening (i.e. removing the flat part at the start)
slow_empty = F_PF_exp - fast_empty;

%Plot slow reaction data (looks like a double exponential) and adjust data ready for a generic double exponential curve fit
% figure(4); plot(t1, slow_empty, '.')
slow_exp_start = in_out.slow_exp_start(match_empty_file);     %%% insert chosen start %%%
slow_exp_begin = find(t1 == slow_exp_start);
slow_fit_y = slow_empty(slow_exp_begin:length(slow_empty));
slow_fit_x = t1(slow_exp_begin:length(slow_empty));
%slow_fit_x = (0:0.2:length(slow_fit_y)/5-0.2)'; %create time vector where the beginning of the exponential is t=0

%fit a generic double exponential
[slow_fit, slow_fit_stats] = fit(slow_fit_x, slow_fit_y, 'exp2');

figure(5)
plot(slow_fit, slow_fit_x, slow_fit_y)
ylabel('F_{PF}')
xlabel('time, s')
title('Slow Reaction Generic Double Exponential Fit (protein free)')

%save coefficients as variables
results_slow_fit = coeffvalues(slow_fit);
c1 = results_slow_fit(1);
c2 = results_slow_fit(2);
c3 = results_slow_fit(3);
c4 = results_slow_fit(4);
%and update table of data
in_out.c1(match_empty_file) = c1;
in_out.c2(match_empty_file) = c2;
in_out.c3(match_empty_file) = c3;
in_out.c4(match_empty_file) = c4;

%% Fit the triple exponential function of Falzone & Accardi (2020) to the proteoliposome data to find foward and backward flipping rates
%This code modifies the method of Falzone & Accardi to account for an
%unknown slow reaction seen in both the protein free and flipping assay
%data sets

PL_directory = "nameofdirectory";
PL_file = "nameoffile";
PL_data = readtable("filepath" + PL_directory + "/" + PL_file + ".xlsx");
match_PL_file = find(file_names == PL_file);

%make arrays of relevant data
F_PL_raw = table2array(PL_data(:,2)); %5 Hz fluorescence
t_raw_PL = table2array(PL_data(:,1));   %5 Hz time

%remove any data past 600s
long = find(t_raw_PL > 600);
F_PL_raw(long) = [];
t_raw_PL(long) = [];

%need to choose which 5s to average over for normalisation separately for each data set. Plot the data and make decision from there.
% figure(6); plot(t_raw_PL, F_PL_raw, '.')
%%% insert chosen end %%%
avg_finish_PL = in_out.avg_finish(match_PL_file);
avg_start_PL = avg_finish_PL - 5;
avg_begin_PL = find(t_raw_PL == avg_start_PL);
avg_end_PL = find(t_raw_PL == avg_finish_PL);

F_max_PL = mean(F_PL_raw(avg_begin_PL:avg_end_PL)); %max fluorescence
F_PL_norm = F_PL_raw./F_max_PL;                %5 Hz normalised fluorescence

%plot normalised data and determine where the decay starts
% figure(7); plot(t_raw_PL, F_PL_norm, '.')
%%% insert chosen start %%%
PL_exp_start = in_out.PL_exp_start(match_PL_file);
PL_exp_begin = find(t_raw_PL == PL_exp_start);
F_PL = F_PL_norm(PL_exp_begin:length(F_PL_norm));
t2 = (0:0.2:length(F_PL)/5-0.2)';               %create time vector where the beginning of the exponential is t=0
t2 = round(t2, 1);

%subtract slow reaction from PL data
PL_slow = c1.*exp(c2.*t2) + c3.*exp(c4.*t2);
F_tot = F_PL - PL_slow;

%Fit exponential
PL_fit_method = fitoptions('Method','NonlinearLeastSquares', 'Lower',[1E-5 1E-5 0], 'Upper',[1 100 1], 'StartPoint',[0.1 0.1 0.5]);
PL_fit_type = fittype('f0.*(Li_PF0 + (1-Li_PF0).*exp(-gamma.*t2)) + (1-f0)./((((-(a+B+gamma-sqrt((a+B+gamma).^2-4.*a.*gamma))/2)+a).*((-(a+B+gamma+sqrt((a+B+gamma).^2-4.*a.*gamma))/2)+B+gamma) - a.*B).*(a+B)).*(a.*((-(a+B+gamma+sqrt((a+B+gamma).^2-4.*a.*gamma))/2)+gamma).*((-(a+B+gamma-sqrt((a+B+gamma).^2-4.*a.*gamma))/2)+a+B).*exp((-(a+B+gamma-sqrt((a+B+gamma).^2-4.*a.*gamma))/2).*t2)+(-(a+B+gamma-sqrt((a+B+gamma).^2-4.*a.*gamma))/2).*B.*((-(a+B+gamma+sqrt((a+B+gamma).^2-4.*a.*gamma))/2)+a+B+gamma).*exp((-(a+B+gamma+sqrt((a+B+gamma).^2-4.*a.*gamma))/2).*t2))', 'dependent',{'F_tot'}, 'independent',{'t2'}, 'problem', {'Li_PF0' 'gamma'}, 'options',PL_fit_method);
[PL_exp, PL_stats] = fit(t2, F_tot, PL_fit_type, 'problem', {Li_PF0; gamma});

%save values of coefficients as variables
results_PL = coeffvalues(PL_exp);
B = results_PL(1);
a = results_PL(2);
f0 = results_PL(3);
PL_fit_ci = confint(PL_exp, 0.95);
%and update table of data
in_out.B(match_PL_file) = B;
in_out.a(match_PL_file) = a;
in_out.f0(match_PL_file) = f0;
in_out.B_ci95low(match_PL_file) = PL_fit_ci(1,1);
in_out.B_ci95high(match_PL_file) = PL_fit_ci(2,1);
in_out.B_rmse(match_PL_file) = PL_stats.rmse;
in_out.a_ci95low(match_PL_file) = PL_fit_ci(1,2);
in_out.a_ci95high(match_PL_file) = PL_fit_ci(2,2);
in_out.a_rmse(match_PL_file) = PL_stats.rmse;
in_out.f0_ci95low(match_PL_file) = PL_fit_ci(1,3);
in_out.f0_ci95high(match_PL_file) = PL_fit_ci(2,3);
in_out.f0_rmse(match_PL_file) = PL_stats.rmse;

figure(8)
plot(PL_exp, t2, F_tot)
ylabel('F_{tot}')
title('"Complex" Triple exponential fit; Scramblase')

%%
writetable(in_out, 'filepathandnameofoutput');