clc;
clear;
close all;
%New FDE function test script

%% Data Collection Toggles:
write_to_excel = 0;
error_analysis = 0;
generate_plots = 0;

start_time = cputime;
%% Input Variables
y0 = [0,0];
t0 = 0;
tfinal = 1;
h_largest = 0.1;
h_list = [h_largest, h_largest/2, h_largest/4, h_largest/8, h_largest/16, h_largest/32];
Num_Sims = length(h_list);
alpha =1.25;
%**************************************************************************
% EQUATION BANK
%**************************************************************************
% % Equation 1:
% fdefun = @(t,y) (40320/gamma(9-alpha))*t^(8-alpha) - 3*(gamma(5+alpha/2)/gamma(5-alpha/2))*t^(4-alpha/2) + 9/4*gamma(alpha+1) + (3/2*t^(alpha/2)-t^4)^3-y^(3/2);
% % Solution to Eq. 1:
%Eq1_Sol = t^8 - 3*t^(4+alpha/2)+9/4*t^alpha;
%***********
%Equation 2:
% fdefun = @(t,y) gamma(4+alpha)/6*t^3+t^(3+alpha) - y;
%Solution to Eq. 2:
%Eq2_Sol = t^(3+alpha);
%***********
% %Equation 3:
fdefun = @(t,y) gamma(5+alpha)/24*t^4+t^(8+2*alpha) - y^2
%Solution to Eq. 2:
%Eq3_Sol = t^(4+alpha);

%**************************************************************************
for AA = 1:Num_Sims
    h = h_list(AA)
    
    %% RUN MF_PCL
    [T,Y_MF_PCL] = MF_PCL(fdefun,t0,tfinal,y0,alpha, h);
    
    %% Get True Solution
    i = 1;
    for t = t0:h:tfinal
        true_solution(i) = t^(4+alpha); %Change depending on Eq being used!
        i = i +1;
    end
    
    %% RUN ABM
    [t,Y_ABM] = fde12(alpha,fdefun,t0,tfinal,y0, h);
    
    %% Generate Plots
    if generate_plots == 1
        time = t0:h:tfinal;
        figure(AA)
        plot(time, true_solution);
        hold on
        plot(time,Y_MF_PCL);
        hold on
        plot(time,Y_ABM)
        xlabel('time')
        ylabel('y')
        legend('Exact Solution', 'MF-PCL', 'ABM')
        savename = strcat('Eq3_alpha=',num2str(alpha),'_N=',num2str(1/h),'.png');
        saveas(gcf, savename)
    end

     %% Run Error Analysis
    if error_analysis == 1
        E2_sum_MF_PCL = 0;
        E2_sum_ABM = 0;
        for tt = 1:length(T)
            Error_ABM(tt) = true_solution(tt) - Y_ABM(tt);
            Error_MF_PCL(tt) = true_solution(tt) - Y_MF_PCL(tt);
            E2_ABM(tt) = (Error_ABM(tt))^2;
            E2_MF_PCL(tt) = (Error_MF_PCL(tt))^2;
            E2_sum_ABM = E2_sum_ABM + E2_ABM(tt);
            E2_sum_MF_PCL = E2_sum_MF_PCL + E2_MF_PCL(tt);
            figure(2)
            plot(Error_ABM)
            hold on 
            plot(Error_MF_PCL)
            legend('ABM', 'MF_PCL')
        end
        EL2_ABM = (h*E2_sum_ABM)^(0.5);
        EL2_MF_PCL = (h*E2_sum_MF_PCL)^(0.5);
        EPT_ABM = abs(Error_ABM(length(T)));
        EPT_MF_PCL = abs(Error_MF_PCL(length(T)));
    end

    %% Write to Excel Spreadsheet:
    if write_to_excel == 1 && error_analysis == 1
        export_info = {'alpha'      , alpha;
                       'h'          , h; 
                       'EL2_MF_PCL' , EL2_MF_PCL;
                       'EL2_ABM'    , EL2_ABM;
                       'EPT_MF_PCL' , EPT_MF_PCL;
                       'EPT_ABM'    , EPT_ABM}
        BB = AA*7 - 6;
        filename = 'ErrorAnalysis.xlsx';
        writetocell = strcat('G',num2str(BB));
        xlswrite(filename, export_info, 'Equation_3',  writetocell);  
    end

end


