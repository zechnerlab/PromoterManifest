function [ MSN2_trace, time ] = MSN2_trace_generator_oneVal( inhibitor_conc, time_input, pulse_parameters )
%MSN2_TRACE_GENERATOR: This function generates MSN2 nuclear localisation
%time traces
%   nuclear localisation: input the 1-NM-PP1 concentration
%   time_input: input the time_input vector
%   pulse_parameters: input the following
%                       number of pulses
%                       pulse durations
%                       pulses intervals
%                       e.g [4 5 5] for 4 5min pulses with 5 min intervals.

%original 0.05;
time = 0:0.2:time_input;
MSN2_trace = zeros(1, length(time));
c1 = 0.60;
c2 = 0.60;

if inhibitor_conc == 100
    nuclear_localisation = 313.2;
    c1 = 1.11;
    c2 = 0.97;
elseif inhibitor_conc == 175
    nuclear_localisation = 535;
elseif inhibitor_conc == 275
    nuclear_localisation = 744.5;
    c1 = 0.61;
    c2 = 0.81;
elseif inhibitor_conc == 413
    nuclear_localisation = 932.3;
elseif inhibitor_conc == 690
    nuclear_localisation = 1107.8;
    c1 = 0.59;
    c2 = 0.57;
elseif inhibitor_conc == 1117
    nuclear_localisation = 1259;
elseif inhibitor_conc == 3000
    nuclear_localisation = 1410.1;
    c1 = 1.07;
    c2 = 0.29;
else
    nuclear_localisation = inhibitor_conc;
    if nuclear_localisation > 1100 
        nuc_int = 1110:10:1400;
        c1_int = 0.59:(1.07-0.59)/length(nuc_int):1.07;
        c2_int = 0.57:-(.57-.287)/length(nuc_int):0.287;
        %Find the closest c1 and c2 values:
        [index index] = min(abs(nuclear_localisation-nuc_int));
        c1 = c1_int(index);
        c2 = c2_int(index);
    elseif nuclear_localisation > 1400
        c1 = 1.07;
        c2 = 0.287;
    end
end

%SINGLE PULSE
if length(pulse_parameters) == 2
    for i=1:length(time)
        if time(i) > 0 && time(i) <= pulse_parameters(2)
            MSN2_trace(i) = nuclear_localisation * (1-exp(-c1*time(i))); %ON
            temp_nuclear_localisation = MSN2_trace(i);
        elseif time(i) > pulse_parameters(2)
            MSN2_trace(i) = temp_nuclear_localisation * exp(-c2*(time(i)-pulse_parameters(2))); %OFF
        end
    end
%SEVERAL PULSES
elseif length(pulse_parameters) == 3
    temp_going_out = 0;
    for k=1:pulse_parameters(1)
        for i=1:length(time)
            if time(i) >= (k-1)*pulse_parameters(2) + (k-1)*pulse_parameters(3) &&  time(i) < k*pulse_parameters(2) + (k-1)*pulse_parameters(3)
                MSN2_trace(i) = nuclear_localisation * (1-exp(-c1*(time(i)-((k-1)*pulse_parameters(2) + (k-1)*pulse_parameters(3))))); %ON
                if MSN2_trace(i) < temp_going_out
                    MSN2_trace(i) = temp_going_out;
                end
                temp_nuclear_localisation = MSN2_trace(i);
            else
                if k == pulse_parameters(1) && time(i) >= k*pulse_parameters(2) + (k-1)*pulse_parameters(3)
                    MSN2_trace(i) = temp_nuclear_localisation * exp(-c2*(time(i)-(k*pulse_parameters(2) + (k-1)*pulse_parameters(3)))); %OFF
                    temp_going_out = MSN2_trace(i);
                elseif time(i) >= k*pulse_parameters(2) + (k-1)*pulse_parameters(3) && time(i) < k*pulse_parameters(2) + k*pulse_parameters(3)
                    MSN2_trace(i) = temp_nuclear_localisation * exp(-c2*(time(i)-(k*pulse_parameters(2) + (k-1)*pulse_parameters(3)))); %OFF
                    temp_going_out = MSN2_trace(i);
                end
            end
        end
    end
end

MSN2_trace = MSN2_trace(end);
end

