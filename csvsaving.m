function data=csvsaving(D)
    
L=length(D.Frequency);

if handles.plot_type==1;
    N=size(D.Power,2);
else
    N=size(D.Amplitude,2);
end

    data=cell(L+13,(N*2)+1);
    dstart=15;
    data{1,1}='Time-frequency analysis toolbox';
    data{2,1}=date;
    data{3,1}=[];
    data{4,1}='PARAMETERS';
    if handles.calc_type==1
        data{5,1}='Analysis type';
        data{5,2}='Wavelet';
        data{6,1}='Wavelet type';
        data{6,2}=D.Wavelet_type;
    else
        data{5,1}='Analysis type';
        data{5,2}='Windowed Fourier';
        data{6,1}='Window type';
        data{6,2}=D.Window_type;
    end    
    data{7,1}='Sampling frequency (Hz)';
    data{7,2}=D.Sampling_frequency;
    data{8,1}='Maximum frequency (Hz)';
    data{8,2}=D.fmax;
    data{9,1}='Minimum frequency (Hz)';
    data{9,2}=D.fmin;
    data{10,1}='Central frequency';
    data{10,2}=D.f0;
    data{11,1}='Preprocessing';
    data{11,2}=D.Preprocessing;
    
    data{12,1}='Cut Edges';
    data{12,2}=D.COI;
    data{13,1}='Time start (s)';
    data{13,2}=min(D.Time);
    data{14,1}='Time end (s)';
    data{14,2}=max(D.Time);


data{dstart,1}='Frequency';
for l=1:L;
data{l+dstart,1}=D.Frequency(l);
end

if handles.plot_type==1
for j=1:N
    data{dstart,j+1}=['Power ',num2str(j)];
       for k=1:L
        data{k+dstart,j+1}=D.Power(k,j);  
       end   
    
end

else
    
for j=1:N
    data{dstart,j+1}=['Amplitude ',num2str(j)];
      for k=1:L
        data{k+dstart,j+1}=D.Amplitude(k,j);  
      end   
    
end
end