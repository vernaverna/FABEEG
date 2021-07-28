

addpath('/m/cs/work/eleppaah/MNE/share/matlab/') % Add MNE path
loadp = '/m/cs/work/eleppaah/siblingData/dir1/'; %The path containing the MEG data

splitData = false; %Split the data into smaller parts
splitTimes = [1:10 15 20 30 60];
if ~splitData; splitTimes = 1; end

subj = dir(loadp); %Subjects in the directory
for i=3:size(subj,1) %First two . and .. %
    filename = [loadp subj(i).name '/spont_ecoh_trans_sss.fif'];
    s = filename((findstr(filename,'/nro')+4):(findstr(filename,'/spont')-1)); %Subject ID
    s = str2num(s);
    
    savefile = fullfile([datapath 'dtxt/s' num2str(s) '.txt']);
    if splitData
        savefile = fullfile([datapath 'dtxt/s' num2str(s) 'split0.1.60.txt']);
    end
    
    fprintf(filename)
    fprintf('\n')
    fprintf(savefile)
    fprintf('\n')
    
    data = fiff_setup_read_raw(filename);
    [X,times] = fiff_read_raw_segment(data,data.first_samp,data.last_samp);
    fprintf(['Initial ' num2str(size(X,2)) ' samples.\n'])
    D = size(X,1);
    N = size(X,2);
    
    nfft=2048;
    t1 = 0; t2 = 180; %Starting and ending times
    
    sfr = data.info.sfreq;
    time1 = [t1*sfr t2*sfr];
    
    X = X(:,time1);
    Nch = 306;
    if splitData; Nch = Nch*2; end
    
    for splitT=splitTimes
        %fprintf(['splitT=' num2str(splitT) '\n']);
        if splitT==splitTimes(1)
            zs = find(X(1,:)==0);
            
            X(:,zs) = [];
            Y = X;
        end
        
        if splitData
            t = size(X,2);
            t = min(splitT*sfr, size(Y,2)/2);
            X = [Y(1:306,1:t); Y(1:306,(end-t+1):end)];
        end
        
        freq = [1 3];
        for f=1:13; freq = [freq; freq(f,2) freq(f,2)+2+f/5]; end
        freq = [freq; 52 54.8];
        for f=15:20; freq = [freq; freq(f,2) freq(f,2)+2+f/5]; end
        
        spectra = zeros(0,size(freq,1));
        
        t = 0;
        N = size(X,2);
        
        while true
            fprintf('.');
            t = max(t)+(1:N);
            if max(t) > size(X,2)
                break;
            else
                %Include zeros, processing in R
                if true %size(intersect(t,zs),2) == 0 || splitData
                    tmp = zeros(Nch,size(freq,1));
                    if welch; tmp = zeros(Nch,nfft/2+1); end
                    for m=1:Nch
                        if welch
                            [power,f]=pwelch(X(m,t),[],[],nfft,data.info.sfreq);
                            tmp(m,:)=power;
                        else
                            for f=1:size(freq,1)
                                tmp(m,f) = bandpower(X(m,t), sfr, freq(f,:));
                            end
                        end
                    end
                    spectra = [spectra; tmp];
                end
            end
        end
        fprintf('\n');
        X = spectra*1e24;
    end
    
    fprintf('%d samples, dim=%d per subject.\n',size(X,1),size(X,2));
    
    dlmwrite(fullfile([datapath 'f.txt']), freq);
    
    if ~splitData
        savefile = fullfile([datapath 'dtxt/s' num2str(s) ...
            '.' num2str(N) '.txt']);
    else
        savefile = fullfile([datapath 'dtxt/s' num2str(s) ...
            'split' num2str(splitT) '.' num2str(N) '.txt']);
    end
    
    dlmwrite(savefile, X);
    fprintf(['Data written in file: ' savefile '\n']);
end


