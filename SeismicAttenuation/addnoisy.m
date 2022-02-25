function[noisydata]=addnoisy(data,k,P);

for i=1:(length(data(1,:)))
    if P==1;
    noisy = awgn(data(:,i),k,'measured');
    noisydata(:,i)=noisy;
    end
    if P==2;
        noisy = addNoise(data(:,i), k);
        noisydata(:,i)=noisy;
    end
end 