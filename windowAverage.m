function [out] = windowAverage(signal,size)

%size must be odd
side = (size-1) / 2;

for i = 1:length(signal)
    if i <= (side)
        out(i) = mean(signal(1:i));
    elseif i > (length(signal) - side)
        out(i) = mean(signal(i:end));
    else
        out(i) = mean(signal((i-side):(i+side)));
    end
end