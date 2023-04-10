function [R2,RMSEa, RMSEr] = ErrorCalculator(Yapprox,Yvoulu)
    RMSEa = sqrt(mean((Yapprox-Yvoulu).^2));
    RMSEr = sqrt(mean(((Yapprox-Yvoulu)./Yvoulu).^2));
    Ymoy = mean(Yvoulu);
    R2 = 1-(sum((Yapprox-Yvoulu).^2)/sum((Yvoulu-Ymoy).^2));
end