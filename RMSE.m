function [estimate_RMSE] = RMSE(estimate,original_Signal,start,finish)
% Calculate the RMSE 
estimate_RMSE= [];
for i=start:finish
    y=estimate(start:i);
    error=(y - original_Signal(start:i)).^2 ;   % Errors
    mean_error= mean(error);   % Mean Squared Error
    RMSE = sqrt(mean_error);  % Root Mean Squared Error
    estimate_RMSE=[estimate_RMSE,RMSE];

end