%% Basic computation of temperature using 15000 points of data
function [Tj,parameters] = determine_Tj(R_on_meas,Tj_k,Ron_k)

    % Least Square Fitting through Pseudo-inverse Method
    Ron_k=Ron_k'*1000;
    Tj_k = Tj_k';
    X=[ones(length(Tj_k),1),Tj_k,Tj_k.^2];
    F=inv(X'*X);
    parameters=F*X'*Ron_k;
    % determine Junction Temperature based on given Rds_on
    % Tj =parameters(1)+R_on_meas*parameters(2)+R_on_meas^2*parameters(3)
    root_part=sqrt(parameters(2)^2-4*parameters(3)*(parameters(1)-R_on_meas)); 
    Tj=(-parameters(2)+root_part)/(2*parameters(3)) ;
end
        


