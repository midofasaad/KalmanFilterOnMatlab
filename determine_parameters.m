function parameters = determine_parameters(x,y)
    % Least Square Fitting through Pseudo-inverse Method
    x=x';
    y =y';
    X=[ones(length(x),1),x,x.^2,x.^3,x.^4,x.^5,x.^6];
    F=inv(X'*X);
    parameters=F*X'*y;
    parameters=parameters'
end