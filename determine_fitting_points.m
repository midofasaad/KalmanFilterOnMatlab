% 1. Determine corresponding set of Junction Temperatures for a given drain current.
% 2. Determine corresponding set of Onstate Resistances for a given drain current. 
function [Ron_k, Tj_k] = determine_fitting_points(Idt,Tj,Rds_on,current)
Tj_k=[];
Ron_k=[];
    for i=1:length(Tj(:,1)) 
     Tj_k=[lin_intrplte(Idt,Tj(i,:),current), Tj_k]; 
 
     Ron_k=[lin_intrplte(Idt,Rds_on(i,:),current), Ron_k];
    end 
     
     end
     
% Interpolation of y_ref corresponding to x_ref, which isn't neceassirly
% included in x and y due to discretization.
function y_ref = lin_intrplte(x,y,x_ref)
    [minValue,i]=min(abs(x-x_ref));
     if x[i]<x_ref
         m=(y(i+1)-y(i))/(x(i+1)-x(i));
         y_ref=m*(x_ref-x(i))+y(i);
     end 
     if x[i]>x_ref 
         m=(y(i)-y(i-1))/(x(i)-x(i-1));
         y_ref=m*(x_ref-x(i-1))+y(i-1);
     
     else 
          y_ref=y(i);
     end 
         
end



     
