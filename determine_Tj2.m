function Tj_2  = determine_Tj2(a_l,b_l,c_l,Idt,current,R_on_meas)
     a=  lin_intrplte(Idt,a_l,current)*1000;
     b=  lin_intrplte(Idt,b_l,current)*1000;
     c=  lin_intrplte(Idt,c_l,current)*1000;
     root_part=sqrt(abs(b^2-4*c*(a-R_on_meas))); 
     Tj_est=(-b+root_part)/(2*c) ;
     end
     
% Interpolation of y_ref corresponding to x_ref, which isn't neceassirly
% included in x and y due to discretization.
function y_ref = lin_intrplte(x,y,x_ref)
    [minValue,i]=min(abs(x-x_ref));
     if x[i]<x_ref
         m=(y(i+1)-y(i))/(x(i+1)-x(i));
         y_ref=m*(x_ref-x(i))+y(i);
     end 
     if x[i]>x_ref && i>=2 
         m=(y(i)-y(i-1))/(x(i)-x(i-1));
         y_ref=m*(x_ref-x(i-1))+y(i-1);
     
     else 
          y_ref=y(i);
     end 
         
end
