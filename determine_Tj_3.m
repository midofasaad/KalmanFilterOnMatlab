%% The lookup tables computed and encapsulated as a six order polynomial to estimate temperature
function Tj_3 = determine_Tj_3(a_l,b_l,c_l,current,R_on_meas)

     I=[1; current; current^2; current^3;current^4;current.^5;current.^6]
     a=  a_l*I*1000 
     b=  b_l*I*1000
     c=  c_l*I*1000
     root_part=sqrt(b^2-4*c*(a-R_on_meas)); 
     Tj_3=(-b+root_part)/(2*c) ;
        
end
