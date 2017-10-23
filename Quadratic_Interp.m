function [Q_low, Q_mid, Q_up] = Quadratic_Interp(t_low,t_mid,t_up)
 
 syms tau
 Q_low = ((tau - t_mid)/(t_low - t_mid)) * ((tau - t_up)/(t_low - t_up));

 Q_mid = ((tau - t_low)/(t_mid - t_low)) * ((tau - t_up)/(t_mid - t_up));
        
 Q_up = ((tau - t_low)/(t_up - t_low)) * ((tau - t_mid)/(t_up - t_mid));

end