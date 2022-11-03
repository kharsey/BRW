function [Phi_tilde] = shrink_Phi(Phi_tilde, Phi_hat, ev_restr) 
% shrink matrix Phi.tilde to matrix Phi.hat
% until largest eigenvalue of Phi.tilde
% is less than ev.restr
% (e.g. Kilian's stationarity adjustment)
  maxeig = max(abs(eig(Phi_tilde)));
  Phi_diff0 = Phi_hat - Phi_tilde;
  delta = 1;
  while (maxeig>ev_restr)
    delta = delta - .01;
    Phi_diff = delta*Phi_diff0;
    Phi_tilde = Phi_hat - Phi_diff;
    maxeig = max(abs(eig(Phi_tilde)));
  end

end