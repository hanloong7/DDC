function f = nLL(dim, b, beta, theta30, theta31, theta32, theta33, X, D)


LL = 0;
[T N] = size(X); 


  EV2 = ev(dim, b, beta(1), beta(2), theta30, theta31, theta32, theta33);
  
  s = (0 : 1 : dim)';
  
  % probability of replacement
  Prob2 = exp( -beta(1)*ones(dim,1) +... 
   b*EV2(:,2) )./(exp( -0.001*beta(2)*s(1:dim)+ b*EV2(:,1))+exp(-beta(1)*ones(dim,1)+b*EV2(:,2))); 

[ dim]
size(Prob2)
mean(X')
pause;
      
  % LL calculation
   for k=1:N
      % for t=1
      if D(1,k) == 1
        if X(1,k) == 0
         LL = log(Prob2(X(1,k)+1)) + log(theta30) + LL;
        elseif X(1,k) == 1
         LL = log(Prob2(X(1,k)+1)) + log(theta31) + LL;
        elseif X(1,k) == 2
         LL = log(Prob2(X(1,k)+1)) + log(theta32) + LL;
        else
         LL = log(Prob2(X(1,k)+1)) + log(theta33) + LL;
        end
         
      elseif D(1,k) == 0
        if X(1,k) == 0
         LL = log(1-Prob2(X(1,k)+1)) + log(theta30) + LL;
        elseif X(1,k) == 1
         LL = log(1-Prob2(X(1,k)+1)) + log(theta31) + LL;
        elseif X(1,k) == 2
         LL = log(1-Prob2(X(1,k)+1)) + log(theta32) + LL;
        else
         LL = log(1-Prob2(X(1,k)+1)) + log(theta33) + LL;
        end
      end
      
      % for t
      for i = 2:T
          
       if D(i,k) == 1
           if D(i-1,k) ==1
             LL = log(Prob2(X(i,k)+1)) + LL;
           elseif D(i-1,k) ==0
            if X(i,k) ==  X(i-1,k)
             LL = log(Prob2(X(i,k)+1))+ log(theta30) + LL;
            elseif X(i,k) ==  X(i-1,k) + 1;
             LL = log(Prob2(X(i,k)+1))+ log(theta31) + LL;
            elseif X(i,k) ==  X(i-1,k) + 2;
             LL = log(Prob2(X(i,k)+1))+ log(theta32) + LL;
            else
             LL = log(Prob2(X(i,k)+1))+ log(theta33) + LL;
            end  
           end
           
       elseif D(i,k) == 0
          if D(i-1,k) == 1
            LL = log(1-Prob2(X(i,k)+1)) + LL;
          elseif D(i-1,k) == 0
           if X(i,k) ==  X(i-1,k)
            LL = log(1-Prob2(X(i,k)+1)) + log(theta30) + LL;
           elseif X(i,k) ==  X(i-1,k) + 1;
            LL = log(1-Prob2(X(i,k)+1)) + log(theta31) + LL;
           elseif X(i,k) ==  X(i-1,k) + 2;
            LL = log(1-Prob2(X(i,k)+1)) + log(theta32) + LL;
           else
            LL = log(1-Prob2(X(i,k)+1)) + log(theta33) + LL;
           end
          end
       end
      end
   end

      f = -LL;
