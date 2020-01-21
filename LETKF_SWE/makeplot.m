function makeplot(nsol,xvals,truth,xbar,analysis,choice,yobs,k,h,RMSE)
       %  Plot the truth, analysis mean, and observations on graph #1
       subplot(2,2,1);
       plot(truth, 'k-');
       hold on
          %axis([1 40 -10 15]);
          plot(xvals,xbar, 'b--');
          plot(choice, yobs, 'r*');
          title(['Time = ',num2str(k*h),',  RMSE = ',num2str(RMSE)])
       hold off

       %  Plot the ensemble and observations on graph #2
       subplot(2,2,2); 
       plot(xvals,analysis(:,1), 'b-');
       hold on
          %axis([1 40 -10 15]);
          for j = 2:nsol
             plot(xvals,analysis(:,j), 'b-');
          end
%%          plot(choice, yobs, 'r*');
          title(['Time = ',num2str(k*h)])
       hold off
