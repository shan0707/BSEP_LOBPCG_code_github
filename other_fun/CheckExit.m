function Exit = CheckExit(y, i, criterion, grad_prev, threshold)
%% Exit when the curve is flat or rises abnormally
y_prev = max([y(i - 1), y(i - 2)]);
Exit = false;
switch criterion
    case 1
        if nargin < 4
          grad_prev = (y(6) - y(1))/5;
        end
        if nargin < 5
          threshold = 0.5;
        end
        grad_current = y(i) - y(i - 1);
      
    case 2
        %% threshold is adaptive
        if nargin < 4
          grad_prev = (y(6) - y(1))/5;
        end
        if nargin < 5
            grad_list = y(2 : 6) - y(1 : 5);
            threshold = 10*std(grad_list);
        end
        grad_current = y(i) - y(i - 1);

    case 3
        %% window size = 5
        if nargin < 4
            grad_prev = (y(6) - y(1))/5;
        end
        if nargin < 5
            threshold = 0.5;
        end
        window_size = 5;
        grad_current  = (y(i) - y(i - window_size))/window_size;
       
    case 4
        %% the grad_current of the current five steps compares with the previous ten steps
        threshold = 0.5;
        window_size1 = 10;
        window_size2 = 5;
        grad_prev = (y(i) - y(i - window_size1))/window_size1;
        grad_current = (y(i) - y(i - window_size2))/window_size2;
end

if grad_current > threshold*grad_prev || y(i) > y_prev
    if grad_current > threshold*grad_prev && y(i) <= y_prev
       disp(criterion);
    end
    Exit = true;
end