function P2 = extrapolate_profile(P1,extrapolate_thresh,extrapolate_mode)

% Extrapolate profile to whole region

extrapolate_all_points = false;

extrapolate_mode_fn2 = @(X,Y) [ones(size(X)), X, Y];
extrapolate_mode_fn3 = @(X,Y) [ones(size(X)), X, Y, X.^2, Y.^2, X.*Y];
% extrapolate_mode_fn4 = @(X,Y) [ones(size(X)), X, Y, 1./(X-2), 1./(Y-2), 1./(X+2), 1./(Y+2)];
extrapolate_mode_fn4 = @(X,Y) [ones(size(X)), X, Y, exp(-1./(1-(4.5^-2)*((X-2).^2+(Y-2).^2))),...
                                                    exp(-1./(1-(4.5^-2)*((X-2).^2+(Y+2).^2))),...
                                                    exp(-1./(1-(4.5^-2)*((X+2).^2+(Y-2).^2))),...
                                                    exp(-1./(1-(4.5^-2)*((X+2).^2+(Y+2).^2)))];

if extrapolate_mode > 0
    
    P1v = P1;
    P1v(abs(P1) <= extrapolate_thresh*max(max(abs(P1)))) = 0;
    [X,Y,P1v] = find(P1v); % Find training data for model

    X = 2*X/size(P1,1)-1;
    Y = 2*Y/size(P1,1)-1;

    if extrapolate_mode == 1
        R = sum(P1v)/size(P1v,1);

        P2 = zeros(size(P1))+R;
        if not(extrapolate_all_points)
            P2(abs(P1) > extrapolate_thresh*max(max(abs(P1)))) = P1(abs(P1) > extrapolate_thresh*max(max(abs(P1))));
        end
        
    elseif extrapolate_mode > 1
        if extrapolate_mode == 2
            A = extrapolate_mode_fn2(X,Y);
        elseif extrapolate_mode == 3
            A = extrapolate_mode_fn3(X,Y);
        elseif extrapolate_mode == 4
            A = extrapolate_mode_fn4(X,Y);
        elseif extrapolate_mode == 5
            A = extrapolate_mode_fn5(X,Y);
        end
    
        R = A\P1v;
        
        if extrapolate_all_points
            P2 = ones(size(P1));
            [X,Y] = find(P2); % Find extrapolation points 

            X = 2*X/size(P1,1)-1;
            Y = 2*Y/size(P1,1)-1;

            if extrapolate_mode == 2
                Rv = extrapolate_mode_fn2(X,Y)*R;
            elseif extrapolate_mode == 3
                Rv = extrapolate_mode_fn3(X,Y)*R;
            elseif extrapolate_mode == 4
                Rv = extrapolate_mode_fn4(X,Y)*R;
            elseif extrapolate_mode == 5
                Rv = extrapolate_mode_fn5(X,Y)*R;
            end            
        
            P2(sub2ind(size(P2),round((X+1)*size(P1,1)/2),round((Y+1)*size(P1,1)/2))) = Rv;
        else
            P2 = zeros(size(P1));
            P2(abs(P1) <= extrapolate_thresh*max(max(abs(P1)))) = 1;
            [X,Y] = find(P2); % Find extrapolation points 

            X = 2*X/size(P1,1)-1;
            Y = 2*Y/size(P1,1)-1;

            if extrapolate_mode == 2
                Rv = extrapolate_mode_fn2(X,Y)*R;
            elseif extrapolate_mode == 3
                Rv = extrapolate_mode_fn3(X,Y)*R;
            elseif extrapolate_mode == 4
                Rv = extrapolate_mode_fn4(X,Y)*R;
            elseif extrapolate_mode == 5
                Rv = extrapolate_mode_fn5(X,Y)*R;
            end            
        
            P2(sub2ind(size(P2),round((X+1)*size(P1,1)/2),round((Y+1)*size(P1,1)/2))) = Rv;

            P2(abs(P1) > extrapolate_thresh*max(max(abs(P1)))) = P1(abs(P1) > extrapolate_thresh*max(max(abs(P1))));
        end
        
    end
else
    P2 = P1;
end

