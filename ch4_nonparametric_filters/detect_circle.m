% ==============================================================

%returns arclength of unit circle intersecting cell

% ==============================================================

function b = detect_circle( x, y, c )

    if (x+c/2 < sqrt(2)/2 && y+c/2 < sqrt(2)/2 )||(x-c/2 > 1  || y-c/2 > 1 )
        b = 0;
    else if (x^2+y^2+1/2*c^2+c*(x+y) >1 && x^2+y^2+1/2*c^2-c*(x+y) <1)
            [x1,y1] =linecirc(inf, x-c/2, 0, 0, 1);
            [x2,y2] =linecirc(inf, x+c/2, 0, 0, 1);
            [x3,y3] =linecirc(0, y-c/2, 0, 0, 1);
            [x4,y4] =linecirc(0, y+c/2, 0, 0, 1);
            %strangely linecirc may return NaN in columns.
            A = [x1(1,:) x2(1,:) x3(1,:) x4(1,:); y1(1,:) y2(1,:) y3(1,:) y4(1,:)];
            B = A(1,:);
            C = A(2,:);
            ind = find(abs(B-x)<=c/2+1e-15 & ~isnan(B) );
            B = B(ind);C = C(ind);
            ind = find(abs(C-y)<=c/2+1e-15 & ~isnan(C) );
            B = B(ind);C = C(ind);
            b = acos([B(1) C(1)]*[B(2) C(2)]');
        else
            b = 0;
        end
    end
  
end

