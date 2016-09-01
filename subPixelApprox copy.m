function [disp] = subPixelApprox(x2,Y1,Y2,Y3)
x1=x2-1;
x3=x2+1;
disp=3200;
  %Pixels must be different to compute the minimum   
  if ((Y1 == Y2) && (Y2 == Y3))
    return
  end
  % Three valid pixels are needed to compute the minimum
  if ((Y1 == 3200) || (Y3 == 3200))

    % If we don't have enough valid SADs to fit a parabola, we still
    % keep the integer disparity if its SAD is less than a valid
    % neighbor (which should always be the case) -- LJE
    if (Y1 ~= 3200)
      if (Y2 <= Y1)
          disp=x2;
          return
      else
          return 
      end
    elseif (Y3 ~= 3200)
      if (Y2 <= Y3) 
          disp=x2;
      else
          return 
      end
    end
  end
  
  
  % Y1 and Y3 should be larger than or equal to Y2, if not something
  % is wrong.
  if ((Y1 < Y2) || (Y3 < Y2))
    return 
  end
  
  
  % Invert the matrix and multiply it by X
  p1 = (x1*x1);
  p2 = (x2*x2);
  p3 = (x3*x3);

  q1 = x1;
  q2 =  x2;
  q3 =  x3;

  y1 = Y1;
  y2 =  Y2;
  y3 =  Y3;

  detM = p1*q2 + q1*p3 + p2*q3 - p3*q2 - q3*p1 - p2*q1;
  
  
  epsilon = 0.00001;
  if ((-epsilon < detM) && (detM < epsilon))
    return
  end

  b = (p3 - p2)*y1 + (p1 - p3)*y2 + (p2 - p1)*y3;
  a = (q2 - q3)*y1 + (q3 - q1)*y2 + (q1 - q2)*y3;


  if (epsilon > (a / detM))
    return 
  end
  
  %The min is where the deriv of y = ax^2+bx+c is 0, so xmin = -0.5*b/a
  xsubpix = -0.5 * (b / a);


  if (xsubpix < q1)
    disp=x1;
    return
  end
  if (xsubpix > q3)
    disp= x3;
    return
  end
  
  %return the subpixel value 
  disp=xsubpix;

end