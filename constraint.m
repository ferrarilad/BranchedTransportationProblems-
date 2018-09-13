function [beq, diffMeasure,a,b] = constraint(n,x,example,dsmt,smooth,M)
    diffMeasure = zeros(n,n);
    switch example
        case 1                                                              % transport from dirac to a segment
            a = [2;.1]; b = [0;.5];                                         % I want to see the phase transition between concentrated and diffuse
            diffMeasure(round(3*n/4), round(n/4):round(3*n/4)) = 1;
            diffMeasure = diffMeasure * n^2 / sum(abs(diffMeasure(:)));     % normalize to probability measure
            diffMeasure(round(n/4), round(n/2)) = -n^2;
        
        case 2                                                              % transport from dirac to a segment
            a = [2;.1]; b = [0;5];                                          % same cost as above but adapted
            diffMeasure(round(3*n/4), round(n/4):round(3*n/4)) = 1;
            diffMeasure = 10*diffMeasure * n^2 / sum(abs(diffMeasure(:)));  % each measure has mass 10
            diffMeasure(round(n/4), round(n/2)) = -10*n^2;
            
        case 3                                                              % transport from dirac to a segment
            a = [2;.1]; b = [0;50];                                         % same cost as above but adapted
            diffMeasure(round(3*n/4), round(n/4):round(3*n/4)) = 1;
            diffMeasure = 100*diffMeasure * n^2 / sum(abs(diffMeasure(:))); % each measure has mass 100
            diffMeasure(round(n/4), round(n/2)) = -100*n^2;
        
        case 4                                                              % four points on the vertices of a rectangle
            a = [10;.05];         b = [0;1];                                % minimize total lenght of the support of sigma 
            diffMeasure(round(n/5),round(n/5)) = 3;
            diffMeasure(round(n/5),round(4*n/5)) = -1;
            diffMeasure(round(3*n/5),round(n/5)) = -1;
            diffMeasure(round(3*n/5),round(4*n/5)) = -1;
            diffMeasure = diffMeasure * n^2 / sum(abs(diffMeasure(:))); % normalize to probability measure

        case 5                                                              % four points on the vertices of a rectangle rotated to show the anisotropy of the method
            a = [10;.05];   b = [0;1];                                      % minimize total lenght of the support of sigma 
            alpha = 2*pi/360*17;
            diffMeasure(round(n*(.5+.3*cos(alpha))),round(n*(.5+.3*sin(alpha)))) = 3;
            diffMeasure(round(n*(.5+.3*cos(alpha+pi/2))),round(n*(.5+.3*sin(alpha+pi/2)))) = -1;
            diffMeasure(round(n*(.5+.3*cos(alpha+pi))),round(n*(.5+.3*sin(alpha+pi)))) = -1;
            diffMeasure(round(n*(.5+.3*cos(alpha+3*pi/2))),round(n*(.5+.3*sin(alpha+3*pi/2)))) = -1;
            diffMeasure = diffMeasure * n^2 / sum(abs(diffMeasure(:)));     % normalize to probability measure
            
        case 6                                                              % 4 points on the vertices of a square one source 3 sinks
            a = [10;.05];    b = [0;1];                                     % minimize total lenght of the support of sigma 
            theta = linspace(0,2*pi,5);
            for j = 1 : numel(theta)-1
                if (j == 1)
                    diffMeasure(round(n*(.5+.42*cos(theta(j)+pi/7))),round(n*(.5+.42*sin(theta(j)+pi/7)))) = n^2*(-1);
                else
                    diffMeasure(round(n*(.5+.42*cos(theta(j)+pi/7))),round(n*(.5+.42*sin(theta(j)+pi/7)))) = n^2;
                end
            end
            
       case 7                                                               % 4 points on the vertices of a square 2 sources 2 sinks
            a = [10;.05];     b = [0;1];                                    % minimize total lenght of the support of sigma 
            theta = linspace(0,2*pi,5);
            for j = 1 : numel(theta)-1
                diffMeasure(round(n*(.5+.42*cos(theta(j)+pi/7))),round(n*(.5+.42*sin(theta(j)+pi/7)))) = n^2*(-1)^j;
            end
           
        case 8                                                              % 3 points on the vertices of an equilateral triangle one source 2 sinks
            a = [10;.05];    b = [0;1];                                     % minimize total lenght of the support of sigma 
            theta = linspace(0,2*pi,4);
            for j = 1 : numel(theta)-1
                diffMeasure(round(n*(.5+.3*cos(theta(j)+pi/7))),round(n*(.5+.3*sin(theta(j)+pi/7)))) = n^2*(-1)^(j);
            end
            
        case 9                                                              % 4 points on the vertices of a square 2 sources 2 sinks
            a = [3;1;.1]; b= [0;2/3;38/40];                                 % three phases
            diffMeasure(round(4*n/15),round(2*n/15)) = 1;
            diffMeasure(round(6*n/15),round(2*n/15)) = 1;
            diffMeasure(round(8*n/15),round(2*n/15)) = 1;
            diffMeasure(round(10*n/15),round(2*n/15)) = 1;
            diffMeasure(round(4*n/15),round(13*n/15)) = -1;
            diffMeasure(round(6*n/15),round(13*n/15)) = -1;
            diffMeasure(round(8*n/15),round(13*n/15)) = -1;
            diffMeasure(round(10*n/15),round(13*n/15)) = -1;
            diffMeasure = diffMeasure * n^2 / sum(abs(diffMeasure(:)));     % normalize to probability measure          
            
         case 10                                                            % D points on a circle and a dirac
            a = [10;.05];    b = [0;1];                                     % minimize total lenght of the support of sigma 
            D = 8;
            theta = linspace(0,2*pi,D);
            for j = 1 : numel(theta)
                diffMeasure(round(n*(.5+.3*cos(theta(j)+pi/7))),round(n*(.5+.3*sin(theta(j)+pi/7)))) = 1;
            end
            diffMeasure = diffMeasure * n^2 / sum(abs(diffMeasure(:)));     % normalize to probability measure            
            diffMeasure(round(.5*n),round(.5*n)) = -n^2;   
    end 
    
    
    if smooth                                                               % Smoothing of the constraint measure
        SmoothKernel = exp(-(x-.5).^2/(2*(dsmt)^2));
        SmoothKernel = fftshift(SmoothKernel) / sum(SmoothKernel);
        diffMeasure = real(ifft2(fft2(diffMeasure).*fft2(SmoothKernel'*SmoothKernel)));
    end
    
    beq = M*diffMeasure(:);                                                 % Project diffMeasure on the space of piecewise P1 polynomials 
end