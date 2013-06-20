function waveTrain = CreateWaveletDouble(iPoint, nPoint)
%CreateWaveletDouble. Create a wavelet used for wavelet transforms.
%   waveTrain = CreateWaveletDouble(iPoint, nPoint)
%   Creates a wavelet of total length nPoint, using iPoint functional taps
%   in the wavelet train.  Wavelet is based on the form [-0.5 1 -0.5],
%   which operates as a double derivative of some data trace.  If nPoint is
%   even, iPoint must also be even; conversely is nPoint is odd. If nPoint
%   is greater than iPoint, the arms of the wavelet are zero-padded.  This
%   is useful as an input into a fast wavelet transform.
%   
%   This wavelet is divided into 4 sections.  In the first section [0 1], the
%   wavelet is constructed using the slope -1.  In the second section [1
%   2], the wavelet is calculated using y = 3.*x - 4.  The wavelet is
%   symmetrical about 2, so that section 3 has slope -3 and section 4 has
%   slope 4.  Positive portion of wavelet can be calculated as
%   (iPoint+1)./3; and I have called this the wavelet width in other
%   functions.
%
%   EXAMPLE:
%   w = CreateWaveletDouble(3,3)
%   Produces:
%    w = [-0.5 1 -0.5]
%
%   w = CreateWaveletDouble(3,5)
%   Produces:
%    w = [0 -0.5 1 -0.5 0]
%
%   w = CreateWaveletDouble(8,8)
%   Produces:
%    w = [-2/9 -4/9 0 6/9 6/9 0 -4/9 -4/9]
%
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software is licenced under the Creative Commons Attribution
%   (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)

%%
if iPoint > nPoint
  error('CreateWaveletDouble:errWaveletLength', ['Total wavelet length is shorter than requested ' ...
    'wavelet function length: aborting.']);
end

if mod(nPoint,2)==0
  if mod(iPoint,2)~=0
    iPoint = iPoint - 1;
    fprintf(1, '%s\n', 'Total wavelet length is even, while requested functional');
    fprintf(1, '%s\n', ['wavelet is odd: new functional wavelet length' ...
      num2str(iPoint) '.']);
  end
else
  if mod(iPoint,2)==0
    iPoint = iPoint - 1;
    fprintf(1, '%s\n', 'Total wavelet length is odd, while requested functional');
    fprintf(1, '%s\n', ['wavelet is even: new functional wavelet length' ...
      num2str(iPoint) '.']);
  end
end

% t = 4./(iPoint + 1).*(1:iPoint);
% waveTrain = -t(t<=1);
% waveTrain = [waveTrain 3.*t((t>1 & t<=2))-4];
t = 2./((iPoint + 1)/2).*(1:(iPoint/2));
waveTrain = interp1(0:4, [0 -1./2 1 -1./2 0], t);
waveTrain = waveTrain./sqrt(sum(waveTrain.^2))./2;



if mod(nPoint,2) == 0
  waveTrain = [zeros(1,(nPoint./2) - numel(waveTrain)), waveTrain];
else
  waveTrain = [zeros(1,((nPoint+1)./2) - numel(waveTrain)), waveTrain];
end  
if mod(iPoint, 2) == 1
  waveTrain = [waveTrain(1:end-1), fliplr(waveTrain)];
else
  waveTrain = [waveTrain, fliplr(waveTrain)];  
end


% t = linspace(-nPoint./2, nPoint./2, nPoint);
% waveTrain =  2./(sqrt(3.*sqrt(iPoint)).*pi^(1./4)).*(1-t.^2./iPoint).*exp(-t.^2./2./iPoint);

