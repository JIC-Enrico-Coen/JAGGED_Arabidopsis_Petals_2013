function [cs,fs,Tgmatrix] = tensorsToComponentsTmatrix( ts, preferredframes )
% [cs,fs] = tensorsToComponents( ts, preferredframes )
%   Convert a set of growth tensors ts given as an N*6 matrix into an N*3
%   matrix cs of the growth amounts and a 3*3*N matrix of the axis frames.
%   preferredframes is an optional argument.  If supplied, then the columns
%   of the matrices in fs will be permuted to lie as close as possible to
%   those of preferredframes, and cs will be permuted likewise.  Otherwise,
%   the values in cs will be in descending order.

    havePreferredFrames = nargin >= 2 && ~isempty( preferredframes );
    numtensors = size( ts, 1 );
    cs = zeros( numtensors, 3 );
    fs = zeros( 3, 3, numtensors );
    for i=1:numtensors
        [c,f,gmatrix] = tensorComponentsTmatrix( ts(i,:) );
        if havePreferredFrames
            perm = alignFrames( preferredframes(:,:,i), f );
            cs(i,:) = c(perm);
            fs(:,:,i) = f(:,perm);
        else
            cs(i,:) = c;
            fs(:,:,i) = f;
        end
        Tgmatrix{i} = gmatrix;
    end
end
