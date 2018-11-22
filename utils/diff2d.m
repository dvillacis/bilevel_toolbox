function [mat, matx, maty]=diff2d(dim, method)
    if length(method>2)
        bdry=method(2);
        method=method(1);
    else
        bdry='d';
    end

    xdim=dim(1);
    ydim=dim(2);

    if bdry~='n' && bdry~='d'
        error('Invalid boundary conditions');
    end
    
    mmx=@(pts) mx(xdim, ydim, pts, bdry);
    mmy=@(pts) my(xdim, ydim, pts, bdry);

    if method=='f'
        M0=m0(xdim, ydim);
        matx=mmx(1)-M0;
        maty=mmy(1)-M0;
    elseif method=='b'
        M0=m0(xdim, ydim);
        matx=M0-mmx(-1);
        maty=M0-mmy(-1);
    elseif method=='c'
        matx=(mmx(1)-mmx(-1))/2;
        maty=(mmy(1)-mmy(-1))/2;
    elseif method=='5'
        %f'(x) ≈ [-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h)]/12h
        matx=(-mmx(2)+8*mmx(1)-8*mmx(-1)+mmx(-2))/12;
        maty=(-mmy(2)+8*mmy(1)-8*mmy(-1)+mmy(-2))/12;
    else
        error('Invalid finite differences method')
    end

    mat=[matx; maty];
end

% Base matrices

% This is just the identity, for the use of the point x in the differences formula
function [M0, vdim]=m0(xdim, ydim)
    vdim=xdim*ydim;
    M0=speye(vdim, vdim);
end

% Generate linear index for all entries (besides fault) of a xdim*ydim matrix
function i0=mki0(xdim, ydim, fault)
    ix0=repmat((1:xdim)', 1, ydim);
    iy0=repmat(1:ydim, xdim, 1);
    ix0(fault)=[];
    iy0(fault)=[];
    i0=sub2ind([xdim, ydim], ix0(:), iy0(:));
end

% This generates matrices for using neighbouring points at distance pts in the x direction
function M=mx(xdim, ydim, pts, bdry)
    vdim=xdim*ydim;
    % Generate indices for the points we want
    ix=repmat(pts+(1:xdim)', 1, ydim);
    iy=repmat((1:ydim), xdim, 1);
    % Boundary conditions
    if bdry=='n'
        % Neumann boundary: reset if index over bounds, reset to closest inside domain
        % in order to extend values.
        ix(ix>xdim)=xdim;
        ix(ix<1)=1;
        fault=[];
    else
        % Dirichlet boundary conditions: remove index, in order to not use point,
        % which will set it to zero.
        fault=(ix>xdim | ix<1);
        ix(fault)=[];
        iy(fault)=[];
    end
    % Generate linear indices of all source and target image coordinates
    i=sub2ind([xdim ydim], ix(:), iy(:));
    i0=mki0(xdim, ydim, fault);
    % Map source coordinates into target coordinates
    M=sparse(i0, i, 1, vdim, vdim);
end

% This generates matrices for using neighbouring points at distance pts in the y direction
function M=my(xdim, ydim, pts, bdry)
    vdim=xdim*ydim;
    % Generate indices for the points we want
    iy=repmat(pts+(1:ydim), xdim, 1);
    ix=repmat((1:xdim)', 1, ydim);
    % Boundary conditions
    if bdry=='n'
        % Neumann boundary: reset if index over bounds, reset to closest inside domain
        % in order to extend values.
        iy(iy>ydim)=ydim;
        iy(iy<1)=1;
        fault=[];
    else
        % Dirichlet boundary conditions: remove index, in order to not use point,
        % which will set it to zero.
        fault=(iy>ydim | iy<1);
        ix(fault)=[];
        iy(fault)=[];
    end
    % Generate linear indices of all source and target image coordinates
    i=sub2ind([xdim ydim], ix(:), iy(:));
    i0=mki0(xdim, ydim, fault);
    % Map source coordinates into target coordinates
    M=sparse(i0, i, 1, vdim, vdim);
end
