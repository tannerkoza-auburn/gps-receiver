function varargout = dimHandle(varargin)
% DESCRIPTION: dimHandle transposes input vectors and matrices into the
% correct dimensions to be used in the solution methods, if needed.
%
% AUTHOR: Tanner Koza

  if nargin ~= nargout
      error('Number of Outputs needs to equal Inputs.')
  end

  numMeas = length(varargin{1});

  for i = 1:nargout

      sz = size(varargin{i});

      if isvector(varargin{i})

         if sz(1) ~= numMeas

             varargout{i} = varargin{i}'; %#ok<*AGROW> 

         else

             varargout{i} = varargin{i};

         end

      else

        if sz(2) ~= numMeas

             varargout{i} = varargin{i}';

        else

             varargout{i} = varargin{i};

        end

      end

  end

end