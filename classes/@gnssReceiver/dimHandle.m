function varargout = dimHandle(varargin)

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