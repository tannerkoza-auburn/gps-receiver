function varargout = dimHandle(varargin)

    switch nargin

        case 2
            numMeas = length(varargin{1});

            sz1 = size(varargin{1});
            sz2 = size(varargin{2});

            if sz1(1) ~= numMeas

                varargout{1} = varargin{1}';

            else

                varargout{1} = varargin{1};

            end

            if sz2(2) ~= numMeas

                varargout{2} = varargin{2}';

            else

                varargout{2} = varargin{2};

            end

        case 4

            numMeas = length(varargin{1});

            sz1 = size(varargin{1});
            sz2 = size(varargin{2});
            sz3 = size(varargin{3});
            sz4 = size(varargin{4});

            if sz1(1) ~= numMeas

                varargout{1} = varargin{1}';

            else

                varargout{1} = varargin{1};

            end

            if sz2(1) ~= numMeas

                varargout{2} = varargin{2}';

            else

                varargout{2} = varargin{2};

            end

            if sz3(2) ~= numMeas

                varargout{3} = varargin{3}';

            else

                varargout{3} = varargin{3};

            end

            if sz4(2) ~= numMeas

                varargout{4} = varargin{4}';

            else

                varargout{4} = varargin{4};

            end

    end   

end