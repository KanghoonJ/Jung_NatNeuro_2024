function Data_cdf = gen_cdf(varargin)

Data = varargin{1};
if(length(varargin)>1)
    xrange = varargin{2};
    Data = Data(~isnan(Data));
    Data_cdf.x = [xrange(1); [sort(Data,'ascend')]; xrange(2)];
    Data_cdf.y = [0; [[1:numel(Data)]'./numel(Data)]; 1];
else
    Data = Data(~isnan(Data));
    Data_cdf.x = [sort(Data,'ascend')];
    Data_cdf.y = [1:numel(Data)]'./numel(Data);

end
