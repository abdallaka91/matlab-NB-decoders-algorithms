function H=readalist(fname)

[fid, msg] = fopen(fname,'r','native');
temp = fscanf(fid, '%d %d\n', 2);
N = temp(1);
M = temp(2);
temp = fscanf(fid, '%d %d\n', 2);
dv = temp(1);
dc = temp(2);
[symDeg,errmsg] = fscanf(fid,"%d",N);
[chkDeg,errmsg] = fscanf(fid,"%d",M);
[symIdx,errmsg] = fscanf(fid,"%d",[dv,N]);



%% N
%% M
%% dv
%% dc
%% symDeg
%% chkDeg
%% symIdx

H=zeros(M,N);
for idx=1:N
    for jdx=1:dv
	if (symIdx(jdx,idx)>0)
	   H(symIdx(jdx,idx),idx)=1;
	end
    end
end

