function theta=findmax_new(areas,Pm,N,resolve1)
zz=[];
Pm=Pm(:);
d_Pm=diff(Pm);
d_Pm=[-1;d_Pm;1];

dpos=find(d_Pm>=0);
jj=1;
for ii=1:length(dpos)-1
    if d_Pm(dpos(ii)+1)<=0;
      zz(jj)= dpos(ii);
      jj=jj+1;
    end
end

theta=zeros(N,1);

if length(zz)<N
    theta=zeros(N,1);
elseif  length(zz)==N
    theta(1:length(zz))=areas(zz); 
else
    [xd,I]=sort(Pm(zz), 'descend');    
    temp_alpha=areas(zz(I));

    theta=[];
    theta(1)=temp_alpha(1);
    ii=2;
    while length(theta)~=N
        if min(abs(theta- temp_alpha(ii)))> resolve1;
            theta=[theta; temp_alpha(ii)];
        end
        ii=ii+1;
    end

end

[theta,I]=sort(theta);   