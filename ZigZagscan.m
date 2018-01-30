function Vect=ZigZagscan(X,dct_co)
% ZigZagscan Transform an matrix to a vector using Zig Zag Scan.
[~, N]=size(X);
Vect=zeros(N,N);
Vect(1)=X(1,1);
s=0;
for k=1:2*N-1
    if s<dct_co
    if k<=N
        if mod(k,2)==0
        j=k;
        for i=1:k
        Vect(i,j)=X(i,j);
        j=j-1;s=s+1;    
        end
        else
        i=k;
        for j=1:k   
        Vect(i,j)=X(i,j);
        i=i-1;s=s+1;
        end
        end
    else
        if mod(k,2)==0
        p=mod(k,N); j=N;
        for i=p+1:N
        Vect(i,j)=X(i,j);
        j=j-1;s=s+1;    
        end
        else
        p=mod(k,N);i=N;
        for j=p+1:N   
        Vect(i,j)=X(i,j);
        i=i-1;s=s+1; 
        end
        end
    end
    end
end
 
        
        
