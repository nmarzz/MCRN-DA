dim = 2;
a = ones(dim);
a(2,1)=2;
a(1,2)=3;
a(2,2)=4;
a
b = reshape(a,dim*dim,1)
c = reshape(b,dim,dim)

i=4
[j,k]=vectomat(i,dim,dim)
i=mattovec(j,k,dim,dim)
