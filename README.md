# Solving-a-Multi-objective-hierarchical-hub-location-problem-with-k-and-M-M-1-queuing-framework

model:

sets:
i1/1..9/;
h/1..2/:s,lambda1; 
c/1..3/:lambda;
t/1..2/;
q/1..2/;
cql1/1..6/;
m1(i1,i1):w,c1,d,t1;
f(i1,c):ftild,ll,cql;
iz(i1,i1,i1):z,g,ff,zz,v,tt,y;
cct(c,c,c):ft,u;
cq(q,c):s1;
endsets

data:

w=   6     6     6     5     6     5     4     6     6
     6     4     6     5     4     5     5     6     4
     4     6     5     4     6     5     6     6     6
     6     6     4     6     4     6     5     5     5
     5     5     6     4     5     6     5     4     4
     4     6     6     4     5     4     4     4     4
     4     4     6     4     6     6     6     4     5
     5     5     6     4     6     5     4     6     5
     6     6     6     6     4     4     5     4     5;

c1=45.6764   54.8840   27.5331   29.0821   37.0542   53.3436   32.2507   35.5948   41.4103
   47.4445   27.3453   36.9935   51.0788   27.2790   39.7259   37.1174   49.6358   33.8896
   38.5162   38.2803   32.7961   42.3911   32.1975   39.6776   27.8936   25.4621   47.3408
   27.5146   28.1996   49.0021   41.4958   28.6996   35.1316   28.9592   26.2907   30.6687
   31.8693   53.8569   37.9424   29.3486   30.5172   52.0016   53.2615   30.0697   45.6033
   52.4001   25.1390   52.3194   50.5909   32.1986   36.0774   53.6840   44.4735   30.5053
   29.5713   48.2473   30.4554   43.6617   37.5180   28.3361   42.2563   46.9517   36.0545
   49.7745   49.5191   32.9141   35.5286   26.4896   48.4076   26.7934   44.4324   43.7686
   41.1503   51.0608   29.3662   40.3975   52.0815   36.6922   32.0434   38.5277   48.4068;

d=  2.5679    7.7234    5.8511    3.5815    8.8582    3.5522    2.2045    5.6480    8.1962
    8.5057    7.5638    6.3573    3.1950    5.0721    2.8219    8.5020    3.6212    8.3930
    7.4300    6.5102    6.1093    3.5937    2.7778    4.0767    7.1123    5.4223    7.5733
    5.4075    4.6503    3.4542    5.0499    3.8065    4.2314    5.4203    6.3684    2.6910
    5.0510    7.6811    4.1087    4.1777    4.8610    4.9692    6.0497    6.7539    3.8331
    5.1275    5.7298    5.2965    8.4637    6.1643    5.5550    3.6610    4.7686    4.3475
    4.1444    4.4551    3.6134    5.0115    3.8355    2.5986    5.2119    4.5721    6.7581
    5.5596    8.5730    7.9102    3.2937    6.2199    3.8374    8.7416    8.9159    2.9559
    5.5754    8.1316    3.3634    8.3342    6.9785    7.6071    5.8276    2.2642    7.0486;

t1= 41.0330   20.9780   45.6657   31.5386   30.3163   32.7578   40.9966   41.5508   27.9941
    39.9902   36.8360   39.3429   37.4896   37.5221   29.3816   39.1559   49.0595   24.6097
    36.1738   46.4560   31.2882   27.5542   23.2331   24.8445   21.0081   35.9400   28.4302
    40.9432   40.0753   25.7277   28.7132   47.1892   25.3630   22.0642   29.7544   33.2026
    39.9958   25.7130   32.8476   38.5127   46.3896   32.6866   29.5880   23.1689   35.8143
    25.3440   31.0675   34.4607   27.9584   44.5328   22.8269   35.9259   38.3288   33.7227
    23.8404   33.8218   23.6183   44.7313   27.8218   37.9557   39.6334   43.3641   46.2611
    49.9724   49.4491   37.6852   49.4799   37.8307   34.1277   32.2286   32.7036   35.5416
    25.1336   24.6921   26.7856   41.9075   20.6754   40.8785   44.5994   22.7247   48.3087;

Ftild=  1607  1698  2306  2213  2318  1954  1561   2484  1553
        2154  1531  2077  2001  2318  1932  1899   1667  2238
        1994  2244  1683  1971  2222  2325  2027   1606  1769;
cql=
     3     5     2     5     2     5     4     5     5
     5     2     2     3     4     5     4     5     5
     2     3     6     6     3     2     5     3     2;

alpha_h=0.6718;alpha_c=0.2891;
alphahat_h=0.2051;alphahat_c=0.5312;
ft=800;p=2;p0=3;mu=1200;kk=5;

enddata

!Objective Function;

min=(@sum(i1(i):@sum(i1(m):@sum(h(j):@sum(c(l):(w(i,m)+w(m,i))*c1(i,j)*d(i,j)*z(i,j,l)))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_h*c1(i,j)*d(i,j)*g(i,j,l))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_c*c1(i,j)*d(i,j)*ff(i,j,l))))+
@sum(c(i):@sum(q(j):Ftild(i,j)*ll(i,j)))+
@sum(c(i):@sum(c(k):@sum(t(r):ft(i,k,r)*u(i,k,r)))));
!min=beta;

obj1=(@sum(i1(i):@sum(i1(m):@sum(h(j):@sum(c(l):(w(i,m)+w(m,i))*c1(i,j)*d(i,j)*z(i,j,l)))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_h*c1(i,j)*d(i,j)*g(i,j,l))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_c*c1(i,j)*d(i,j)*ff(i,j,l))))+
@sum(c(i):@sum(q(j):Ftild(i,j)*ll(i,j)))+
@sum(c(i):@sum(c(k):@sum(t(r):ft(i,k,r)*u(i,k,r)))));


@for(i1(i):@sum(h(j):@sum(c(k):z(i,j,k)))=1); !eq3;


@for(i1(i):@for(h(j):@for(c(k):z(i,j,k)<=z(j,j,k)))); !eq4;


@for(h(j):@for(c(i):@sum(h(k):z(j,i,k))<=z(i,i,i))); !eq5;


@sum(h(i):@sum(c(k):z(i,i,k)))=p; !eq6;


@sum(c(i):z(i,i,i))=p0; !eq7;


@for(c(i):z(i,i,i)<=@sum(q(j):ll(i,j))); !eq8;


@for(c(i):@for(c(j):@sum(t(r):u(i,j,r))<=z(j,j,j))); !eq9;


@for(c(i):@for(c(j):@sum(t(r):u(i,j,r))<=z(i,i,i))); !eq10;


@for(c(i):@for(c(j):z(i,i,i)+z(j,j,j)-@sum(t(r):u(i,j,r))<=1)); !eq11;


@for(i1(i):@for(c(j):(@sum(c(k):ff(i,j,k))-@sum(c(k):ff(i,k,j)))=@sum(i1(o):w(i,o)*(@sum(h(o1):z(i,o1,o)-z(o,o1,j)))))); !eq12;


@for(i1(i):@for(h(j):@for(c(k):g(i,j,k)>=@sum(i1(l1):(w(i,l1)+w(l1,i))*(z(i,j,k)-z(l1,j,k)))))); !eq13;


@for(h(i):@for(c(j):zz(i,j,i)=0)); !eql4;


@for(i1(i):@for(i1(j):@for(c(o):@for(c(o1):@for(q(r):@sum(h(jj):(t1(i,jj)+s(jj)+alphahat_h*t1(jj,o)+s1(r,o))*z)+
  @sum(t(rr):alphahat_c*t1(rr,o)*v(rr,o,o1)+alphahat_c*tt(rr,o,o1))+s1(r,o)*z(o1,o1,o1)+
  @sum(h(n):(alphahat_h*t1(o1,n)+s(n)+t1(j,n))*z)+@sum(h(nn):s(nn)*y(j,nn,o1))<=beta))))); !eq15;


@for(c(k):@for(c(i):@for(t(r):v(r,i,k)<=u(r,i,k)))); !eq16;


@for(c(k):@for(c(i):@for(t(r):v(r,i,k)<=z(k,k,k)))); !eq17;


@for(c(k):@for(c(i):@for(t(r):u(r,i,k)+z(k,k,k)-v(r,i,k)<=1))); !eq18;


@for(i1(m):@for(h(n):@for(c(k):y(m,n,k)<z(m,n,k)))); !eq19;


@for(i1(m):@for(h(n):@for(c(k):y(m,n,k)<z(n,n,k)))); !eq20;


@for(i1(m):@for(h(n):@for(c(k):z(m,n,k)+z(n,n,k)-y(m,n,k)<=1))); !eq21;


@for(c(i):@for(q(r):s1(r,i)=@sum(q(j):(( (lambda(i)/mu)^cql(r,i) )*mu)/( (@exp(@lgm(cql(r,i)) ))*((cql(r,i)*mu-lambda(i))^2))*
  (((1+@sum(cql1(j):((lambda(i)/mu)^j)* (1/@exp(@lgm(j))) )+(((lambda(i)/mu)^cql(r,i))*(1/@exp(@lgm(cql(r,i)))))* (lambda(i)/(mu-lambda(i)))))^-1)*ll(r,i)))); !eq22;

@for(h(j):s(j)= (lambda1(j)/mu) /(1-lambda1(j)/mu)-(((1+kk)*((lambda1(j)/mu)^(kk+1)))/(1-(lambda1(j)/mu)^(kk+1)))); !eq23;


@for(c(i):lambda(i)=@sum(i1(j):@sum(i1(m):(w(j,m)+w(m,j))*z(j,i,i)))+ @sum(i1(k):@sum(h(ii):g(k,ii,i))) ); !eq24;


@for(h(jj):lambda(jj)=@sum(i1(k):@sum(c(j):g(k,jj,j)))); !eq25;



@for(i1(k):@for(i1(i):@for(i1(r):@bin(v(k,r,i)))));
@for(i1(k):@for(c(i):@bin(ll(k,i))));
@for(i1(k):@for(i1(i):@for(i1(r):@bin(u(k,r,i)))));
@for(i1(k):@for(c(i):@for(i1(r):@bin(z(k,r,i)))));
@for(i1(k):@for(i1(i):@for(i1(r):@bin(zz(k,r,i)))));
@for(i1(k):@for(i1(i):@for(i1(r):@bin(y(k,r,i)))));
end
model:

sets:
i1/1..9/;
h/1..2/:s,lambda1; 
c/1..3/:lambda;
t/1..2/;
q/1..2/;
cql1/1..6/;
m1(i1,i1):w,c1,d,t1;
f(i1,c):ftild,ll,cql;
iz(i1,i1,i1):z,g,ff,zz,v,tt,y;
cct(c,c,c):ft,u;
cq(q,c):s1;
endsets


!Objective Function;
min=(((-bestobjctv1- (@sum(i1(i):@sum(i1(m):@sum(h(j):@sum(c(l):(w(i,m)+w(m,i))*c1(i,j)*d(i,j)*z(i,j,l)))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_h*c1(i,j)*d(i,j)*g(i,j,l))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_c*c1(i,j)*d(i,j)*ff(i,j,l))))+
@sum(c(i):@sum(q(j):Ftild(i,j)*ll(i,j)))+
@sum(c(i):@sum(c(k):@sum(t(r):ft(i,k,r)*u(i,k,r)))))) /(-bestobjctv1))^2+((-bestobjctv2-beta)/(-bestobjctv2))^2)^(1/2);

obj1=(@sum(i1(i):@sum(i1(m):@sum(h(j):@sum(c(l):(w(i,m)+w(m,i))*c1(i,j)*d(i,j)*z(i,j,l)))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_h*c1(i,j)*d(i,j)*g(i,j,l))))+
@sum(i1(i):@sum(h(j):@sum(c(l):alpha_c*c1(i,j)*d(i,j)*ff(i,j,l))))+
@sum(c(i):@sum(q(j):Ftild(i,j)*ll(i,j)))+
@sum(c(i):@sum(c(k):@sum(t(r):ft(i,k,r)*u(i,k,r)))));


@for(i1(i):@sum(h(j):@sum(c(k):z(i,j,k)))=1); !eq3;


@for(i1(i):@for(h(j):@for(c(k):z(i,j,k)<=z(j,j,k)))); !eq4;


@for(h(j):@for(c(i):@sum(h(k):z(j,i,k))<=z(i,i,i))); !eq5;


@sum(h(i):@sum(c(k):z(i,i,k)))=p; !eq6;


@sum(c(i):z(i,i,i))=p0; !eq7;


@for(c(i):z(i,i,i)<=@sum(q(j):ll(i,j))); !eq8;


@for(c(i):@for(c(j):@sum(t(r):u(i,j,r))<=z(j,j,j))); !eq9;


@for(c(i):@for(c(j):@sum(t(r):u(i,j,r))<=z(i,i,i))); !eq10;


@for(c(i):@for(c(j):z(i,i,i)+z(j,j,j)-@sum(t(r):u(i,j,r))<=1)); !eq11;


@for(i1(i):@for(c(j):(@sum(c(k):ff(i,j,k))-@sum(c(k):ff(i,k,j)))=@sum(i1(o):w(i,o)*(@sum(h(o1):z(i,o1,o)-z(o,o1,j)))))); !eq12;


@for(i1(i):@for(h(j):@for(c(k):g(i,j,k)>=@sum(i1(l1):(w(i,l1)+w(l1,i))*(z(i,j,k)-z(l1,j,k)))))); !eq13;


@for(h(i):@for(c(j):zz(i,j,i)=0)); !eql4;


@for(i1(i):@for(i1(j):@for(c(o):@for(c(o1):@for(q(r):@sum(h(jj):(t1(i,jj)+s(jj)+alphahat_h*t1(jj,o)+s1(r,o))*z)+
  @sum(t(rr):alphahat_c*t1(rr,o)*v(rr,o,o1)+alphahat_c*tt(rr,o,o1))+s1(r,o)*z(o1,o1,o1)+
  @sum(h(n):(alphahat_h*t1(o1,n)+s(n)+t1(j,n))*z)+@sum(h(nn):s(nn)*y(j,nn,o1))<=beta))))); !eq15;


@for(c(k):@for(c(i):@for(t(r):v(r,i,k)<=u(r,i,k)))); !eq16;


@for(c(k):@for(c(i):@for(t(r):v(r,i,k)<=z(k,k,k)))); !eq17;


@for(c(k):@for(c(i):@for(t(r):u(r,i,k)+z(k,k,k)-v(r,i,k)<=1))); !eq18;


@for(i1(m):@for(h(n):@for(c(k):y(m,n,k)<z(m,n,k)))); !eq19;


@for(i1(m):@for(h(n):@for(c(k):y(m,n,k)<z(n,n,k)))); !eq20;


@for(i1(m):@for(h(n):@for(c(k):z(m,n,k)+z(n,n,k)-y(m,n,k)<=1))); !eq21;


@for(c(i):@for(q(r):s1(r,i)=@sum(q(j):(( (lambda(i)/mu)^cql(r,i) )*mu)/( (@exp(@lgm(cql(r,i)) ))*((cql(r,i)*mu-lambda(i))^2))*
  (((1+@sum(cql1(j):((lambda(i)/mu)^j)* (1/@exp(@lgm(j))) )+(((lambda(i)/mu)^cql(r,i))*(1/@exp(@lgm(cql(r,i)))))* (lambda(i)/(mu-lambda(i)))))^-1)*ll(r,i)))); !eq22;

@for(h(j):s(j)= (lambda1(j)/mu) /(1-lambda1(j)/mu)-(((1+kk)*((lambda1(j)/mu)^(kk+1)))/(1-(lambda1(j)/mu)^(kk+1)))); !eq23;


@for(c(i):lambda(i)=@sum(i1(j):@sum(i1(m):(w(j,m)+w(m,j))*z(j,i,i)))+ @sum(i1(k):@sum(h(ii):g(k,ii,i))) ); !eq24;


@for(h(jj):lambda(jj)=@sum(i1(k):@sum(c(j):g(k,jj,j)))); !eq25;



@for(i1(k):@for(i1(i):@for(i1(r):@bin(v(k,r,i)))));
@for(i1(k):@for(c(i):@bin(ll(k,i))));
@for(i1(k):@for(i1(i):@for(i1(r):@bin(u(k,r,i)))));
@for(h(k):@for(h(i):@for(h(r):@bin(z(k,r,i)))));
@for(i1(k):@for(i1(i):@for(i1(r):@bin(y(k,r,i)))));
end


