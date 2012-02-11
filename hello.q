/the value of PI
PI:{2*asin 1}[];
/generate two independent normal distribution series
genNorm:{((cos;sin)@\:2*PI*y)*\:sqrt -2*log x};
/generate two price time series - `T`dt`sigma`s0
simPrice:{[a] steps:`int$a[`T]%a`dt;a[`price]+sums a[`sigma]*a[`price]*(sqrt a[`dt])*genNorm[steps?1f;steps?1f]};

/fit an autoregressive time series model to the data by OLS, returns the parameter vector
arOLS:{[x;p;i] X:{[x;p;y] p _ y xprev x}[x;p;]@/:1+til p;Y:p _ x;if[i;X,:(count Y)#1f];Y lsq X};

/generate log normal distribution
logNorm:{[m;v;x] mu:log[(m2)%sqrt[v+m2:m*m]]; sigma:sqrt[log 1+v%m2];:exp(mu+sigma*(sqrt[-2*log x?1f])*cos(2*PI*x?1f))};

/algorithm poisson random number (Knuth)
knuth:{[x;s]system"S ",string s;i:count x;k:i#0;p:i#1f;L:exp neg[x];while[any g:p>L;p:@[p;h:where g;*;(sum g)?1f];k:@[k;h;+;1]];:k-1};

/ar 1
{[yhat;rho;e]e+rho*yhat}\'[0;];