## 常量

```c
lx     = single(200);        % x length of computational domain
ly     = single(150);        % y length of computational domain
lxy    = lx*ly;
cc     = single(1);          % lattice speeds
c_squ  = cc*cc/3;            % square of sound speed
ex     = single([0 1 0 -1 0 1 -1 -1 1]);                % velocity scheme
ey     = single([0 0 1 0 -1 1 1 -1 -1]);
M      = single([ 1  1  1  1  1  1  1  1  1;
                         -4 -1 -1 -1 -1  2  2  2  2;
                          4 -2 -2 -2 -2  1  1  1  1;
                          0  1  0 -1  0  1 -1 -1  1;
                          0 -2  0  2  0  1 -1 -1  1;
                          0  0  1  0 -1  1  1 -1 -1;
                          0  0 -2  0  2  1  1 -1 -1;
                          0  1 -1  1 -1  0  0  0  0;
                          0  0  0  0  0  1 -1  1 -1 ]);
Minv   = single([ 4 -4  4  0  0  0  0  0  0;   % inverse matrix
                             4 -1 -2  6 -6  0  0  9  0;
                             4 -1 -2  0  0  6 -6 -9  0;
                             4 -1 -2 -6  6  0  0  9  0;
                             4 -1 -2  0  0 -6  6 -9  0;
                             4  2  1  6  3  6  3  0  9;
                             4  2  1 -6 -3  6  3  0 -9;
                             4  2  1 -6 -3 -6 -3  0  9;
                             4  2  1  6  3 -6 -3  0 -9])/36;
sigm   = single(1.2);
```


## 变量

## 函数
