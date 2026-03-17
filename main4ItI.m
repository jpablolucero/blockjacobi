clear;

k = 2;
eta = k;
[u_exact, rhs, c0, ~, ~, IBC, ~, ~] = DtNTest1(k);
% [u_exact,rhs,c0,IBC,poincareSteklovOperator] = ItITest1(k);
poincareSteklovOperator = "ItI";

for div = 5:5
    sSW = Subdomain(div, rhs, 0.0, 0.5, 0.0, 0.5, eta, c0, poincareSteklovOperator, @get_sem);
    sSE = Subdomain(div, rhs, 0.5, 1.0, 0.0, 0.5, eta, c0, poincareSteklovOperator, @get_sem);
    sNW = Subdomain(div, rhs, 0.0, 0.5, 0.5, 1.0, eta, c0, poincareSteklovOperator, @get_sem);
    sNE = Subdomain(div, rhs, 0.5, 1.0, 0.5, 1.0, eta, c0, poincareSteklovOperator, @get_sem);

    m = 2^div - 1;
    cBL = 4*m+1; cBR = 4*m+2; cTL = 4*m+3; cTR = 4*m+4;

    % --- mc: mass matrix entry at interior corner (ItI only) ---

    % --- Apply boundary conditions ---
    sSW.setBoundaryCondition(IBC{1}(sSW.px(sSW.idx_left),   sSW.py(sSW.idx_left)),   1);
    sSW.setBoundaryCondition(IBC{3}(sSW.px(sSW.idx_bottom), sSW.py(sSW.idx_bottom)), 3);
    sSW.setBoundaryCondition([0.5*(IBC{1}(sSW.px(sSW.idx_boundary(cBL)),sSW.py(sSW.idx_boundary(cBL))) + IBC{3}(sSW.px(sSW.idx_boundary(cBL)),sSW.py(sSW.idx_boundary(cBL))));
                               0.5* IBC{3}(sSW.px(sSW.idx_boundary(cBR)),sSW.py(sSW.idx_boundary(cBR)));
                               0.5* IBC{1}(sSW.px(sSW.idx_boundary(cTL)),sSW.py(sSW.idx_boundary(cTL)));
                               0], 5);

    sSE.setBoundaryCondition(IBC{2}(sSE.px(sSE.idx_right),  sSE.py(sSE.idx_right)),  2);
    sSE.setBoundaryCondition(IBC{3}(sSE.px(sSE.idx_bottom), sSE.py(sSE.idx_bottom)), 3);
    sSE.setBoundaryCondition([0.5* IBC{3}(sSE.px(sSE.idx_boundary(cBL)),sSE.py(sSE.idx_boundary(cBL)));
                               0.5*(IBC{2}(sSE.px(sSE.idx_boundary(cBR)),sSE.py(sSE.idx_boundary(cBR))) + IBC{3}(sSE.px(sSE.idx_boundary(cBR)),sSE.py(sSE.idx_boundary(cBR))));
                               0;
                               0.5* IBC{2}(sSE.px(sSE.idx_boundary(cTR)),sSE.py(sSE.idx_boundary(cTR)))], 5);

    sNW.setBoundaryCondition(IBC{1}(sNW.px(sNW.idx_left),   sNW.py(sNW.idx_left)),   1);
    sNW.setBoundaryCondition(IBC{4}(sNW.px(sNW.idx_top),    sNW.py(sNW.idx_top)),    4);
    sNW.setBoundaryCondition([0.5* IBC{1}(sNW.px(sNW.idx_boundary(cBL)),sNW.py(sNW.idx_boundary(cBL)));
                               0;
                               0.5*(IBC{1}(sNW.px(sNW.idx_boundary(cTL)),sNW.py(sNW.idx_boundary(cTL))) + IBC{4}(sNW.px(sNW.idx_boundary(cTL)),sNW.py(sNW.idx_boundary(cTL))));
                               0.5* IBC{4}(sNW.px(sNW.idx_boundary(cTR)),sNW.py(sNW.idx_boundary(cTR)))], 5);

    sNE.setBoundaryCondition(IBC{2}(sNE.px(sNE.idx_right),  sNE.py(sNE.idx_right)),  2);
    sNE.setBoundaryCondition(IBC{4}(sNE.px(sNE.idx_top),    sNE.py(sNE.idx_top)),    4);
    sNE.setBoundaryCondition([0;
                               0.5* IBC{2}(sNE.px(sNE.idx_boundary(cBR)),sNE.py(sNE.idx_boundary(cBR)));
                               0.5* IBC{4}(sNE.px(sNE.idx_boundary(cTL)),sNE.py(sNE.idx_boundary(cTL)));
                               0.5*(IBC{2}(sNE.px(sNE.idx_boundary(cTR)),sNE.py(sNE.idx_boundary(cTR))) + IBC{4}(sNE.px(sNE.idx_boundary(cTR)),sNE.py(sNE.idx_boundary(cTR))))], 5);

    N = 8*m + 12;

    A = [ %                 1                                   2                                  3                                   4                               5                                    6                                   7                                  8                                   9                                             10                                            11                                        12                         13                                             14                                             15                                    16                                    17                                      18                                          19                                   20                                          
          sSW.T{2,2}                        ,  sSW.T{2,4}                       ,  speye(m)                        ,  sparse(m,m)                      ,  sparse(m,m)                     , sparse(m,m)                      ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sparse(m,1)                         ,  sparse(m,1)                          ,  sSW.T{2,5}(:,4)                      ,  0.5*sSW.T{2,5}(:,2)                 ,  sparse(m,1)                           ,  sparse(m,1)                          ,  0.5*sSW.T{2,5}(:,3)                   ,  sparse(m,1)                           ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  sparse(m,1)                          ;
          sSW.T{4,2}                        ,  sSW.T{4,4}                       ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sparse(m,m)                     , speye(m)                         ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sparse(m,1)                         ,  sparse(m,1)                          ,  sSW.T{4,5}(:,4)                      ,  0.5*sSW.T{4,5}(:,2)                 ,  sparse(m,1)                           ,  sparse(m,1)                          ,  0.5*sSW.T{4,5}(:,3)                   ,  sparse(m,1)                           ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  sparse(m,1)                          ;
          speye(m)                          ,  sparse(m,m)                      ,  sSE.T{1,1}                      ,  sSE.T{1,4}                       ,  sparse(m,m)                     , sparse(m,m)                      ,  sparse(m,m)                     ,  sparse(m,m)                      ,  0.5*sSE.T{1,5}(:,1)                 ,  sparse(m,1)                          ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sSE.T{1,5}(:,3)                       ,  sparse(m,1)                          ,  sparse(m,1)                           ,  sparse(m,1)                           ,  sparse(m,1)                         ,  sparse(m,1)                           ,  0.5*sSE.T{1,5}(:,4)                  ,  sparse(m,1)                          ;
          sparse(m,m)                       ,  sparse(m,m)                      ,  sSE.T{4,1}                      ,  sSE.T{4,4}                       ,  sparse(m,m)                     , sparse(m,m)                      ,  sparse(m,m)                     ,  speye(m)                         ,  0.5*sSE.T{4,5}(:,1)                 ,  sparse(m,1)                          ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sSE.T{4,5}(:,3)                       ,  sparse(m,1)                          ,  sparse(m,1)                           ,  sparse(m,1)                           ,  sparse(m,1)                         ,  sparse(m,1)                           ,  0.5*sSE.T{4,5}(:,4)                  ,  sparse(m,1)                          ;
          sparse(m,m)                       ,  sparse(m,m)                      ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sNW.T{2,2}                      , sNW.T{2,3}                       ,  speye(m)                        ,  sparse(m,m)                      ,  sparse(m,1)                         ,  0.5*sNW.T{2,5}(:,1)                  ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  sparse(m,1)                           ,  sNW.T{2,5}(:,2)                       ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  0.5*sNW.T{2,5}(:,4)                  ;
          sparse(m,m)                       ,  speye(m)                         ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sNW.T{3,2}                      , sNW.T{3,3}                       ,  sparse(m,m)                     ,  sparse(m,m)                      ,  sparse(m,1)                         ,  0.5*sNW.T{3,5}(:,1)                  ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  sparse(m,1)                           ,  sNW.T{3,5}(:,2)                       ,  sparse(m,1)                         ,  sparse(m,1)                           ,  sparse(m,1)                          ,  0.5*sNW.T{3,5}(:,4)                  ;
          sparse(m,m)                       ,  sparse(m,m)                      ,  sparse(m,m)                     ,  sparse(m,m)                      ,  speye(m)                        , sparse(m,m)                      ,  sNE.T{1,1}                      ,  sNE.T{1,3}                       ,  sparse(m,1)                         ,  sparse(m,1)                          ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sparse(m,1)                           ,  0.5*sNE.T{1,5}(:,2)                  ,  sparse(m,1)                           ,  sparse(m,1)                           ,  0.5*sNE.T{1,5}(:,3)                 ,  sNE.T{1,5}(:,1)                       ,  sparse(m,1)                          ,  sparse(m,1)                          ;
          sparse(m,m)                       ,  sparse(m,m)                      ,  sparse(m,m)                     ,  speye(m)                         ,  sparse(m,m)                     , sparse(m,m)                      ,  sNE.T{3,1}                      ,  sNE.T{3,3}                       ,  sparse(m,1)                         ,  sparse(m,1)                          ,  sparse(m,1)                          ,  sparse(m,1)                         ,  sparse(m,1)                           ,  0.5*sNE.T{3,5}(:,2)                  ,  sparse(m,1)                           ,  sparse(m,1)                           ,  0.5*sNE.T{3,5}(:,3)                 ,  sNE.T{3,5}(:,1)                       ,  sparse(m,1)                          ,  sparse(m,1)                          ;
          sSW.T{5,2}(2,:)                   ,  sSW.T{5,4}(2,:)                  ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sparse(1,m)                     , sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  1                                   ,  0                                    ,  sSW.T{5,5}(2,4)                      ,  0.5*sSW.T{5,5}(2,2)+0.5             ,  0                                     ,  0                                    ,  0.5*sSW.T{5,5}(2,3)                   ,  0                                     ,  0                                   ,  0                                     ,  0                                    ,  0                                    ;
          sSW.T{5,2}(3,:)                   ,  sSW.T{5,4}(3,:)                  ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sparse(1,m)                     , sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0                                   ,  1                                    ,  sSW.T{5,5}(3,4)                      ,  0.5*sSW.T{5,5}(3,2)                 ,  0                                     ,  0                                    ,  0.5*sSW.T{5,5}(3,3)+0.5               ,  0                                     ,  0                                   ,  0                                     ,  0                                    ,  0                                    ;
          -sSW.T{5,2}(4,:)/sSW.Mb(cTR,cTR)  ,  -sSW.T{5,4}(4,:)/sSW.Mb(cTR,cTR) ,  sSE.T{5,1}(3,:)/sSE.Mb(cTL,cTL) ,  sSE.T{5,4}(3,:)/sSE.Mb(cTL,cTL)  ,  sparse(1,m)                     , sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0.5*sSE.T{5,5}(3,1)/sSE.Mb(cTL,cTL) ,  0                                    ,  (1-sSW.T{5,5}(4,4))/sSW.Mb(cTR,cTR)  ,  -0.5*sSW.T{5,5}(4,2)/sSW.Mb(cTR,cTR),  (-1+sSE.T{5,5}(3,3))/sSE.Mb(cTL,cTL)  ,  0                                    ,  -0.5*sSW.T{5,5}(4,3)/sSW.Mb(cTR,cTR)  ,  0                                     ,  0                                   ,  0                                     ,  0.5*sSE.T{5,5}(3,4)/sSE.Mb(cTL,cTL)  ,  0                                    ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sSE.T{5,1}(1,:)                 ,  sSE.T{5,4}(1,:)                  ,  sparse(1,m)                     , sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0.5*sSE.T{5,5}(1,1)+0.5             ,  0                                    ,  0                                    ,  1                                   ,  sSE.T{5,5}(1,3)                       ,  0                                    ,  0                                     ,  0                                     ,  0                                   ,  0                                     ,  0.5*sSE.T{5,5}(1,4)                  ,  0                                    ;
          -sSW.T{5,2}(4,:)/sSW.Mb(cTR,cTR)  ,  -sSW.T{5,4}(4,:)/sSW.Mb(cTR,cTR) ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sNW.T{5,2}(2,:)/sNW.Mb(cBR,cBR) , sNW.T{5,3}(2,:)/sNW.Mb(cBR,cBR)  ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0                                   ,  0.5*sNW.T{5,5}(2,1)/sNW.Mb(cBR,cBR)  ,  (1-sSW.T{5,5}(4,4))/sSW.Mb(cTR,cTR)  ,  -0.5*sSW.T{5,5}(4,2)/sSW.Mb(cTR,cTR),  0                                     ,  0                                    ,  -0.5*sSW.T{5,5}(4,3)/sSW.Mb(cTR,cTR)  ,  (-1+sNW.T{5,5}(2,2))/sNW.Mb(cBR,cBR)  ,  0                                   ,  0                                     ,  0                                    ,  0.5*sNW.T{5,5}(2,4)/sNW.Mb(cBR,cBR)  ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sSE.T{5,1}(4,:)                 ,  sSE.T{5,4}(4,:)                  ,  sparse(1,m)                     , sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0.5*sSE.T{5,5}(4,1)                 ,  0                                    ,  0                                    ,  0                                   ,  sSE.T{5,5}(4,3)                       ,  1                                    ,  0                                     ,  0                                     ,  0                                   ,  0                                     ,  0.5*sSE.T{5,5}(4,4)+0.5              ,  0                                    ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sNW.T{5,2}(1,:)                 , sNW.T{5,3}(1,:)                  ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0                                   ,  0.5*sNW.T{5,5}(1,1)+0.5              ,  0                                    ,  0                                   ,  0                                     ,  0                                    ,  1                                     ,  sNW.T{5,5}(1,2)                       ,  0                                   ,  0                                     ,  0                                    ,  0.5*sNW.T{5,5}(1,4)                  ;
          -sSW.T{5,2}(4,:)/sSW.Mb(cTR,cTR)  ,  -sSW.T{5,4}(4,:)/sSW.Mb(cTR,cTR) ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sparse(1,m)                     , sparse(1,m)                      ,  sNE.T{5,1}(1,:)/sNE.Mb(cBL,cBL) ,  sNE.T{5,3}(1,:)/sNE.Mb(cBL,cBL)  ,  0                                   ,  0                                    ,  (1-sSW.T{5,5}(4,4))/sSW.Mb(cTR,cTR)  ,  -0.5*sSW.T{5,5}(4,2)/sSW.Mb(cTR,cTR),  0                                     ,  0.5*sNE.T{5,5}(1,2)/sNE.Mb(cBL,cBL)  ,  -0.5*sSW.T{5,5}(4,3)/sSW.Mb(cTR,cTR)  ,  0                                     ,  0.5*sNE.T{5,5}(1,3)/sNE.Mb(cBL,cBL) ,  (-1+sNE.T{5,5}(1,1))/sNE.Mb(cBL,cBL)  ,  0                                    ,  0                                    ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sNW.T{5,2}(4,:)                 , sNW.T{5,3}(4,:)                  ,  sparse(1,m)                     ,  sparse(1,m)                      ,  0                                   ,  0.5*sNW.T{5,5}(4,1)                  ,  0                                    ,  0                                   ,  0                                     ,  0                                    ,  0                                     ,  sNW.T{5,5}(4,2)                       ,  1                                   ,  0                                     ,  0                                    ,  0.5*sNW.T{5,5}(4,4)+0.5              ;
          sSW.T{5,2}(4,:)                   ,  sSW.T{5,4}(4,:)                  ,  sSE.T{5,1}(3,:)                 ,  sSE.T{5,4}(3,:)                  ,  sNW.T{5,2}(2,:)                 , sNW.T{5,3}(2,:)                  ,  sNE.T{5,1}(1,:)                 ,  sNE.T{5,3}(1,:)                  ,  0.5*sSE.T{5,5}(3,1)                 ,  0.5*sNW.T{5,5}(2,1)                  ,  1+sSW.T{5,5}(4,4)                    ,  0.5*sSW.T{5,5}(4,2)                 ,  1+sSE.T{5,5}(3,3)                     ,  0.5*sNE.T{5,5}(1,2)                  ,  0.5*sSW.T{5,5}(4,3)                   ,  1+sNW.T{5,5}(2,2)                     ,  0.5*sNE.T{5,5}(1,3)                 ,  1+sNE.T{5,5}(1,1)                     ,  0.5*sSE.T{5,5}(3,4)                  ,  0.5*sNW.T{5,5}(2,4)                  ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sparse(1,m)                     , sparse(1,m)                      ,  sNE.T{5,1}(2,:)                 ,  sNE.T{5,3}(2,:)                  ,  0                                   ,  0                                    ,  0                                    ,  0                                   ,  0                                     ,  0.5*sNE.T{5,5}(2,2)+0.5              ,  0                                     ,  0                                     ,  0.5*sNE.T{5,5}(2,3)                 ,  sNE.T{5,5}(2,1)                       ,  1                                    ,  0                                    ;
          sparse(1,m)                       ,  sparse(1,m)                      ,  sparse(1,m)                     ,  sparse(1,m)                      ,  sparse(1,m)                     , sparse(1,m)                      ,  sNE.T{5,1}(3,:)                 ,  sNE.T{5,3}(3,:)                  ,  0                                   ,  0                                    ,  0                                    ,  0                                   ,  0                                     ,  0.5*sNE.T{5,5}(3,2)                  ,  0                                     ,  0                                     ,  0.5*sNE.T{5,5}(3,3)+0.5             ,  sNE.T{5,5}(3,1)                       ,  0                                    ,  1                                    ;
        ];

    b = [
        -sSW.h{2};
        -sSW.h{4};
        -sSE.h{1};
        -sSE.h{4};
        -sNW.h{2};
        -sNW.h{3};
        -sNE.h{1};
        -sNE.h{3};
        -(sSW.h{5}(2) - 0.5* IBC{3}(sSW.px(sSW.idx_boundary(cBR)),sSW.py(sSW.idx_boundary(cBR))));
        -(sSW.h{5}(3) - 0.5* IBC{1}(sSW.px(sSW.idx_boundary(cTL)),sSW.py(sSW.idx_boundary(cTL))));
        sSW.h{5}(4) / sSW.Mb(cTR,cTR) - sSE.h{5}(3) / sSE.Mb(cTL,cTL);
        -(sSE.h{5}(1) - 0.5* IBC{3}(sSE.px(sSE.idx_boundary(cBL)),sSE.py(sSE.idx_boundary(cBL))));
        sSW.h{5}(4) / sSW.Mb(cTR,cTR) - sNW.h{5}(2) / sNW.Mb(cBR,cBR);
        -(sSE.h{5}(4) - 0.5* IBC{2}(sSE.px(sSE.idx_boundary(cTR)),sSE.py(sSE.idx_boundary(cTR))));
        -(sNW.h{5}(1) - 0.5* IBC{1}(sNW.px(sNW.idx_boundary(cBL)),sNW.py(sNW.idx_boundary(cBL))));
        sSW.h{5}(4) / sSW.Mb(cTR,cTR) - sNE.h{5}(1) / sNE.Mb(cBL,cBL);
        -(sNW.h{5}(4) - 0.5* IBC{4}(sNW.px(sNW.idx_boundary(cTR)),sNW.py(sNW.idx_boundary(cTR))));
        -(sSW.h{5}(4) + sSE.h{5}(3) + sNW.h{5}(2) + sNE.h{5}(1));
        -(sNE.h{5}(2) - 0.5* IBC{2}(sNE.px(sNE.idx_boundary(cBR)),sNE.py(sNE.idx_boundary(cBR))));
        -(sNE.h{5}(3) - 0.5* IBC{4}(sNE.px(sNE.idx_boundary(cTL)),sNE.py(sNE.idx_boundary(cTL))));
    ];

    x = A \ b;

    ooSW = cell2mat(sSW.T) * [...
        zeros(m,1); 
        x(1:m);
        zeros(m,1);
        x(m+1:2*m); 
        0;             
        0.5*x(8*m+4);  
        0.5*x(8*m+7); 
        x(8*m+3)       ] + cell2mat(sSW.h);
    ooSE = cell2mat(sSE.T) * [...
	x(2*m+1:3*m); 
        zeros(m,1);   
        zeros(m,1); 
        x(3*m+1:4*m);
        0.5*x(8*m+1); 
        0;           
        x(8*m+5);    
        0.5*x(8*m+11)  ] + cell2mat(sSE.h);
    ooNW = cell2mat(sNW.T) * [...
        zeros(m,1); 
        x(4*m+1:5*m);
        x(5*m+1:6*m);
        zeros(m,1); 
        0.5*x(8*m+2); 
        x(8*m+8);   
        0;           
        0.5*x(8*m+12)  ] + cell2mat(sNW.h);
    ooNE = cell2mat(sNE.T) * [...
	x(6*m+1:7*m);
        zeros(m,1); 
        x(7*m+1:8*m);
        zeros(m,1);
        x(8*m+10);    
        0.5*x(8*m+6);
        0.5*x(8*m+9); 
        0             ] + cell2mat(sNE.h);

    ubSW = ([IBC{1}(sSW.px(sSW.idx_left),sSW.py(sSW.idx_left)); ...
             x(1:m); ...
             IBC{3}(sSW.px(sSW.idx_bottom),sSW.py(sSW.idx_bottom)); ...
             x(m+1:2*m); ...
             0.5*(IBC{1}(sSW.px(sSW.idx_boundary(cBL)),sSW.py(sSW.idx_boundary(cBL)))+IBC{3}(sSW.px(sSW.idx_boundary(cBL)),sSW.py(sSW.idx_boundary(cBL)))); ...
             0.5*IBC{3}(sSW.px(sSW.idx_boundary(cBR)),sSW.py(sSW.idx_boundary(cBR)))+0.5*x(8*m+4); ...
             0.5*IBC{1}(sSW.px(sSW.idx_boundary(cTL)),sSW.py(sSW.idx_boundary(cTL)))+0.5*x(8*m+7); ...
             x(8*m+3)] - ooSW) / (2i * eta);
    ubSE = ([x(2*m+1:3*m); ...
             IBC{2}(sSE.px(sSE.idx_right),sSE.py(sSE.idx_right)); ...
             IBC{3}(sSE.px(sSE.idx_bottom),sSE.py(sSE.idx_bottom)); ...
             x(3*m+1:4*m); ...
             0.5*IBC{3}(sSE.px(sSE.idx_boundary(cBL)),sSE.py(sSE.idx_boundary(cBL)))+0.5*x(8*m+1); ...
             0.5*(IBC{2}(sSE.px(sSE.idx_boundary(cBR)),sSE.py(sSE.idx_boundary(cBR)))+IBC{3}(sSE.px(sSE.idx_boundary(cBR)),sSE.py(sSE.idx_boundary(cBR)))); ...
             x(8*m+5); ...
             0.5*IBC{2}(sSE.px(sSE.idx_boundary(cTR)),sSE.py(sSE.idx_boundary(cTR)))+0.5*x(8*m+11)] - ooSE) / (2i * eta);
    ubNW = ([IBC{1}(sNW.px(sNW.idx_left),sNW.py(sNW.idx_left)); ...
             x(4*m+1:5*m); ...
             x(5*m+1:6*m); ...
             IBC{4}(sNW.px(sNW.idx_top),sNW.py(sNW.idx_top)); ...
             0.5*IBC{1}(sNW.px(sNW.idx_boundary(cBL)),sNW.py(sNW.idx_boundary(cBL)))+0.5*x(8*m+2); ...
             x(8*m+8); ...
             0.5*(IBC{1}(sNW.px(sNW.idx_boundary(cTL)),sNW.py(sNW.idx_boundary(cTL)))+IBC{4}(sNW.px(sNW.idx_boundary(cTL)),sNW.py(sNW.idx_boundary(cTL)))); ...
             0.5*IBC{4}(sNW.px(sNW.idx_boundary(cTR)),sNW.py(sNW.idx_boundary(cTR)))+0.5*x(8*m+12)] - ooNW) / (2i * eta);
    ubNE = ([x(6*m+1:7*m); ...
             IBC{2}(sNE.px(sNE.idx_right),sNE.py(sNE.idx_right)); ...
             x(7*m+1:8*m); ...
	     IBC{4}(sNE.px(sNE.idx_top),sNE.py(sNE.idx_top)); ...
             x(8*m+10); ...
             0.5*IBC{2}(sNE.px(sNE.idx_boundary(cBR)),sNE.py(sNE.idx_boundary(cBR)))+0.5*x(8*m+6); ...
             0.5*IBC{4}(sNE.px(sNE.idx_boundary(cTL)),sNE.py(sNE.idx_boundary(cTL)))+0.5*x(8*m+9); ...
             0.5*(IBC{2}(sNE.px(sNE.idx_boundary(cTR)),sNE.py(sNE.idx_boundary(cTR)))+IBC{4}(sNE.px(sNE.idx_boundary(cTR)),sNE.py(sNE.idx_boundary(cTR))))] - ooNE) / (2i * eta);

    uSW = zeros(size(sSW.rhs));
    uSE = zeros(size(sSE.rhs));
    uNW = zeros(size(sNW.rhs));
    uNE = zeros(size(sNE.rhs));

    uSW(sSW.idx_interior) = sSW.A \ (sSW.rhs(sSW.idx_interior) - sSW.B * ubSW);
    uSE(sSE.idx_interior) = sSE.A \ (sSE.rhs(sSE.idx_interior) - sSE.B * ubSE);
    uNW(sNW.idx_interior) = sNW.A \ (sNW.rhs(sNW.idx_interior) - sNW.B * ubNW);
    uNE(sNE.idx_interior) = sNE.A \ (sNE.rhs(sNE.idx_interior) - sNE.B * ubNE);

    uSW(sSW.idx_boundary) = ubSW;
    uSE(sSE.idx_boundary) = ubSE;
    uNW(sNW.idx_boundary) = ubNW;
    uNE(sNE.idx_boundary) = ubNE;

    uexSW = u_exact(sSW.px, sSW.py);
    uexSE = u_exact(sSE.px, sSE.py);
    uexNW = u_exact(sNW.px, sNW.py);
    uexNE = u_exact(sNE.px, sNE.py);

    num = norm(uSW - uexSW)^2 + norm(uSE - uexSE)^2 + norm(uNW - uexNW)^2 + norm(uNE - uexNE)^2;
    den = norm(uexSW)^2 + norm(uexSE)^2 + norm(uexNW)^2 + norm(uexNE)^2;

    fprintf("div: %i\n", div);
    fprintf("Rel L2 error = %.15e\n", sqrt(num / den));

    n = 2^div + 1;

    figure;
    hold on;
    surf(reshape(sSW.px, n, n), reshape(sSW.py, n, n), reshape(real(uSW), n, n));
    surf(reshape(sSE.px, n, n), reshape(sSE.py, n, n), reshape(real(uSE), n, n));
    surf(reshape(sNW.px, n, n), reshape(sNW.py, n, n), reshape(real(uNW), n, n));
    surf(reshape(sNE.px, n, n), reshape(sNE.py, n, n), reshape(real(uNE), n, n));
    view(3);
    hold off;
end
