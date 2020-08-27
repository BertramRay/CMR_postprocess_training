%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 心肌高亮梗死区域图像分割任务 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ！！！代码运行快速上手步骤！！！
% 1.配置【导入数据】部分的文件地址，注意图像数据和Mask数据的区别
% 2.点击运行后将看到梗死区域图像顺次播放
% 3.可在【参数设置】和【图像修正方法开关】部分对程序进行精确调控
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 清空工作区
clear all;close all;clc;

%%%%%%%%%%%%%%%%%%%
%    导入数据     %
%%%%%%%%%%%%%%%%%%%
% 导入图像数据与Mask数据 A：图像数据 B：Mask数据
% 实验对象C139
A = load('D:\【学习】\SRT\PostProcessing_data_demo\PostProcessing_data_demo\PSIR\VTA_Pig_C139_2_days_post_MI0_WIP_LHi-INAV-PSIR+T25ms_SENSE_41_1');
B = load('D:\【学习】\SRT\PostProcessing_data_demo\PostProcessing_data_demo\ROI_Seg3D_VTA_Pig_C139_2_days_post_MI0_WIP_LHi-INAV-PSIR+T25ms_SENSE_41_1\LV');
% 实验对象B306
% A = load('D:\【学习】\SRT\PostProcessing_data(1)\PostProcessing_data\20110811 VTA Pig B306 AcuteMI 3 Days Post\PSIR\VTA_Pig_B306_v2_WIP_Hi-INAV-PSIR+T25ms_SENSE_36_1');
% B = load('D:\【学习】\SRT\PostProcessing_data(1)\PostProcessing_data\20110811 VTA Pig B306 AcuteMI 3 Days Post\PSIR\XOR_MaskLayer');

%%%%%%%%%%%%%%%%%%%
%   全局变量定义  %
%%%%%%%%%%%%%%%%%%%
global flags;
global Mask;
global level;
global xlength;
global ylength;

%%%%%%%%%%%%%%%%%%%
%    参数设置     %
%%%%%%%%%%%%%%%%%%%
% 图像的长，宽
xlength = 320;
ylength = 320;
% n-SD法常数
n = 5;
% 生长法去噪阈值
groupnum = 60;
% 边缘试探法扩展方块数量
level = 4;
% n-SD法所得的参考区域平均值与标准差
valSD = 12;
valMean = 75;
% 种子点坐标
seedx = 190;
seedy = 142;
% [seedx,seedy] = ginput(1);
% 种子点作用半径
seedRadius = 50;
% 梯度作用半径
gradientRadius = 50;
% 选择Mask的起始与结束层面,默认1和60
picstart = 1;
picend = 60;
% 梗死体积与ROI体积初始化
MaskVolumn = 0;
LGEVolumn = 0;

%%%%%%%%%%%%%%%%%%%
% 图像修正方法开关%
%%%%%%%%%%%%%%%%%%%
% 注：1=on,0=off
% 梯度信息加权
method0 = 1;
% 生长法去噪
method1 = 1;
% 闭合法去黑核
method2 = 1;
% 半闭合法去黑核
method3 = 1;
% 种子点法消除远端伪影
method4 = 1;

%%%%%%%%%%%%%%%%%%%
%  正式处理部分   %
%%%%%%%%%%%%%%%%%%%
for picnum=picstart:picend
    %%%%%%%%%%%%%%%%%%%
    %   ROI区域提取   %
    %%%%%%%%%%%%%%%%%%%
    % 图像数据，MASK数据提取
    Image = A.I.M(:,:,picnum); 
    Mask = B.scirunnrrd.data(:,:,picnum);
    MaskCopy = Mask;
    MaskArea = 0;%ROI区域面积
    for i=1:xlength
        for j=1:ylength
            if Mask(i,j) == 1 
                MaskArea = MaskArea + 1;%计算ROI面积
            else
                Image(i,j) = 0;%非ROI区域填充为黑色背景区域
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%
    %  梯度信息提取   %
    %%%%%%%%%%%%%%%%%%%
    if method0 == 1
        gradient = zeros(xlength,ylength);
        %计算每个像素点与上下左右四个像素点的绝对值的差的和
        for i=1:xlength
            for j=1:ylength
                if Mask(i,j) == 1
                    if Mask(i+1,j) == 1
                        gradient(i,j) = gradient(i,j)+abs(Image(i+1,j)-Image(i,j));
                    end
                    if Mask(i-1,j) == 1
                        gradient(i,j) = gradient(i,j)+abs(Image(i-1,j)-Image(i,j));
                    end
                    if Mask(i,j+1) == 1
                        gradient(i,j) = gradient(i,j)+abs(Image(i,j+1)-Image(i,j));
                    end
                    if Mask(i,j-1) == 1
                        gradient(i,j) = gradient(i,j)+abs(Image(i,j-1)-Image(i,j));
                    end
                end
            end
        end
        % 将像素点与种子点的距离纳入梯度信息加权的考量范畴，离种子点越近，梯度作用越大
        for i=1:xlength
            for j=1:ylength
                dis = abs(seedx-i)+abs(seedy-j);
                if dis<gradientRadius
                    Image(i,j) = (gradientRadius-dis)/gradientRadius*gradient(i,j)+Image(i,j);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%
    %  n-SD法初步分割 %
    %%%%%%%%%%%%%%%%%%%
    for i=1:xlength
        for j=1:ylength
            if Mask(i,j) == 1
                if Image(i,j) <= valMean + n*valSD
                        Mask(i,j) = 0;%将低灰度值区域填充为黑色背景区域
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%
    %    图像修正     %
    %%%%%%%%%%%%%%%%%%%
    % 1.生长法去除边缘高亮噪声
    if method1 == 1
        flags = zeros(xlength,ylength);
        for i=1:xlength
            for j=1:ylength
                if flags(i,j) == 0 %未扩展
                    if Mask(i,j) == 1
                        sum = growup(i,j);
                        if sum< groupnum
                            delete(i,j);
                        end
                    end
                end
            end
        end
    end
    % 2.闭合法去黑核
    if method2 == 1
        fillNucleus();
    end
    % 3.半闭合法去黑核
    if method3 == 1
        flags = zeros(xlength,ylength);
        for i=1:xlength
            for j=1:ylength
                if Mask(i,j) == 1
                    flags(i,j) = 1;
                end
            end
        end
        for i=1:xlength
            for j=1:ylength
                if flags(i,j) == 1
                    flat(i,j);
                end
            end
        end
        fillNucleus();
    end
    % 4.种子点法去除远端伪影
    if method4 == 1
        for i=1:xlength
            for j=1:ylength
                if abs(seedx-i)+abs(seedy-j)>seedRadius
                    Mask(i,j) = 0;
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%
    %    结果展示     %
    %%%%%%%%%%%%%%%%%%%
    %不同层面的梗死区域
    h = figure(1);
    Mask = Mask+MaskCopy;
    imshow(Mask,[],'InitialMagnification','fit');
    
    %%%%%%%%%%%%%%%%%%%
    %    数据累计     %
    %%%%%%%%%%%%%%%%%%%
    MaskVolumn = MaskVolumn + MaskArea;
    LGEArea = 0;
    for i=1:xlength
        for j=1:ylength
            if Mask(i,j) == 1
                LGEArea = LGEArea+1;
            end
        end
    end
    LGEVolumn = LGEVolumn + LGEArea;
end

% 梗死比例计算
LGEratio = LGEVolumn/MaskVolumn;



% 区域生长函数
function sum = growup(i,j)
global flags;
global Mask;
    sum = 1;
    if flags(i-1,j) == 0 %未扩展
        flags(i-1,j) = 1;
        if Mask(i-1,j) == 1
            sum = sum + growup(i-1,j);
        end
    end
    if flags(i+1,j) == 0 %未扩展
        flags(i+1,j) = 1;
        if Mask(i+1,j) == 1
            sum = sum + growup(i+1,j);
        end
    end
    if flags(i,j-1) == 0 %未扩展
        flags(i,j-1) = 1;
        if Mask(i,j-1) == 1
            sum = sum + growup(i,j-1);
        end
    end
    if flags(i,j+1) == 0 %未扩展
        flags(i,j+1) = 1;
        if Mask(i,j+1) == 1
            sum = sum + growup(i,j+1);
        end
    end
end

% 区域删除函数
function delete(i,j)
global Mask;
    Mask(i,j) = 0;
    if Mask(i-1,j) == 1
        delete(i-1,j);
    end
    if Mask(i+1,j) == 1
        delete(i+1,j);
    end
    if Mask(i,j-1) == 1
        delete(i,j-1);
    end
    if Mask(i,j+1) == 1
        delete(i,j+1);
    end
end

% 闭合法去黑核函数
function fillNucleus()
global Mask;
global flags;
global xlength;
global ylength;
    flags = zeros(xlength,ylength);
    flags(1,1) = 1;
    %标记黑色背景
    for k=1:3%扩展三次，确保完全覆盖
        for i=1:xlength
            for j=1:ylength
                if i-1>=1
                    if flags(i-1,j) == 1
                        if Mask(i-1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end        
            end
        end
        for j=1:ylength
            for i=1:xlength
                if i-1>=1
                    if flags(i-1,j) == 1
                        if Mask(i-1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end        
            end
        end
        for i=xlength:-1:1
            for j=ylength:-1:1
                if i-1>=1
                    if flags(i-1,j) == 1
                        if Mask(i-1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end        
            end
        end
        for j=ylength:-1:1
            for i=xlength:-1:1
                if i-1>=1
                    if flags(i-1,j) == 1
                        if Mask(i-1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%黑色背景 
                            flags(i,j) = 1;
                        end
                    end
                end        
            end
        end
    end
    for i=1:xlength
        for j=1:ylength
            if flags(i,j) == 0
                if Mask(i,j) == 0 %未标记黑色区域，即黑核
                    Mask(i,j) = 1;
                end
            end
        end
    end
end
% 专家阈值去黑核扩展函数
function flat(i,j)
global Mask;
global flags;
global level;
    if flags(i+1,j) == 0 %0为未标记黑块，1为白块，2为扩展的边缘白块,3连接的边缘白块
        flags(i+1,j) = 2;
        for b=j-1:j+1
            if flags(i+2,b) == 1
                Mask(i+1,j) = 1;
                flags(i+1,j) = 3;
                break;
            end
            if flags(i+2,b) == 2
                Mask(i+1,j) = 1;
                Mask(i+2,b) = 1;
                flags(i+1,j) = 3;
                flags(i+2,b) = 3;
                break;
            end
        end
        for len = 3:level
            if flags(i+len,j) == 2
                for idx = 1:len
                    Mask(i+idx,j) = 1;
                    flags(i+idx,j) = 3;
                end
                break;
            end
        end
    end
    if flags(i-1,j) == 0 %0为未标记黑块，1为白块，2为扩展的边缘白块
        flags(i-1,j) = 2;
        for b=j-1:j+1
            if flags(i-2,b) == 1
                Mask(i-1,j) = 1;
                flags(i-1,j) = 3;
                break;
            end
            if flags(i-2,b) == 2
                Mask(i-1,j) = 1;
                Mask(i-2,b) = 1;
                flags(i-1,j) = 3;
                flags(i-2,b) = 3;
                break;
            end
        end
        for len = 3:level
            if flags(i-len,j) == 2
                for idx = 1:len
                    Mask(i-idx,j) = 1;
                    flags(i-idx,j) = 3;
                end
                break;
            end
        end
    end
    if flags(i,j+1) == 0 %0为未标记黑块，1为白块，2为扩展的边缘白块
        flags(i,j+1) = 2;
        for a=i-1:i+1
            if flags(a,j+2) == 1
                Mask(i,j+1) = 1;
                flags(i,j+1) = 3;
                break;
            end
            if flags(a,j+2) == 2
                Mask(i,j+1) = 1;
                Mask(a,j+2) = 1;
                flags(i,j+1) = 3;
                flags(a,j+2) = 3;
                break;
            end
        end
        for len = 3:level
            if flags(i,j+len) == 2
                for idx = 1:len
                    Mask(i,j+idx) = 1;
                    flags(i,j+idx) = 3;
                end
                break;
            end
        end
    end
    if flags(i,j-1) == 0 %0为未标记黑块，1为白块，2为扩展的边缘白块
        flags(i,j-1) = 2;
        for a=i-1:i+1
            if flags(a,j-2) == 1
                Mask(i,j-1) = 1;
                flags(i,j-1) = 3;
                break;
            end
            if flags(a,j-2) == 2
                Mask(i,j-1) = 1;
                Mask(a,j-2) = 1;
                flags(i,j-1) = 3;
                flags(a,j-2) = 3;
                break;
            end
        end
        for len = 3:level
            if flags(i,j-len) == 2
                for idx = 1:len
                    Mask(i,j-idx) = 1;
                    flags(i,j-idx) = 3;
                end
                break;
            end
        end
    end
end