%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �ļ�������������ͼ��ָ����� %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������������п������ֲ��裡����
% 1.���á��������ݡ����ֵ��ļ���ַ��ע��ͼ�����ݺ�Mask���ݵ�����
% 2.������к󽫿�����������ͼ��˳�β���
% 3.���ڡ��������á��͡�ͼ�������������ء����ֶԳ�����о�ȷ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ��չ�����
clear all;close all;clc;

%%%%%%%%%%%%%%%%%%%
%    ��������     %
%%%%%%%%%%%%%%%%%%%
% ����ͼ��������Mask���� A��ͼ������ B��Mask����
% ʵ�����C139
A = load('D:\��ѧϰ��\SRT\PostProcessing_data_demo\PostProcessing_data_demo\PSIR\VTA_Pig_C139_2_days_post_MI0_WIP_LHi-INAV-PSIR+T25ms_SENSE_41_1');
B = load('D:\��ѧϰ��\SRT\PostProcessing_data_demo\PostProcessing_data_demo\ROI_Seg3D_VTA_Pig_C139_2_days_post_MI0_WIP_LHi-INAV-PSIR+T25ms_SENSE_41_1\LV');
% ʵ�����B306
% A = load('D:\��ѧϰ��\SRT\PostProcessing_data(1)\PostProcessing_data\20110811 VTA Pig B306 AcuteMI 3 Days Post\PSIR\VTA_Pig_B306_v2_WIP_Hi-INAV-PSIR+T25ms_SENSE_36_1');
% B = load('D:\��ѧϰ��\SRT\PostProcessing_data(1)\PostProcessing_data\20110811 VTA Pig B306 AcuteMI 3 Days Post\PSIR\XOR_MaskLayer');

%%%%%%%%%%%%%%%%%%%
%   ȫ�ֱ�������  %
%%%%%%%%%%%%%%%%%%%
global flags;
global Mask;
global level;
global xlength;
global ylength;

%%%%%%%%%%%%%%%%%%%
%    ��������     %
%%%%%%%%%%%%%%%%%%%
% ͼ��ĳ�����
xlength = 320;
ylength = 320;
% n-SD������
n = 5;
% ������ȥ����ֵ
groupnum = 60;
% ��Ե��̽����չ��������
level = 4;
% n-SD�����õĲο�����ƽ��ֵ���׼��
valSD = 12;
valMean = 75;
% ���ӵ�����
seedx = 190;
seedy = 142;
% [seedx,seedy] = ginput(1);
% ���ӵ����ð뾶
seedRadius = 50;
% �ݶ����ð뾶
gradientRadius = 50;
% ѡ��Mask����ʼ���������,Ĭ��1��60
picstart = 1;
picend = 60;
% ���������ROI�����ʼ��
MaskVolumn = 0;
LGEVolumn = 0;

%%%%%%%%%%%%%%%%%%%
% ͼ��������������%
%%%%%%%%%%%%%%%%%%%
% ע��1=on,0=off
% �ݶ���Ϣ��Ȩ
method0 = 1;
% ������ȥ��
method1 = 1;
% �պϷ�ȥ�ں�
method2 = 1;
% ��պϷ�ȥ�ں�
method3 = 1;
% ���ӵ㷨����Զ��αӰ
method4 = 1;

%%%%%%%%%%%%%%%%%%%
%  ��ʽ������   %
%%%%%%%%%%%%%%%%%%%
for picnum=picstart:picend
    %%%%%%%%%%%%%%%%%%%
    %   ROI������ȡ   %
    %%%%%%%%%%%%%%%%%%%
    % ͼ�����ݣ�MASK������ȡ
    Image = A.I.M(:,:,picnum); 
    Mask = B.scirunnrrd.data(:,:,picnum);
    MaskCopy = Mask;
    MaskArea = 0;%ROI�������
    for i=1:xlength
        for j=1:ylength
            if Mask(i,j) == 1 
                MaskArea = MaskArea + 1;%����ROI���
            else
                Image(i,j) = 0;%��ROI�������Ϊ��ɫ��������
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%
    %  �ݶ���Ϣ��ȡ   %
    %%%%%%%%%%%%%%%%%%%
    if method0 == 1
        gradient = zeros(xlength,ylength);
        %����ÿ�����ص������������ĸ����ص�ľ���ֵ�Ĳ�ĺ�
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
        % �����ص������ӵ�ľ��������ݶ���Ϣ��Ȩ�Ŀ������룬�����ӵ�Խ�����ݶ�����Խ��
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
    %  n-SD�������ָ� %
    %%%%%%%%%%%%%%%%%%%
    for i=1:xlength
        for j=1:ylength
            if Mask(i,j) == 1
                if Image(i,j) <= valMean + n*valSD
                        Mask(i,j) = 0;%���ͻҶ�ֵ�������Ϊ��ɫ��������
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%
    %    ͼ������     %
    %%%%%%%%%%%%%%%%%%%
    % 1.������ȥ����Ե��������
    if method1 == 1
        flags = zeros(xlength,ylength);
        for i=1:xlength
            for j=1:ylength
                if flags(i,j) == 0 %δ��չ
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
    % 2.�պϷ�ȥ�ں�
    if method2 == 1
        fillNucleus();
    end
    % 3.��պϷ�ȥ�ں�
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
    % 4.���ӵ㷨ȥ��Զ��αӰ
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
    %    ���չʾ     %
    %%%%%%%%%%%%%%%%%%%
    %��ͬ����Ĺ�������
    h = figure(1);
    Mask = Mask+MaskCopy;
    imshow(Mask,[],'InitialMagnification','fit');
    
    %%%%%%%%%%%%%%%%%%%
    %    �����ۼ�     %
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

% ������������
LGEratio = LGEVolumn/MaskVolumn;



% ������������
function sum = growup(i,j)
global flags;
global Mask;
    sum = 1;
    if flags(i-1,j) == 0 %δ��չ
        flags(i-1,j) = 1;
        if Mask(i-1,j) == 1
            sum = sum + growup(i-1,j);
        end
    end
    if flags(i+1,j) == 0 %δ��չ
        flags(i+1,j) = 1;
        if Mask(i+1,j) == 1
            sum = sum + growup(i+1,j);
        end
    end
    if flags(i,j-1) == 0 %δ��չ
        flags(i,j-1) = 1;
        if Mask(i,j-1) == 1
            sum = sum + growup(i,j-1);
        end
    end
    if flags(i,j+1) == 0 %δ��չ
        flags(i,j+1) = 1;
        if Mask(i,j+1) == 1
            sum = sum + growup(i,j+1);
        end
    end
end

% ����ɾ������
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

% �պϷ�ȥ�ں˺���
function fillNucleus()
global Mask;
global flags;
global xlength;
global ylength;
    flags = zeros(xlength,ylength);
    flags(1,1) = 1;
    %��Ǻ�ɫ����
    for k=1:3%��չ���Σ�ȷ����ȫ����
        for i=1:xlength
            for j=1:ylength
                if i-1>=1
                    if flags(i-1,j) == 1
                        if Mask(i-1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%��ɫ���� 
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
                        if Mask(i-1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%��ɫ���� 
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
                        if Mask(i-1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%��ɫ���� 
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
                        if Mask(i-1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if i+1<=xlength
                    if flags(i+1,j) == 1
                        if Mask(i+1,j) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j-1>=1
                    if flags(i,j-1) == 1
                        if Mask(i,j-1) == 0%��ɫ���� 
                            flags(i,j) = 1;
                        end
                    end
                end
                if j+1<=ylength
                    if flags(i,j+1) == 1
                        if Mask(i,j+1) == 0%��ɫ���� 
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
                if Mask(i,j) == 0 %δ��Ǻ�ɫ���򣬼��ں�
                    Mask(i,j) = 1;
                end
            end
        end
    end
end
% ר����ֵȥ�ں���չ����
function flat(i,j)
global Mask;
global flags;
global level;
    if flags(i+1,j) == 0 %0Ϊδ��Ǻڿ飬1Ϊ�׿飬2Ϊ��չ�ı�Ե�׿�,3���ӵı�Ե�׿�
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
    if flags(i-1,j) == 0 %0Ϊδ��Ǻڿ飬1Ϊ�׿飬2Ϊ��չ�ı�Ե�׿�
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
    if flags(i,j+1) == 0 %0Ϊδ��Ǻڿ飬1Ϊ�׿飬2Ϊ��չ�ı�Ե�׿�
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
    if flags(i,j-1) == 0 %0Ϊδ��Ǻڿ飬1Ϊ�׿飬2Ϊ��չ�ı�Ե�׿�
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