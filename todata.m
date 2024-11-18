function [e_num,m,n,band,o2,o2_3d,m_turth,filename,RGBband] = todata(datanaumber)
number = datanaumber;
e_num = 4;
hdr = read_envihdr('./data/2019/roi'+string(number)+'.hdr');
Image = multibandread('./data/2019/roi'+string(number)+'.dat',hdr.size,[hdr.format '=>double'],hdr.header_offset,hdr.interleave,hdr.machine);
[n,m,band]=size(Image);
o2_3d_1=Image(:,:,1);o2_3d_1=mapminmax(o2_3d_1,0,1);
o2_3d_2=Image(:,:,2);o2_3d_2=mapminmax(o2_3d_2,0,1);
o2_3d_3=Image(:,:,3);o2_3d_3=mapminmax(o2_3d_3,0,1);
o2_3d_4=Image(:,:,4);o2_3d_4=mapminmax(o2_3d_4,0,1);
o2_3d_5=Image(:,:,5);o2_3d_5=mapminmax(o2_3d_5,0,1);
o2_3d_6=Image(:,:,6);o2_3d_6=mapminmax(o2_3d_6,0,1);
o2_3d_7=Image(:,:,7);o2_3d_7=mapminmax(o2_3d_7,0,1);
o2_3d(:,:,1)=o2_3d_1;o2_3d(:,:,2)=o2_3d_2;o2_3d(:,:,3)=o2_3d_3;o2_3d(:,:,4)=o2_3d_4;
o2_3d(:,:,5)=o2_3d_5;o2_3d(:,:,6)=o2_3d_6;o2_3d(:,:,7)=o2_3d_7;
o2 = hyperConvert2d(o2_3d);
%------ground truth-----
m_turth = 0;   %no ground truth
%------parameter------
band = 7;
filename = 'Beijing';
RGBband = [3,4,5];
end