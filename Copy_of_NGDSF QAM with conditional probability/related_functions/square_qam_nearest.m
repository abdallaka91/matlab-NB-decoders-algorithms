function [d, cndt] = square_qam_nearest(s, cnstl1)


Iw = min(real(cnstl1)):2:max(real(cnstl1));
Qw = (min(imag(cnstl1)):2:max(imag(cnstl1)));


r1 = real(s);
im1 = imag(s);

lbnd = [];
rbnd = [];
tbnd = [];
bbnd = [];



if r1 >= Iw(1)
    lbnd = Iw(find(Iw<=r1,1,'last'));
end


if im1 >= Qw(1)
    bbnd = Qw(find(Qw<=im1,1,'last'));
end

%----

if r1 < Iw(end)
    rbnd = Iw(find(Iw>=r1,1,'first'));
end


if im1 < Qw(end)
    tbnd = Qw(find(Qw>=im1,1,'first'));
end

p1 = [];p2=[];p3=[];p4=[];
if numel(rbnd)+numel(tbnd)==2
    p1 = rbnd + tbnd*1i;
end
if numel(rbnd)+numel(bbnd)==2
    p2 = rbnd + bbnd*1i;
end
if numel(lbnd)+numel(tbnd)==2
    p3 = lbnd + tbnd*1i;
end
if numel(lbnd)+numel(bbnd)==2
    p4 = lbnd + bbnd*1i;
end

cndt = [p1 p2 p3 p4];
d = abs(s-cndt);

% plot(cnstl1,'+')
% hold on
% plot(s,'*r')
% plot(cndt,'ok')
% xlim([-5 5]*2)
% ylim([-5 5]*2)
% hold off

