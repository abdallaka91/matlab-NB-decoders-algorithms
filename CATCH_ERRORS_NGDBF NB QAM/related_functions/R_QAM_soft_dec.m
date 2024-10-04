function [HD_bin, soft_dec, d, cndt, HD]  = R_QAM_soft_dec(y, norm_const, gray_bin_words, I, Q)
p=log2(length(norm_const));

ry = real(y);
iy = imag(y);

[d, cndt] = square_qam_nearest(y, norm_const);
[~,i1] = min(d);
HD = cndt(i1);
HD_bin = gray_bin_words(HD==norm_const,:);

soft_dec = nan(1,p);

if ry>=I(1) && ry<=I(end) && iy>=Q(1) && iy<=Q(end)
    cx = 0;
    if ry>real(HD)
        if real(HD)~=I(end)
            x_ngbr = HD+2;
        else
            x_ngbr=HD-2;
            cx = 1;
        end
    else
        if real(HD)~=I(1)
            x_ngbr = HD-2;
        else
            cx = 2;
            x_ngbr = HD+2;
        end
    end
    x_ngbr_bin = gray_bin_words(x_ngbr==norm_const,:);
    bit_dff = find(HD_bin~=x_ngbr_bin);
    if ~cx
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = abs(ry-real(HD))-1;
        else
            soft_dec(bit_dff) = 1-abs(ry-real(HD));
        end
    elseif cx==1
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = 1-abs((ry-real(HD)));
        else
            soft_dec(bit_dff) = -(1-abs((ry-real(HD))));
        end
    end

    cy = 0;
    if iy>imag(HD)
        if imag(HD)~=Q(end)
            y_ngbr = HD+2i;
        else
            y_ngbr=HD-2i;
            cy = 1;
        end
    else
        if imag(HD)~=Q(1)
            y_ngbr = HD-2i;
        else
            cy = 2;
            y_ngbr = HD+2i;
        end
    end


    HD_bin = gray_bin_words(HD==norm_const,:);
    y_ngbr_bin = gray_bin_words(y_ngbr==norm_const,:);
    bit_dff = find(HD_bin~=y_ngbr_bin);
    if ~cy
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = abs(iy-imag(HD))-1;
        else
            soft_dec(bit_dff) = 1-abs(iy-imag(HD));
        end
    elseif cy==1
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = 1-abs((iy-imag(HD)));
        else
            soft_dec(bit_dff) = -(1-abs((iy-imag(HD))));
        end
    end
    soft_dec(isnan(soft_dec)) =  HD_bin(isnan(soft_dec)) ;


elseif ry>=I(1) && ry<=I(end)
    cx = 0;
    if ry>real(HD)
        if real(HD)~=I(end)
            x_ngbr = HD+2;
        else
            x_ngbr=HD-2;
            cx = 1;
        end
    else
        if real(HD)~=I(1)
            x_ngbr = HD-2;
        else
            cx = 2;
            x_ngbr = HD+2;
        end
    end
    x_ngbr_bin = gray_bin_words(x_ngbr==norm_const,:);
    bit_dff = find(HD_bin~=x_ngbr_bin);
    if ~cx
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = abs(ry-real(HD))-1;
        else
            soft_dec(bit_dff) = 1-abs(ry-real(HD));
        end
    elseif cx==1
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = 1-abs((ry-real(HD)));
        else
            soft_dec(bit_dff) = -(1-abs((ry-real(HD))));
        end
    end
    soft_dec(isnan(soft_dec)) =  HD_bin(isnan(soft_dec)) ;

elseif iy>=Q(1) && iy<=Q(end)
    cy = 0;
    if iy>imag(HD)
        if imag(HD)~=Q(end)
            y_ngbr = HD+2i;
        else
            y_ngbr=HD-2i;
            cy = 1;
        end
    else
        if imag(HD)~=Q(1)
            y_ngbr = HD-2i;
        else
            cy = 2;
            y_ngbr = HD+2i;
        end
    end
    y_ngbr_bin = gray_bin_words(y_ngbr==norm_const,:);
    bit_dff = find(HD_bin~=y_ngbr_bin);
    if ~cy
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = abs(iy-imag(HD))-1;
        else
            soft_dec(bit_dff) = 1-abs(iy-imag(HD));
        end
    elseif cy==1
        if HD_bin(bit_dff)==-1
            soft_dec(bit_dff) = 1-abs((iy-imag(HD)));
        else
            soft_dec(bit_dff) = -(1-abs((iy-imag(HD))));
        end
    end
else
    soft_dec(isnan(soft_dec)) =  HD_bin(isnan(soft_dec)) ;
end
soft_dec(isnan(soft_dec)) =  HD_bin(isnan(soft_dec)) ;