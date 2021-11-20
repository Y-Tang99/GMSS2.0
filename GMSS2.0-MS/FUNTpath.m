function [Tpath]=FUNTpath(R)

if R<=50
    Tpath=1+0.18*R;
else
    Tpath=10+0.05*(R-50);
end

end