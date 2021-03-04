function dms = rt_prttrace(dm, dims, sind)

dms = 0;
Base1 = eye(dims(1));
Base2 = eye(dims(2));
for i = 1:dims(sind)
    if sind == 1
        Vec = kron(Base1(:,i),Base2);
    else
        Vec = kron(Base1,Base2(:,i));
    end
    dms = dms + Vec'*dm*Vec;
end

end