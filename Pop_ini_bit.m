function [pop_ini] = Pop_ini_bit(pop_ini_size,pd,n)


G_lib=(n*(n-1)*0.5);
pop_ini=zeros(pop_ini_size,G_lib);
bit_represent=zeros(G_lib,1);
for k=1:pop_ini_size
    ini=1;
    for i =1:n-1
        count=ini:ini+(n-i-1);
        bit_represent(count)= random(pd,1, n-i);
        ini=count(end)+1;
    end
    bit_represent=bit_represent-2;
    pop_ini(k,:)=bit_represent;
end

end