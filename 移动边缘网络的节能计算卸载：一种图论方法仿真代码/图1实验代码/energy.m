function y =energy()
E=zeros(10,10);
for i=1:10
    E(i,:)= 7+6.*rand(1,10);
    y=E;
end