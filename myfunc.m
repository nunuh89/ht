function dy = myfunc(t,y)
Sc=.55;C25=2;C26=2;
dy=zeros(5,1);
dy(1)=y(2);
dy(2)=y(3);
dy(3)=-y(1)*y(3)/2;
dy(4)=Sc*y(1)/C26;
dy(5)=-exp(-y(4))/C25;
end

