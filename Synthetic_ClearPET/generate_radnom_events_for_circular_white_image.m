function [x,y]=generate_radnom_events_for_circular_white_image(N,R)
NOE = round(4/pi*N*1.2); %number of events
x =  2*R*rand(NOE,1)-R;
y =  2*R*rand(NOE,1)-R;

ind = x.^2+y.^2>=R^2;
x(ind)=[];
y(ind)=[];
x(round(N)+1:end)=[];
y(round(N)+1:end)=[];
end