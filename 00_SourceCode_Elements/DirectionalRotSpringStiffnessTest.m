

direction = 1;
MVfactor = 100;
sprRotK = 1;

tempK=zeros(100,1);
theta=zeros(100,1);

for i=1:1000

    theta(i)=i/1000*2*pi;

        if theta(i)<pi && direction==0
            tempK(i) = sprRotK;
        %elseif theta(i)<pi+0.1*pi && direction==0
            %tempK(i) = sprRotK+sprRotK*((theta(i)-pi)/(0.1*pi))*(MVfactor-1);

        elseif theta(i)>pi && direction==1
            tempK(i) = sprRotK;
        %elseif theta(i)>pi-0.1*pi && direction==1
            %tempK(i) = sprRotK+sprRotK*((pi-theta(i))/(0.1*pi))*(MVfactor-1);
        
        else
            tempK(i) = sprRotK*MVfactor;
        end
end

figure
plot(theta/pi, tempK)


for i=1:1000

    theta(i)=i/1000*2*pi;
    if direction==0    
        factor=(tanh((theta(i)-pi)*10)+1)/2*MVfactor+1;
        tempK(i) = sprRotK*factor;
    else
        factor=(tanh(-(theta(i)-pi)*10)+1)/2*MVfactor+1;
        tempK(i) = sprRotK*factor;
    end

end

figure
plot(theta/pi, tempK)