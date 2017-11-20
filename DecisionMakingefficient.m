%% Recurrent network of IF neurons with NMDA, AMPA, and GABA channels
% Probabilistic Decision Making by Slow Reverberation in Cortical Circuits (Wang, 2002)%
clear all; clf;

%Percent of total neuron population that is...
%All further variables with lowercase e & i indicates it to be a pyradimal
%or interneuron measurement respectively Constants
N=1000; Pyramidal=0.8; Inter=1-Pyramidal; Vl=-70; 

%Setting Connectivity, assuming fully connected network plus 800 external inputs
%could add  acessin variales i.e. Inhibitory neurons Ni: Ne+1:(Ne+Ni)
Ne=ceil(N*Pyramidal); Ni=ceil(N*Inter); 

%Setting rest potential, and synaptic current in mV & mA
VN=zeros(N,1)-55+randn(N,1); Isyn=zeros(N,1); 

%Leakage Conductances (g) in mS and receptor conductances 
g_e=25; g_i=20; Cm_e=500; Cm_i=200;
s_Aext=zeros(N,1); s_Arec=zeros(N,Ne); 
s_N=zeros(N,Ne); x=zeros(N,Ne); s_G=zeros(N,Ni);

%objects are represented by assemblie o, which should not overlap 
%Note that in this scheme all objects recruit the same number of
%neurons to represent them, different
f=0.15; p=2; o=f*Ne; 
%o=1; f=o/Ne; %use these numbers when N=10
%Weighting equations
w=ones(Ne,Ne);  
w_plus=1.7; w_minus=1-f*(w_plus-1)/(1-f); 
%Nonselective to selective weights
w(o*p+1:Ne,o*p+1:Ne)=1;
w(1:o*p,1:Ne)=w_minus;
%Selective populations
w(1:o,1:o)=w_plus; w(o+1:2*o,o+1:2*o)=w_plus; 
%B&W, set w_plus=2.1, o=2.1, f=0.1 and p=4 
%w(2*o+1:3*o,2*o+1:3*o)=w_plus; w(3*o+1:4*o,3*o+1:4*o)=w_plus;
%Setting time step index in ms
dt=0.1; %time(1:2)=0;

%Start generating Poisson spike train at 2.4 KHz, setting stimulus intensity
lambda=zeros(N,1)+2400/(1000/dt); stimA=64; stimB=36; co=80; mnot=40; 
meu_A=mnot+(mnot/100)*co; meu_B= mnot-(mnot/100)*co; SD=10; dt=0.1;



%buffer to account for latency
i=1; syn_latency=round(0.5/dt); trefe=round(2/dt); trefi=round(1/dt);
firings=zeros(N,syn_latency); fired=[];
refractory=zeros(N,1);

for t=0:dt:1500;
    
    i=i+1
         
    %% Integrate-and-fire
    %Ohm's Law:
    VN(1:Ne)=VN(1:Ne)+dt*(-g_e/Cm_e*(VN(1:Ne)+70)-Isyn(1:Ne,1)/Cm_e);
    VN(Ne+1:N)=VN(Ne+1:N)+dt*(-g_i/Cm_i*(VN(Ne+1:N)+70)-Isyn(Ne+1:N)/Cm_i);
    
    %Stimulus Protocol
    %Sampling mean & determining stimulus occurence
    if ceil(t/50)==floor(t/50)&&t>=200&&t<=1200;
        stimA=SD*randn+meu_A;
    else stimA=stimA;
    end 
    if t>=1200; stimA=0; end

    if ceil(t/50)==floor(t/50)&&t>=200&&t<=1200;
        stimB=SD*randn+meu_B;
    else stimB=stimB;
    end 
    if t>=1200; stimB=0; end
   
    %Updating lambda
      lambda(1:o)=(2400+stimA)/(1000/dt);
      lambda(o+1:o*2)=(2400+stimB)/(1000/dt);
    
    %%External Input
   delta_ext=zeros(N,1);
   
   %vector of random numbers
   chance=rand(N,1);
   n=1;
   while 0<sum(chance<=(lambda.^n).*exp(-lambda)/factorial(n));
       delta_ext(find(chance<=(lambda.^n).*exp(-lambda)/factorial(n)))=n;
       n=n+1;
   end
    
    %Refractory & Reset, first checks if theres been a firing then 
    %goes back through the firings history tref timesteps back to
    %check if there has been a firing if there has, VN is held at reset
    refractory=refractory-1; refractory(find(refractory<0))=0;
    VN(find(refractory>0))=-55;
    refractory(find(VN(1:Ne)>-50))=trefe; refractory(find(VN(Ne+1:N)>-50)+Ne)=trefi;
     
    
    spikes=find(VN>-50);
    fired=[fired; t+0*spikes, spikes];
    firings=[firings,VN>-50];  firings=firings(:,2:syn_latency+1);
    
    %inputs
    delta_ext=delta_ext/dt;
    delta_e=firings(1:Ne,1)/dt;
    delta_i=firings(Ne+1:N,1)/dt;
    
    
    %% Update Synaptic conductance matrices and currents
    s_Aext=s_Aext+dt*(-0.5*s_Aext+delta_ext);
    s_Arec=s_Arec+dt*(-0.5*s_Arec+repmat(delta_e',N,1));
    x=x+dt*(-0.5*x+repmat(delta_e',N,1));
    s_N=s_N+dt*(-0.01*s_N+0.5*x.*(1-s_N));
    s_G=s_G+dt*(-0.2*s_G+repmat(delta_i',N,1));
    
    %Currents  For Wang 2002    
     Isyn(1:Ne)= (2.1)*VN(1:Ne).*s_Aext(1:Ne) ...
                +(0.05)*VN(1:Ne).*sum(w.*s_Arec(1:Ne,:),2) ...
                +(0.165)*VN(1:Ne).*sum(w.*s_N(1:Ne,:),2) ./ (1+1*exp(-0.062*VN(1:Ne)/3.57)) ...
                +(1.3)*(VN(1:Ne)+70).*sum(s_G(1:Ne,:),2);
            
            
     Isyn(Ne+1:N)= (1.62)* VN(Ne+1:N).*s_Aext(Ne+1:N) ...
                    +(0.04)* VN(Ne+1:N).*sum(s_Arec(Ne+1:N,:),2) ...
                    +(0.13)*VN(Ne+1:N).*sum(s_N(Ne+1:N,:),2) ./ ((1+1*exp((-0.062*VN(Ne+1:N))/3.57))) ...
                    +(1.0)*(VN(Ne+1:N)+70).*sum(s_G(Ne+1:N,:),2);
                
      %Currents  For Brunell & Wang 2001           
%      Isyn(1:Ne)= (2.08)*VN(1:Ne).*s_Aext(1:Ne) ...
%                +(0.104)*VN(1:Ne).*sum(w.*s_Arec(1:Ne,:),2) ...
%                +(0.327)*VN(1:Ne).*sum(w.*s_N(1:Ne,:),2) ./ (1+1*exp(-0.062*VN(1:Ne)/3.57)) ...
%                +(1.25)*(VN(1:Ne)+70).*sum(s_G(1:Ne,:),2);
%
%            
%     Isyn(Ne+1:N)= (1.62)* VN(Ne+1:N).*s_Aext(Ne+1:N) ...
%                    +(0.081)* VN(Ne+1:N).*sum(s_Arec(Ne+1:N,:),2) ...
%                    +(0.258)*VN(Ne+1:N).*sum(s_N(Ne+1:N,:),2) ./ ((1+1*exp((-0.062*VN(Ne+1:N))/3.57))) ...
%                    +(0.973)*(VN(Ne+1:N)+70).*sum(s_G(Ne+1:N,:),2);           
  

end

%Firing rates, note time must always be multiple of 5 or this will give an
%error
%Calculating firing rates
timestep=5; window=50; c=0;
Fr_A=[]; Fr_B=[];Fr_I=[]; Fr_NS=[];

fA=fired(find(fired(:,2)<=o)); 
fB=fired(find(fired(:,2)<=o*2));%fB=fB(find(fired(find(fired(:,2)<=o))):end); 
fI=fired(find(fired(:,2)<=N));%fI=fI(find(fired(find(fired(:,2)<=Ne+1))):end); 
fNS=fired(find(fired(:,2)<=Ne));%fNS=fNS(find(fired(find(fired(:,2)<=o*2+1))):end);
    for step=0:timestep:t;
    
   sumA=sum(fA<(window/2+timestep*c))-sum(fA<(timestep*c-window/2));
   instA=sumA/(o*(window/1000)); Fr_A=[Fr_A; timestep*c, instA];
   sumB=sum(fB<(window/2+timestep*c))-sum(fB<(timestep*c-window/2))-sumA; if sumB<0; sumB=0; end;
   instB=sumB/(o*(window/1000)); Fr_B=[Fr_B; timestep*c, instB]; 
   sumNS=sum(fNS<(window/2+timestep*c))-sum(fNS<(timestep*c-window/2))-sumA-sumB; if sumNS<0; sumNS=0; end;
   instNS=sumNS/((Ne-o*2)*(window/1000)); Fr_NS=[Fr_NS; timestep*c, instNS];
   sumI=sum(fI<(window/2+timestep*c))-sum(fI<(timestep*c-window/2))-sumA-sumB-sumNS; if sumI<0; sumI=0; end;
   instI=sumI/(Ni*(window/1000)); Fr_I=[Fr_I; timestep*c, instI];
   c=c+1;
end

subplot(2,2,1:2);
plot(fired(:,1),fired(:,2),'.','Markersize',2); xlabel('Time (ms)'); title('Raster Plot')
subplot(2,2,3); plot(Fr_A(:,1),Fr_A(:,2),'r-'); hold on
plot(Fr_B(:,1),Fr_B(:,2),'b-'); plot(Fr_NS(:,1),Fr_NS(:,2),'g-'); plot(Fr_I(:,1),Fr_I(:,2),'k-');
xlabel('Time (ms)'); ylabel ('Firaing rate (Hz)'); title('Population Average Firing rates (Hz)');
legend('A Selective','B Selecetive','Non-Selective','Inhibitory','Location','NorthWest')
subplot(2,2,4); plot(Fr_B(:,2),Fr_A(:,2),'r-');  hold on
plot(1:(t*0+50),1:(i*0+50),'k-.'); xlabel('Firaing rate B (Hz)'); ylabel ('Firaing rate A (Hz)'); title('Decision Space')

%For plotting input:
%subplot(2,2,4); plot(trackA(:,1),trackA(:,2),'r-'); hold on 
%plot(trackB(:,1),trackB(:,2), 'b-');
%xlabel('Time (ms)'); ylabel ('Mean Firing Rate (Hz)'); title('Stimulus Input');legend('Evidence for A','Evidence for B');

