clear all;

%By Xiaoquan Sun
%Build output files.
%trajectory.xyz is the trajectory file without PBC(Periodic Boundary Condtion).
%trajectory_imaged.xyz is the trajectory file with PBC.
%md.out is the output file contains information of Step, Temperature, Energies(Kinetic, Potential, Total).
%msd_r.out is MSD(Mean Square Displacement) for diffusion coefficient calculation.
s=warning;
warning off all;
delete('trajectory.xyz','trajectory_imaged.xyz','md.out','msd_r.out');
warning(s);
file_trajectory=fopen('trajectory.xyz','a');
file_trajectory_imaged=fopen('trajectory_imaged.xyz','a');
file_energy=fopen('md.out','a');
file_msd=fopen('msd_r.out','a');

%Parameters.
kb=1;%Boltzmann constant.
natom=100;%Number of atoms.
boxsize=14;%The side-length of the cubic box, unit angstrom.
mass=1;%Mass of one atom.
dt=1e-3;%Time step. Unit second.
temp=300;%Temperature, unit K.
step=100000;%Total running step.
outstep=50;%Output frequency.
epsilon=1;%Parameter for LJ potential, the depth of the potential well.
sigma=2.5;%Parameter for LJ potential, the distance at which the inter-particle potential is 0.
rc=2.5*sigma;%Cutoff distance, the LJ potential is truncated at the cutoff distance.
ecut=4*epsilon*((sigma/rc)^12 - (sigma/rc)^6);%LJ potential at the cutoff distance.

tic; %Timer
%Initilization position.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=floor(power(natom,1/3))+1;
for i=1:power(n,3)
    p(1,i)=mod(i,n);
    if (mod(i,n)==0)
        p(1,i)=n;
    end
end
    j=1;
for i=1:n:power(n,3)
    j=mod(j,n);
    if (mod(j,n)==0)
        j=n;
    end
    p(2,i:i+(n-1))=j;
    j=j+1;
end
    j=1;
for i=1:power(n,2):power(n,3)
    p(3,i:i+(power(n,2)-1))=j;
    j=j+1;
end
p=p*boxsize/n;%Adjust the lattice to the box.
p=p(:,randperm(power(n,3),natom));%p is the position for all atoms.
p_ref=p; %Choose initial position as reference for diffusion coefficient calculation.

%Calculate potential and force for the initial condition.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pen=0;%pen is the potential energy.
f=zeros(3,natom);%f is the force.
for i=1:(natom-1)
    for j=(i+1):natom
        r1=p(:,i)-p(:,j);
        r1=r1-boxsize*round(r1./boxsize); %Consider PBC to get the correct distance between two atoms.
        r=sqrt(sum(r1.^2));
        if (r <= rc)
            ff=24*epsilon*sigma^6*(2*sigma^6-r^6)/(r^14);
            f(:,i)=f(:,i)+ff*r1;
            f(:,j)=f(:,j)-ff*r1;
            pen=pen+4*epsilon*((sigma/r)^12-(sigma/r)^6)-ecut;
        end
    end
end
acce=f/mass;%Calculate acceleration on each atom.

%Initilization velocity.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=rand(3,natom); %Generate random numbers as initial velocities.
ave_v=sum(v,2)/natom; %Calculate the average value of velocities.
v=(v-repmat(ave_v,[1,natom])); %Scale velocities to make sure the total momentum is zero.
ave_v=sum(v,2)/natom;
ave_v2=sum(sum(v.^2,2)/natom);
fa=sqrt(3*kb*temp/(mass*ave_v2)); %Velocity scale factor to specific temperature.
v=(v-repmat(ave_v,[1,natom]))*fa; %Scale velocites to make sure the total momentum is zero and fit the temperature.
ken=sum(sum(0.5*mass*v.^2,2)); %ken is the kinetic energy.

%Integration.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:step
    p=p+v*dt+0.5*acce*dt^2; %Update position.
 
    %Calculate force and potential according to the new position.
    pen=0;
    f=zeros(3,natom);
    for i=1:(natom-1)
        for j=(i+1):natom
            r1=p(:,i)-p(:,j);
            r1=r1-boxsize*round(r1./boxsize);
            r=sqrt(sum(r1.^2));
            if (r <= rc)
                ff=24*epsilon*sigma^6*(2*sigma^6-r^6)/(r^14);
                f(:,i)=f(:,i)+ff*r1;
                f(:,j)=f(:,j)-ff*r1;
                pen=pen+4*epsilon*((sigma./r).^12-(sigma./r).^6)-ecut;
            end
        end
    end

    v=v+0.5*acce*dt+0.5*f/mass*dt^2;%Update velocity.
    ken=sum(sum(0.5*mass*v.^2,2));%Calculate kinetic energy.
    realt=mass*sum(sum(v.^2,2)/natom)/(3*kb);%Calculate real temperature.
    
    %Output according to the frequency setting.
    if (mod(m,outstep)==0) %When the running step fit the output frequency, output the following information.
        %Output MSD result.
        msd_r=sum(sum((p-p_ref).^2,2))/natom; %calculte MSD
        fprintf(file_msd,[num2str(m),'\t',num2str(msd_r,'%5.3e'),'\n']);
        
        %Output trajectory file without PBC.
        fprintf(file_trajectory,[num2str(natom),'\n']);
        fprintf(file_trajectory,[num2str(m),' step \n']);
        for i = 1:natom
            fprintf(file_trajectory,['C','\t']);
            for j =1:3
                fprintf(file_trajectory,[num2str(p(j,i)'),'\t']);
            end
            fprintf(file_trajectory,'\n');
        end
        
        %Output trajectory file with PBC.
        pi=p-boxsize*round(p./boxsize); %pi is the postion with PBC.
        fprintf(file_trajectory_imaged,[num2str(natom),'\n']);
        fprintf(file_trajectory_imaged,[num2str(m),' step \n']);
        for i = 1:natom
            fprintf(file_trajectory_imaged,['C','\t']);
            for j =1:3
                fprintf(file_trajectory_imaged,[num2str(pi(j,i)'),'\t']);
            end
            fprintf(file_trajectory_imaged,'\n');
        end
        
        %Output step and energies.
        fprintf(file_energy,[num2str(m),'\t','T= ',num2str(realt,'%5.2f'),'\t']);
        fprintf(file_energy,['Kinetic= ',num2str(ken,'%5.3e'),'\t']);
        fprintf(file_energy,['Potential= ',num2str(pen,'%5.3e'),'\t']);
        fprintf(file_energy,['Total= ',num2str((ken+pen),'%5.3e'),'\n']);
    end
end
fprintf(file_energy,['#', num2str(step),' steps MD simulations cost ',num2str(toc),' seconds. \n']); %Output the total running time
fprintf(file_energy,['#', 'The speed is about ',num2str(step/toc,'%8.3f'),' step/second. \n']); %Output the calculation speed.
