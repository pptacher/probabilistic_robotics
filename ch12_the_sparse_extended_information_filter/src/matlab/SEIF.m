% ==============================================================

%Sparse extended information filter
%written by Pierre-Paul TACHER (ptacher@gmail.com)

% =============================================================
function SEIF

fileName1 = 'aa3_lsr2.mat' ;
fileName2 = 'aa3_dr.mat' ;


%file1 = 'landmarks_20000.mat' ;
%file2 = 'poses_20000.mat' ;
%file2b = 'poses_20000.mat' ;
%file3 = 'm0_20000.mat' ;
%file4 = 'xi_20000.mat' ;
%file5 = 'omega_20000.mat' ;
%file6 = 'G_20000.mat' ;


dt = 25e-3;
Mask13 = uint16(2^13 -1) ;

global G;
G=sparse(1,1);

global  Q;
Q= diag([5 0.02]);
load(fileName1) ;
L = size(LASER);L=L(1) ;
timeLsr = double(TLsr) ; clear TLsr; j=1;

load(fileName2) ;

global m xi O
m=zeros(3,1);
xi=zeros(3,1);
O=10e4*eye(3);
m0 = zeros(1,0);

N = 20;%max active landmarks

% figure(1) ;clf ;
% zoom on ;
%hold on;
%hhh2=plot(G) ;
%axis([-50,50,0,75]);
%hold off ;

% figure(1) ;clf ;
% zoom on ;
% hhh =plot(0,0,'.','erasemode','xor') ;   %laser
% hold on;
% hhh2=plot(0,0,'+','erasemode','xor','MarkerSize',2) ;  % landmarks centers
% hhh3=plot(0,0,'o','erasemode','xor','MarkerSize',2) ;   % approx. landm. circles
% axis([-200,200,-50,350]);
% hold off ;

% G=importdata(file6);
% O=importdata(file5);
% xi=importdata(file4);
% m=importdata(file1);
% m0=importdata(file3);
% poses=importdata(file2);
% poses1=importdata(file2b);
% set(hhh2,'XData',[poses(1,:) poses1(1,:)],'YData',[poses(2,:) poses1(2,:)]) ;
% clearvars poses poses1;

poses=zeros(3,5000);
pause(2);

tstart1=tic;
stindex=20000;

for i=stindex:size(time);
    
    j=find(timeLsr>=time(i),1);
    tstart2=tic;
    
    if i>stindex && mod(i,5000)==0
        save(['poses_' num2str(i) '.mat'],'poses');
        clearvars poses;
        poses=zeros(3,5000);
        
        save(['landmarks_' num2str(i) '.mat'],'m');
        save(['xi_' num2str(i) '.mat'],'xi');
        save(['omega_' num2str(i) '.mat'],'O');
        save(['m0_' num2str(i) '.mat'],'m0');
        save(['G_' num2str(i) '.mat'],'G');
    end
    
    SEIF_motion(speed(i),steering(i),dt,m0);
    
    SEIF_estimate(m0);
    
    while j <= L && timeLsr(j) < time(i)+dt*1000
        m3 = zeros(1,0);c=zeros(1,0);
        raw = double(  bitand( Mask13,LASER(j,:)) ) ;
        raw = raw/100;
        z = detectTreesI16(raw);
        if ~isempty(z)
            z = z(1:2,:);
            c = SEIF_correspondence(z,m0);
            SEIF_measurement(z,c);
        end
        j=j+1;
        m3 = union(m3,c);
        
        if ~isempty(m3)
            n=(size(xi,1)-3)/2;
            
            m0a = setdiff(m0,m3);
            q = [m0a unique(m3)];
            
            tmp=q(1,max(1,size(q,2)-N+1):end);
            
            m1 = setdiff(m0,tmp);
            m1 = union(m1,setdiff(m3,tmp));
            m0=tmp;
            m2 = setdiff(setdiff(1:n,m1),m0);
            if ~isempty(m1)
                SEIF_sparsification(m0,m1,m2);
            end
        end
    end
    
    titer1=toc(tstart1);titer2=toc(tstart2);clc;
    fprintf('iter: %d\n', i);
    fprintf('iter time: %.5f\n',titer2 );
    fprintf('avg: %.5f\n', titer1/(i-stindex+1));
    poses(:,1+mod(i,5000))=m(1:3);
    
    %         if mod(i,200)==0
    %             set(hhh3,'XData', m(2:2:end)','YData',m(3:2:end)') ;
    %         end
    %         title( num2str(i) );
    %         set(hhh2,'XData',[get(hhh2,'XData') m(1)],'YData',[get(hhh2,'YData') m(2)]) ;
    %         drawnow ;
    
end
end
